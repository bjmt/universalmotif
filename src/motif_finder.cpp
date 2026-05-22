// motif_finder.cpp -- yamtk-aligned _de novo_ motif discovery.
//
// Faithful port of /Users/gok24jef/yamtk/src/yamme.c. Pipeline:
//   1. For each motif width w in [min_w, max_w], discover up to nmotifs:
//      a. Enumerate candidate seeds (k-mer Fisher's exact log-p, per-seq
//         canonical-strand presence).
//      b. Try each seed: align Hamming-1 neighbours to build initial PPM.
//      c. Convert PPM -> integer PWM, build score CDF, choose threshold
//         by hit p-value.
//      d. Refine REFINE_PASSES times (re-scan, rebuild PPM, re-threshold).
//      e. Evaluate on positives + negatives; accept if Fisher's per-seq
//         presence p-value below stop.pvalue.
//      f. On accept, mask covered positions in the per-thread positives
//         buffer and move to next iteration.
//   2. Merge per-thread results in deterministic (width, discovery_seq)
//      order. Re-evaluate every accepted motif on the pristine
//      (unmasked) positives so coverage bitmasks are comparable.
//   3. Cross-width dedup: drop higher-pvalue motifs whose coverage
//      overlap with a lower-pvalue motif exceeds dedup_overlap.
//   4. Trim low-IC flanks. Return per-motif PPMs + stats.
//
// Threading: per-width parallelism via RcppThread::parallelFor, matching
// yamtk's worker pool. State per thread (positives copy, CDF scratch,
// local result list) is fully isolated. No shared mutable state inside
// the parallel region.

#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <climits>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// -----------------------------------------------------------------------------
// Constants (mirror yamme.c)
// -----------------------------------------------------------------------------
static constexpr int      HAMMING_MISMATCH    = 1;
static constexpr int      REFINE_PASSES       = 2;
static constexpr int      TOP_K_SEEDS         = 4;
static constexpr int      MAX_MOTIF_WIDTH     = 30;
static constexpr int      MIN_REFINE_HITS     = 20;
static constexpr double   MIN_IC_BITS         = 0.5;
static constexpr double   CONSENSUS_THRESHOLD = 0.75;
static constexpr double   PWM_INT_MULTIPLIER  = 1000.0;
static constexpr int      AMBIGUITY_SCORE     = -10000000;
static constexpr uint64_t MAX_CDF_SIZE        = 2097152;
static constexpr int      COMP4[4]            = {3, 2, 1, 0};

// Lookup: char -> {A=0, C=1, G=2, T/U=3, else=4}.  Built at startup.
static unsigned char CHAR2IDX[256];
static bool          CHAR2IDX_INIT = false;
static void init_char2idx() {
  if (CHAR2IDX_INIT) return;
  for (int i = 0; i < 256; i++) CHAR2IDX[i] = 4;
  CHAR2IDX[(unsigned char) 'A'] = 0; CHAR2IDX[(unsigned char) 'a'] = 0;
  CHAR2IDX[(unsigned char) 'C'] = 1; CHAR2IDX[(unsigned char) 'c'] = 1;
  CHAR2IDX[(unsigned char) 'G'] = 2; CHAR2IDX[(unsigned char) 'g'] = 2;
  CHAR2IDX[(unsigned char) 'T'] = 3; CHAR2IDX[(unsigned char) 't'] = 3;
  CHAR2IDX[(unsigned char) 'U'] = 3; CHAR2IDX[(unsigned char) 'u'] = 3;
  CHAR2IDX_INIT = true;
}

// -----------------------------------------------------------------------------
// Motif representation (port of motif_t).  PWMs are stored row-major,
// flattened as pwm[pos*5 + idx], where idx=4 (N) holds AMBIGUITY_SCORE.
// -----------------------------------------------------------------------------
struct Motif {
  std::vector<int>    pwm;            // size = w*5
  std::vector<int>    pwm_rc;         // size = w*5
  std::vector<std::array<double, 4>> pwm_probs;  // size = w  (pseudocount-adj)
  std::vector<double> cdf;            // null score CDF (tail prob)
  int                 threshold;      // accept-hit score cutoff
  int                 mn;             // min per-position PWM entry
  int                 mx;             // max per-position PWM entry
  int                 max_score;      // sum of per-position max
  int                 min_score;      // sum of per-position min
  int                 cdf_max;        // mx - mn
  int                 cdf_offset;     // mn * w  (CDF index offset)
  uint64_t            cdf_size;       // w * cdf_max + 1
  uint64_t            size;           // motif width
  uint64_t            nsites_actual;  // refined hit count

  void init(uint64_t w) {
    size = w;
    pwm.assign(w * 5, 0);
    pwm_rc.assign(w * 5, 0);
    for (uint64_t i = 0; i < w; i++) {
      pwm[i * 5 + 4]    = AMBIGUITY_SCORE;
      pwm_rc[i * 5 + 4] = AMBIGUITY_SCORE;
    }
    pwm_probs.assign(w, {0.0, 0.0, 0.0, 0.0});
    cdf.clear();
    threshold = 0;
    mn = mx = 0; max_score = 0; min_score = 0;
    cdf_max = 0; cdf_offset = 0; cdf_size = 0;
    nsites_actual = 0;
  }
};

static inline void set_score(Motif &m, int idx, uint64_t pos, int s) {
  m.pwm[idx + pos * 5] = s;
}
static inline void set_score_rc(Motif &m, int idx, uint64_t pos, int s) {
  m.pwm_rc[idx + pos * 5] = s;
}
static inline int get_score(const Motif &m, unsigned char ch, uint64_t pos) {
  return m.pwm[CHAR2IDX[ch] + pos * 5];
}
static inline int get_score_rc(const Motif &m, unsigned char ch, uint64_t pos) {
  return m.pwm_rc[CHAR2IDX[ch] + pos * 5];
}
static inline int get_score_i(const Motif &m, int idx, uint64_t pos) {
  return m.pwm[idx + pos * 5];
}

// Parameterised log-odds score (no global state).
static int calc_score_ns(double prob_i, double bkg_i, int ns, int pc) {
  double x = prob_i * ns;
  x += pc / 4.0;
  x /= (ns + pc);
  if (x <= 0.0) x = 1e-12;
  return (int) (std::log2(x / bkg_i) * PWM_INT_MULTIPLIER);
}

static void fill_pwm_rc(Motif &m) {
  // pwm_rc[base, pos] = pwm[comp(base), size-1-pos]
  for (uint64_t pos = 0; pos < m.size; pos++) {
    for (int idx = 0; idx < 4; idx++) {
      set_score_rc(m, idx, m.size - 1 - pos,
                   get_score_i(m, COMP4[idx], pos));
    }
  }
}

static int get_pwm_min(const Motif &m) {
  int mn = 0;
  for (uint64_t pos = 0; pos < m.size; pos++)
    for (int let = 0; let < 4; let++) {
      int v = get_score_i(m, let, pos);
      if (v < mn) mn = v;
    }
  return mn;
}
static int get_pwm_max(const Motif &m) {
  int mx = 0;
  for (uint64_t pos = 0; pos < m.size; pos++)
    for (int let = 0; let < 4; let++) {
      int v = get_score_i(m, let, pos);
      if (v > mx) mx = v;
    }
  return mx;
}

// -----------------------------------------------------------------------------
// Score CDF (port of fill_cdf + set_threshold_pval).
// -----------------------------------------------------------------------------
static int fill_cdf(Motif &m, const std::array<double, 4> &bkg) {
  if (m.cdf_size == 0 || m.cdf_size > MAX_CDF_SIZE) return 1;
  m.cdf.assign(m.cdf_size, 1.0);
  std::vector<double> tmp_pdf(m.cdf_size);
  for (uint64_t i = 0; i < m.size; i++) {
    uint64_t max_step = i * (uint64_t) m.cdf_max;
    std::copy(m.cdf.begin(), m.cdf.end(), tmp_pdf.begin());
    std::fill(m.cdf.begin(),
              m.cdf.begin() + std::min<uint64_t>(max_step + m.cdf_max + 1, m.cdf_size),
              0.0);
    for (int j = 0; j < 4; j++) {
      uint64_t s = (uint64_t) (get_score_i(m, j, i) - m.mn);
      for (uint64_t k = 0; k <= max_step; k++)
        m.cdf[k + s] += tmp_pdf[k] * bkg[j];
    }
  }
  double pdf_sum = 0.0;
  for (uint64_t i = 0; i < m.cdf_size; i++) pdf_sum += m.cdf[i];
  if (std::fabs(pdf_sum - 1.0) > 0.0001) {
    for (uint64_t i = 0; i < m.cdf_size; i++) m.cdf[i] /= pdf_sum;
  }
  // Convert PDF to tail probability (CDF from the top).
  for (uint64_t i = m.cdf_size - 2; i < (uint64_t) -1; i--)
    m.cdf[i] += m.cdf[i + 1];
  return 0;
}

static void set_threshold_pval(Motif &m, double pval) {
  uint64_t threshold_i = m.cdf_size;
  for (uint64_t i = 0; i < m.cdf_size; i++) {
    if (m.cdf[i] < pval) { threshold_i = i; break; }
  }
  m.threshold -= m.mn;
  m.threshold *= (int) m.size;
  m.threshold = (int) threshold_i - m.threshold;
  m.max_score = 0; m.min_score = 0;
  for (uint64_t i = 0; i < m.size; i++) {
    int mx = get_score_i(m, 0, i), mn = mx;
    for (int j = 1; j < 4; j++) {
      int t = get_score_i(m, j, i);
      if (t > mx) mx = t;
      if (t < mn) mn = t;
    }
    m.max_score += mx;
    m.min_score += mn;
  }
  double min_pval = m.cdf[m.max_score - m.cdf_offset];
  if (min_pval / pval > 1.0001) m.threshold = INT_MAX;
}

// -----------------------------------------------------------------------------
// PPM -> motif (port of convert_ppm_to_motif).
// -----------------------------------------------------------------------------
static void convert_ppm_to_motif(Motif &m,
                                 const int ppm[][4],
                                 uint64_t w, int nsites,
                                 const std::array<double, 4> &bkg,
                                 int pseudocount, double hit_pval) {
  m.init(w);
  m.nsites_actual = (uint64_t) nsites;

  for (uint64_t j = 0; j < w; j++) {
    double col = (double) (ppm[j][0] + ppm[j][1] + ppm[j][2] + ppm[j][3]);
    if (col <= 0.0) col = 1.0;
    for (int i = 0; i < 4; i++) {
      double raw = ppm[j][i] / col;
      m.pwm_probs[j][i] = (raw * nsites + pseudocount / 4.0) /
                          (nsites + pseudocount);
    }
    set_score(m, 0, j, calc_score_ns(ppm[j][0] / col, bkg[0], nsites, pseudocount));
    set_score(m, 1, j, calc_score_ns(ppm[j][1] / col, bkg[1], nsites, pseudocount));
    set_score(m, 2, j, calc_score_ns(ppm[j][2] / col, bkg[2], nsites, pseudocount));
    set_score(m, 3, j, calc_score_ns(ppm[j][3] / col, bkg[3], nsites, pseudocount));
  }

  m.mn         = get_pwm_min(m);
  m.mx         = get_pwm_max(m);
  m.cdf_offset = m.mn * (int) w;
  fill_pwm_rc(m);
  m.cdf_max    = m.mx - m.mn;
  m.cdf_size   = (uint64_t) w * (uint64_t) m.cdf_max + 1;

  if (m.cdf_max == 0 || m.cdf_size > MAX_CDF_SIZE) {
    m.threshold = INT_MAX;
    return;
  }

  if (fill_cdf(m, bkg) != 0) {
    m.threshold = INT_MAX;
    return;
  }

  m.threshold = 0; m.max_score = 0; m.min_score = 0;
  set_threshold_pval(m, hit_pval);

  if (m.threshold == INT_MAX) {
    // Fallback: relax to top 1%
    m.threshold = 0; m.max_score = 0; m.min_score = 0;
    set_threshold_pval(m, 0.01);
  }
}

// -----------------------------------------------------------------------------
// K-mer encoding (2-bit, 4-base). Returns UINT64_MAX on N.
// -----------------------------------------------------------------------------
static inline uint64_t kmer_encode4(const std::string &seq, uint64_t w, uint64_t off) {
  uint64_t kmer = 0;
  for (uint64_t i = 0; i < w; i++) {
    uint8_t idx = CHAR2IDX[(unsigned char) seq[off + i]];
    if (idx >= 4) return UINT64_MAX;
    kmer = kmer * 4 + idx;
  }
  return kmer;
}
static uint64_t rc_kmer4(uint64_t kmer, uint64_t w) {
  uint64_t rc = 0;
  for (uint64_t i = 0; i < w; i++) {
    rc = rc * 4 + (uint64_t) COMP4[kmer & 3];
    kmer >>= 2;
  }
  return rc;
}
static uint64_t hamming4(uint64_t a, uint64_t b, uint64_t w) {
  uint64_t mm = 0;
  for (uint64_t i = 0; i < w; i++) {
    if ((a & 3) != (b & 3)) mm++;
    a >>= 2; b >>= 2;
  }
  return mm;
}

// -----------------------------------------------------------------------------
// Stats (port of log_choose, fishers_exact_log_greater, bh_qvalues).
// -----------------------------------------------------------------------------
static double log_choose(double n, double k) {
  return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}
static double fishers_exact_log_greater(uint64_t a, uint64_t b,
                                        uint64_t c, uint64_t d) {
  double r1 = (double) (a + b), r2 = (double) (c + d);
  double c1 = (double) (a + c), N = (double) (a + b + c + d);
  double log_denom = log_choose(N, c1);
  uint64_t x_max = (uint64_t) std::min(r1, c1);
  double log_sum = log_choose(r1, a) + log_choose(r2, c1 - a) - log_denom;
  for (uint64_t x = a + 1; x <= x_max; x++) {
    double term = log_choose(r1, x) + log_choose(r2, c1 - x) - log_denom;
    if (term > log_sum) log_sum = term + std::log1p(std::exp(log_sum - term));
    else                log_sum = log_sum + std::log1p(std::exp(term - log_sum));
  }
  double p = std::exp(log_sum);
  return p > 1.0 ? 1.0 : p;
}
static void bh_qvalues(const std::vector<double> &p, std::vector<double> &q) {
  uint64_t n = p.size();
  if (n == 0) { q.clear(); return; }
  q.assign(n, 1.0);
  std::vector<std::pair<double, uint64_t>> srt(n);
  for (uint64_t i = 0; i < n; i++) srt[i] = {p[i], i};
  std::sort(srt.begin(), srt.end(),
            [](const std::pair<double, uint64_t> &a,
               const std::pair<double, uint64_t> &b) {
              return a.first < b.first;
            });
  double prev = 1.0;
  for (uint64_t i = n; i > 0;) {
    i--;
    double qv = srt[i].first * (double) n / (double) (i + 1);
    if (qv > prev) qv = prev;
    prev = qv;
    q[srt[i].second] = qv > 1.0 ? 1.0 : qv;
  }
}

// -----------------------------------------------------------------------------
// Scoring helpers (port of score_subseq / score_subseq_rc).
// -----------------------------------------------------------------------------
static inline int score_subseq(const Motif &m, const std::string &seq, uint64_t off) {
  int s = 0;
  for (uint64_t i = 0; i < m.size; i++)
    s += get_score(m, (unsigned char) seq[i + off], i);
  return s;
}
static inline void score_subseq_rc(const Motif &m, const std::string &seq,
                                    uint64_t off, int *s_out, int *src_out) {
  int s = 0, src = 0;
  for (uint64_t i = 0; i < m.size; i++) {
    s   += get_score(m,    (unsigned char) seq[i + off], i);
    src += get_score_rc(m, (unsigned char) seq[i + off], i);
  }
  *s_out = s; *src_out = src;
}

// -----------------------------------------------------------------------------
// Seed enumeration (port of enumerate_seeds).
// -----------------------------------------------------------------------------
struct Seed { uint64_t kmer; double pval; };

static int enumerate_seeds(const std::vector<std::string> &pos_seqs,
                           const std::vector<std::string> &neg_seqs,
                           uint64_t w, bool RC,
                           std::vector<Seed> &out) {
  out.clear();
  std::unordered_map<uint64_t, uint64_t> pos_h, neg_h;

  auto count_set = [&](const std::vector<std::string> &seqs,
                       std::unordered_map<uint64_t, uint64_t> &h) {
    for (const auto &seq : seqs) {
      uint64_t seqlen = seq.size();
      if (seqlen < w) continue;
      std::unordered_set<uint64_t> seen;
      seen.reserve(seqlen);
      for (uint64_t off = 0; off + w <= seqlen; off++) {
        uint64_t kmer = kmer_encode4(seq, w, off);
        if (kmer == UINT64_MAX) continue;
        uint64_t canon = kmer;
        if (RC) {
          uint64_t rc = rc_kmer4(kmer, w);
          canon = kmer <= rc ? kmer : rc;
        }
        seen.insert(canon);
      }
      for (uint64_t k : seen) h[k]++;
    }
  };

  count_set(pos_seqs, pos_h);
  count_set(neg_seqs, neg_h);

  uint64_t n_pos = pos_seqs.size(), n_neg = neg_seqs.size();
  if (pos_h.empty()) return 0;

  std::vector<Seed> cands;
  cands.reserve(pos_h.size());
  for (const auto &kv : pos_h) {
    uint64_t kmer = kv.first;
    uint64_t pc   = kv.second;
    uint64_t nc   = 0;
    auto it = neg_h.find(kmer);
    if (it != neg_h.end()) nc = it->second;
    Seed s;
    s.kmer = kmer;
    s.pval = fishers_exact_log_greater(pc, nc,
                                       n_pos > pc ? n_pos - pc : 0,
                                       n_neg > nc ? n_neg - nc : 0);
    cands.push_back(s);
  }
  std::sort(cands.begin(), cands.end(),
            [](const Seed &a, const Seed &b) { return a.pval < b.pval; });
  if (cands.size() > (uint64_t) TOP_K_SEEDS) cands.resize(TOP_K_SEEDS);
  out = std::move(cands);
  return 0;
}

// -----------------------------------------------------------------------------
// Build PPM from seed via Hamming alignment (port of build_ppm_from_seed).
// `ppm` is int[MAX_MOTIF_WIDTH][4], zero-initialised by caller.
// -----------------------------------------------------------------------------
static int build_ppm_from_seed(const std::vector<std::string> &pos_seqs,
                               uint64_t seed_kmer, uint64_t w, bool RC,
                               int ppm[][4]) {
  uint64_t rc_seed = RC ? rc_kmer4(seed_kmer, w) : UINT64_MAX;
  int nsites = 0;
  for (const auto &seq : pos_seqs) {
    uint64_t seqlen = seq.size();
    if (seqlen < w) continue;
    for (uint64_t off = 0; off + w <= seqlen; off++) {
      uint64_t kmer = kmer_encode4(seq, w, off);
      if (kmer == UINT64_MAX) continue;
      if (hamming4(kmer, seed_kmer, w) <= (uint64_t) HAMMING_MISMATCH) {
        for (uint64_t i = 0; i < w; i++) {
          uint8_t idx = CHAR2IDX[(unsigned char) seq[off + i]];
          if (idx < 4) ppm[i][idx]++;
        }
        nsites++;
      } else if (rc_seed != UINT64_MAX &&
                 hamming4(kmer, rc_seed, w) <= (uint64_t) HAMMING_MISMATCH) {
        // RC match: add reverse-complement of the window
        for (uint64_t i = 0; i < w; i++) {
          uint8_t idx = CHAR2IDX[(unsigned char) seq[off + w - 1 - i]];
          if (idx < 4) ppm[i][COMP4[idx]]++;
        }
        nsites++;
      }
    }
  }
  return nsites;
}

// -----------------------------------------------------------------------------
// Refine motif (port of refine_motif). Returns new nsites or 0 if bailed.
// -----------------------------------------------------------------------------
static int refine_motif(const std::vector<std::string> &pos_seqs,
                        Motif &m, uint64_t w, bool RC,
                        const std::array<double, 4> &bkg,
                        int pseudocount, double hit_pval) {
  if (m.threshold == INT_MAX) return 0;
  int ppm[MAX_MOTIF_WIDTH][4];
  std::memset(ppm, 0, sizeof(ppm));
  int nsites = 0;
  const int thr = m.threshold - 1;

  for (const auto &seq : pos_seqs) {
    uint64_t seqlen = seq.size();
    if (seqlen < w) continue;
    for (uint64_t off = 0; off + w <= seqlen; off++) {
      if (RC) {
        int s, src;
        score_subseq_rc(m, seq, off, &s, &src);
        if (s > thr) {
          for (uint64_t i = 0; i < w; i++) {
            uint8_t idx = CHAR2IDX[(unsigned char) seq[off + i]];
            if (idx < 4) ppm[i][idx]++;
          }
          nsites++;
        }
        if (src > thr) {
          for (uint64_t i = 0; i < w; i++) {
            uint8_t idx = CHAR2IDX[(unsigned char) seq[off + w - 1 - i]];
            if (idx < 4) ppm[i][COMP4[idx]]++;
          }
          nsites++;
        }
      } else {
        int s = score_subseq(m, seq, off);
        if (s > thr) {
          for (uint64_t i = 0; i < w; i++) {
            uint8_t idx = CHAR2IDX[(unsigned char) seq[off + i]];
            if (idx < 4) ppm[i][idx]++;
          }
          nsites++;
        }
      }
    }
  }

  if (nsites < MIN_REFINE_HITS) return 0;
  convert_ppm_to_motif(m, ppm, w, nsites, bkg, pseudocount, hit_pval);
  return nsites;
}

// -----------------------------------------------------------------------------
// Evaluate motif (port of evaluate_motif). Fills coverage bitmask over pos.
// -----------------------------------------------------------------------------
static double evaluate_motif(const std::vector<std::string> &pos_seqs,
                             const std::vector<std::string> &neg_seqs,
                             const Motif &m, bool RC,
                             uint64_t &sites_pos_out, uint64_t &sites_neg_out,
                             uint64_t &seqs_pos_out,  uint64_t &seqs_neg_out,
                             std::vector<std::vector<uint8_t>> *covered) {
  const int thr = m.threshold - 1;
  uint64_t seqs_pos = 0, seqs_neg = 0;
  uint64_t sites_pos = 0, sites_neg = 0;
  uint64_t w = m.size;

  for (uint64_t si = 0; si < pos_seqs.size(); si++) {
    const auto &seq = pos_seqs[si];
    uint64_t seqlen = seq.size();
    if (seqlen < w) continue;
    int has_hit = 0;
    for (uint64_t off = 0; off + w <= seqlen; off++) {
      if (RC) {
        int s, src;
        score_subseq_rc(m, seq, off, &s, &src);
        if (s > thr) {
          has_hit = 1; sites_pos++;
          if (covered)
            for (uint64_t p = off; p < off + w; p++)
              (*covered)[si][p / 8] |= (uint8_t) (1u << (p % 8));
        }
        if (src > thr) {
          has_hit = 1; sites_pos++;
          if (covered)
            for (uint64_t p = off; p < off + w; p++)
              (*covered)[si][p / 8] |= (uint8_t) (1u << (p % 8));
        }
      } else {
        int s = score_subseq(m, seq, off);
        if (s > thr) {
          has_hit = 1; sites_pos++;
          if (covered)
            for (uint64_t p = off; p < off + w; p++)
              (*covered)[si][p / 8] |= (uint8_t) (1u << (p % 8));
        }
      }
    }
    if (has_hit) seqs_pos++;
  }
  for (const auto &seq : neg_seqs) {
    uint64_t seqlen = seq.size();
    if (seqlen < w) continue;
    int has_hit = 0;
    for (uint64_t off = 0; off + w <= seqlen; off++) {
      if (RC) {
        int s, src;
        score_subseq_rc(m, seq, off, &s, &src);
        if (s > thr)   { has_hit = 1; sites_neg++; }
        if (src > thr) { has_hit = 1; sites_neg++; }
      } else {
        int s = score_subseq(m, seq, off);
        if (s > thr)   { has_hit = 1; sites_neg++; }
      }
    }
    if (has_hit) seqs_neg++;
  }

  sites_pos_out = sites_pos; sites_neg_out = sites_neg;
  seqs_pos_out  = seqs_pos;  seqs_neg_out  = seqs_neg;

  uint64_t n_pos = pos_seqs.size(), n_neg = neg_seqs.size();
  return fishers_exact_log_greater(seqs_pos, seqs_neg,
           n_pos > seqs_pos ? n_pos - seqs_pos : 0,
           n_neg > seqs_neg ? n_neg - seqs_neg : 0);
}

// -----------------------------------------------------------------------------
// Mask covered positions in the per-thread positives buffer (port of
// mask_positions). Writes 'N' over covered bases so the next iteration
// does not re-discover the same motif.
// -----------------------------------------------------------------------------
static void mask_positions(std::vector<std::string> &mut_seqs,
                           const std::vector<std::vector<uint8_t>> &covered) {
  for (uint64_t si = 0; si < mut_seqs.size(); si++) {
    if (covered[si].empty()) continue;
    uint64_t seqlen = mut_seqs[si].size();
    uint64_t nbytes = (seqlen + 7) / 8;
    for (uint64_t b = 0; b < nbytes; b++) {
      if (!covered[si][b]) continue;
      for (int bit = 0; bit < 8; bit++) {
        if (covered[si][b] & (1u << bit)) {
          uint64_t pos = b * 8 + bit;
          if (pos < seqlen) mut_seqs[si][pos] = 'N';
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
// IUPAC consensus (port of build_consensus / consensus_char).
// -----------------------------------------------------------------------------
static char iupac_from_mask(int mask) {
  switch (mask) {
    case 0x1: return 'A'; case 0x2: return 'C'; case 0x4: return 'G';
    case 0x8: return 'T';
    case 0x3: return 'M'; case 0x5: return 'R'; case 0x9: return 'W';
    case 0x6: return 'S'; case 0xA: return 'Y'; case 0xC: return 'K';
    case 0x7: return 'V'; case 0xB: return 'H'; case 0xD: return 'D';
    case 0xE: return 'B'; default:  return 'N';
  }
}
static char consensus_char(const std::array<double, 4> &probs) {
  int order[4] = {0, 1, 2, 3};
  // Insertion sort descending by probability.
  for (int i = 1; i < 4; i++) {
    int j = i;
    while (j > 0 && probs[order[j]] > probs[order[j - 1]]) {
      std::swap(order[j], order[j - 1]);
      j--;
    }
  }
  double cum = 0.0;
  int mask = 0;
  for (int k = 0; k < 4; k++) {
    cum += probs[order[k]];
    mask |= 1 << order[k];
    if (cum >= CONSENSUS_THRESHOLD) break;
  }
  return iupac_from_mask(mask);
}
static std::string build_consensus(const Motif &m) {
  std::string out(m.size, 'N');
  for (uint64_t j = 0; j < m.size; j++) out[j] = consensus_char(m.pwm_probs[j]);
  return out;
}

// -----------------------------------------------------------------------------
// IC trim (port of ic_trim_motif). Trims low-IC flanks in-place.
// -----------------------------------------------------------------------------
static double column_ic(const std::array<double, 4> &probs) {
  double ic = 2.0;  // log2(4)
  for (int i = 0; i < 4; i++) {
    if (probs[i] > 0) ic += probs[i] * std::log2(probs[i]);
  }
  return ic < 0 ? 0 : ic;
}
static void ic_trim_motif(Motif &m) {
  uint64_t left = 0;
  while (left < m.size && column_ic(m.pwm_probs[left]) < MIN_IC_BITS) left++;
  uint64_t right = m.size;
  while (right > left && column_ic(m.pwm_probs[right - 1]) < MIN_IC_BITS) right--;
  uint64_t new_size = right - left;
  if (new_size < 3) return;
  if (left == 0 && right == m.size) return;
  for (uint64_t i = 0; i < new_size; i++) m.pwm_probs[i] = m.pwm_probs[i + left];
  m.pwm_probs.resize(new_size);
  m.size = new_size;
}

// -----------------------------------------------------------------------------
// Discovery result (port of disc_result_t).
// -----------------------------------------------------------------------------
struct DiscResult {
  Motif    motif;
  std::vector<std::vector<uint8_t>> covered;
  double   pvalue;
  double   qvalue;
  uint64_t sites_pos, sites_neg;
  uint64_t seqs_pos,  seqs_neg;
  uint64_t width;
  uint64_t discovery_seq;
  bool     dropped;
};

// -----------------------------------------------------------------------------
// Per-width discovery loop (port of discover_for_width). Operates on a
// per-thread mutable copy of positives; appends to local_results.
// -----------------------------------------------------------------------------
static void discover_for_width(const std::vector<std::string> &pristine_pos,
                               const std::vector<std::string> &neg_seqs,
                               uint64_t w, int nmotifs, bool RC,
                               const std::array<double, 4> &bkg,
                               int pseudocount, double hit_pval, double stop_pval,
                               std::vector<DiscResult> &local_results) {
  std::vector<std::string> mut_pos = pristine_pos;  // mutable copy for masking

  uint64_t intra = 0;
  for (int ki = 0; ki < nmotifs; ki++) {
    std::vector<Seed> seeds;
    if (enumerate_seeds(mut_pos, neg_seqs, w, RC, seeds) != 0) return;
    if (seeds.empty()) break;

    bool accepted = false;
    for (uint64_t si = 0; si < seeds.size() && !accepted; si++) {
      int ppm[MAX_MOTIF_WIDTH][4];
      std::memset(ppm, 0, sizeof(ppm));
      int nsites = build_ppm_from_seed(mut_pos, seeds[si].kmer, w, RC, ppm);
      if (nsites < MIN_REFINE_HITS) continue;

      Motif m;
      convert_ppm_to_motif(m, ppm, w, nsites, bkg, pseudocount, hit_pval);
      if (m.threshold == INT_MAX) continue;

      for (int ri = 0; ri < REFINE_PASSES; ri++) {
        if (refine_motif(mut_pos, m, w, RC, bkg, pseudocount, hit_pval) == 0)
          break;
      }
      if (m.threshold == INT_MAX) continue;

      // Allocate coverage bitmask
      std::vector<std::vector<uint8_t>> covered(mut_pos.size());
      for (uint64_t sj = 0; sj < mut_pos.size(); sj++) {
        uint64_t nbytes = (mut_pos[sj].size() + 7) / 8;
        covered[sj].assign(nbytes, 0);
      }

      uint64_t sp, sn, ep, en;
      double pval = evaluate_motif(mut_pos, neg_seqs, m, RC, sp, sn, ep, en, &covered);
      if (pval > stop_pval) continue;

      DiscResult r;
      r.motif         = std::move(m);
      r.covered       = std::move(covered);
      r.pvalue        = pval;
      r.qvalue        = 1.0;
      r.sites_pos     = sp;
      r.sites_neg     = sn;
      r.seqs_pos      = ep;
      r.seqs_neg      = en;
      r.width         = w;
      r.discovery_seq = intra;
      r.dropped       = false;
      local_results.push_back(std::move(r));

      mask_positions(mut_pos, local_results.back().covered);
      accepted = true;
      intra++;
    }
    if (!accepted) break;
  }
}

// -----------------------------------------------------------------------------
// Rebuild coverage (port of rebuild_coverage). After per-thread discovery,
// re-evaluate each motif on PRISTINE positives so cross-width dedup
// compares apples-to-apples.
// -----------------------------------------------------------------------------
static void rebuild_coverage(const std::vector<std::string> &pristine_pos,
                             const std::vector<std::string> &neg_seqs,
                             bool RC, std::vector<DiscResult> &results) {
  for (auto &r : results) {
    for (uint64_t si = 0; si < pristine_pos.size(); si++) {
      uint64_t nbytes = (pristine_pos[si].size() + 7) / 8;
      r.covered[si].assign(nbytes, 0);
    }
    uint64_t sp, sn, ep, en;
    double pval = evaluate_motif(pristine_pos, neg_seqs, r.motif, RC,
                                 sp, sn, ep, en, &r.covered);
    r.pvalue    = pval;
    r.sites_pos = sp;
    r.sites_neg = sn;
    r.seqs_pos  = ep;
    r.seqs_neg  = en;
  }
}

// -----------------------------------------------------------------------------
// Cross-width dedup (port of dedup_results). Drops motifs whose coverage
// overlaps too much with a lower-pvalue motif.
// -----------------------------------------------------------------------------
static void dedup_results(std::vector<DiscResult> &results, double dedup_overlap) {
  if (results.size() < 2) return;
  std::sort(results.begin(), results.end(),
            [](const DiscResult &a, const DiscResult &b) {
              return a.pvalue < b.pvalue;
            });
  for (uint64_t i = 0; i < results.size(); i++) {
    if (results[i].dropped) continue;
    for (uint64_t j = i + 1; j < results.size(); j++) {
      if (results[j].dropped) continue;
      uint64_t inter = 0, sz_i = 0, sz_j = 0;
      uint64_t nseqs = results[i].covered.size();
      for (uint64_t si = 0; si < nseqs; si++) {
        uint64_t nbytes = results[i].covered[si].size();
        for (uint64_t b = 0; b < nbytes; b++) {
          uint8_t ci = results[i].covered[si][b];
          uint8_t cj = results[j].covered[si][b];
          inter += (uint64_t) __builtin_popcount(ci & cj);
          sz_i  += (uint64_t) __builtin_popcount(ci);
          sz_j  += (uint64_t) __builtin_popcount(cj);
        }
      }
      uint64_t sz_min = sz_i < sz_j ? sz_i : sz_j;
      if (sz_min > 0 && (double) inter / sz_min > dedup_overlap)
        results[j].dropped = true;
    }
  }
}

// -----------------------------------------------------------------------------
// Top-level Rcpp entry: runs the full pipeline.
// -----------------------------------------------------------------------------
// [[Rcpp::export(rng = false)]]
Rcpp::List motif_finder_cpp(const std::vector<std::string> &pos_seqs,
                            const std::vector<std::string> &neg_seqs,
                            int min_w, int max_w, int nmotifs,
                            double hit_pvalue, double stop_pvalue,
                            double dedup_overlap, bool RC,
                            Rcpp::NumericVector bkg_in,
                            int pseudocount, int nthreads) {
  init_char2idx();

  if (min_w < 3 || max_w > MAX_MOTIF_WIDTH || min_w > max_w)
    Rcpp::stop("Invalid width range: need 3 <= min.width <= max.width <= %d.",
               MAX_MOTIF_WIDTH);
  if (bkg_in.size() != 4)
    Rcpp::stop("bkg must be a length-4 numeric vector (A,C,G,T).");

  std::array<double, 4> bkg = {bkg_in[0], bkg_in[1], bkg_in[2], bkg_in[3]};

  uint64_t n_widths = (uint64_t) (max_w - min_w + 1);
  // Per-width local results; threading picks widths off this.
  std::vector<std::vector<DiscResult>> per_width(n_widths);

  RcppThread::parallelFor(0, n_widths,
    [&](uint64_t wi) {
      uint64_t w = (uint64_t) min_w + wi;
      discover_for_width(pos_seqs, neg_seqs, w, nmotifs, RC, bkg,
                         pseudocount, hit_pvalue, stop_pvalue,
                         per_width[wi]);
    }, nthreads);

  // Merge in deterministic (width, discovery_seq) order
  std::vector<DiscResult> merged;
  for (uint64_t wi = 0; wi < n_widths; wi++)
    for (auto &r : per_width[wi]) merged.push_back(std::move(r));
  per_width.clear();

  // Sort by (width, discovery_seq)
  std::sort(merged.begin(), merged.end(),
            [](const DiscResult &a, const DiscResult &b) {
              if (a.width != b.width) return a.width < b.width;
              return a.discovery_seq < b.discovery_seq;
            });

  if (merged.empty()) {
    return Rcpp::List::create();
  }

  // Rebuild coverage on pristine positives
  rebuild_coverage(pos_seqs, neg_seqs, RC, merged);

  // Cross-width dedup (this also sorts by pvalue ascending)
  dedup_results(merged, dedup_overlap);

  // Collect survivors in p-value order; BH q-values across survivors only
  std::vector<DiscResult> survivors;
  for (auto &r : merged) if (!r.dropped) survivors.push_back(std::move(r));
  std::vector<double> pvals(survivors.size());
  for (uint64_t i = 0; i < survivors.size(); i++) pvals[i] = survivors[i].pvalue;
  std::vector<double> qvals;
  bh_qvalues(pvals, qvals);
  for (uint64_t i = 0; i < survivors.size(); i++) survivors[i].qvalue = qvals[i];

  // IC-trim flanks on the final survivors
  for (auto &r : survivors) ic_trim_motif(r.motif);

  // Build the Rcpp::List return: one element per motif.
  Rcpp::List out(survivors.size());
  for (uint64_t i = 0; i < survivors.size(); i++) {
    const DiscResult &r = survivors[i];
    const Motif &m = r.motif;
    Rcpp::NumericMatrix ppm((int) m.size, 4);
    for (uint64_t j = 0; j < m.size; j++) {
      for (int k = 0; k < 4; k++) ppm(j, k) = m.pwm_probs[j][k];
    }
    out[i] = Rcpp::List::create(
      Rcpp::Named("ppm")           = ppm,
      Rcpp::Named("width")         = (int) m.size,
      Rcpp::Named("nsites")        = (double) m.nsites_actual,
      Rcpp::Named("consensus")     = build_consensus(m),
      Rcpp::Named("seqs_pos")      = (double) r.seqs_pos,
      Rcpp::Named("seqs_neg")      = (double) r.seqs_neg,
      Rcpp::Named("sites_pos")     = (double) r.sites_pos,
      Rcpp::Named("sites_neg")     = (double) r.sites_neg,
      Rcpp::Named("n_pos")         = (double) pos_seqs.size(),
      Rcpp::Named("n_neg")         = (double) neg_seqs.size(),
      Rcpp::Named("pvalue")        = r.pvalue,
      Rcpp::Named("qvalue")        = r.qvalue
    );
  }
  return out;
}
