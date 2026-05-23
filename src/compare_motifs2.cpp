#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>
#include "types.h"

/* compare_motifs2.cpp
 * --------------------
 * yamtk cmp-aligned motif comparison primitives.
 *
 * This file is a port of the algorithm in yamtk/src/yamcmp.c (Tremblay 2026)
 * adapted to operate on Rcpp inputs in a thread-safe way (pre-extract
 * Rcpp objects in serial, parallelise over queries with RcppThread, wrap
 * in serial).
 *
 * Two exports:
 *   compare_motifs2_align_cpp:
 *     given a flat list of motif PPMs and per-pair query/target index
 *     vectors, find the best PCC alignment per pair and return
 *     score/offset/strand/overlap/n_tested/q_start_oriented.
 *   compare_motifs2_pvalue_cpp:
 *     given the same motif list plus per-pair alignment metadata, build
 *     the yamtk-style null PMF (empirical or parametric Dirichlet-
 *     Multinomial), convolve per-(start, L) for each unique query, and
 *     return Bonferroni-corrected p-values per pair.
 *
 * The conventions:
 *   - Motifs are PPMs, position-major (mat[i][j] = probability of base j
 *     at position i). 4-base alphabet only (DNA/RNA, treated as ACGT).
 *   - PCC is centred on a uniform 0.25 background (yamcmp.c:1028-1039).
 *   - PCC bins: 51, spanning [-1, 1] with step 0.04 (PCC_BINS).
 *   - Parametric grid: K=5 simplex, 56 columns (ENUM_N_COLS); Dirichlet
 *     concentration alpha = 4 * bkg (ENUM_DIRICHLET_N).
 *   - Consensus threshold for IUPAC: 0.75 (CONSENSUS_THRESHOLD).
 */

static const int    PCC_BINS              = 51;
static const int    ENUM_GRID_K           = 5;
static const int    ENUM_DIRICHLET_N      = 4;
static const int    ENUM_N_COLS           = 56;   /* C(K+3, 3) for K=5 */
static const double CONSENSUS_THRESHOLD   = 0.75;

/* ============================================================
 * Column PCC and discretization (yamcmp.c:1028-1079)
 * ============================================================ */

static inline double pcc_column(const double *q, const double *t) {
  double qd[4], td[4];
  double sqq = 0.0, stt = 0.0, sqt = 0.0;
  for (int i = 0; i < 4; ++i) { qd[i] = q[i] - 0.25; td[i] = t[i] - 0.25; }
  for (int i = 0; i < 4; ++i) {
    sqq += qd[i] * qd[i];
    stt += td[i] * td[i];
    sqt += qd[i] * td[i];
  }
  if (sqq < 1e-12 || stt < 1e-12) return 0.0;
  return sqt / std::sqrt(sqq * stt);
}

/* PCC value in [-1, 1] -> bin index in [0, PCC_BINS-1] */
static inline int score_to_bin_pcc(double s) {
  double t = (s + 1.0) / 2.0;                 /* [-1,1] -> [0,1] */
  if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
  int bin = (int) std::floor(t * (PCC_BINS - 1) + 0.5);
  if (bin < 0) bin = 0;
  if (bin > PCC_BINS - 1) bin = PCC_BINS - 1;
  return bin;
}

/* Sum-of-L PCCs -> bin in [0, L*(PCC_BINS-1)]. yamcmp.c:1072-1079 */
static inline int sum_score_to_bin_pcc(double s, int L) {
  double per_col = (s / (double)L + 1.0) / 2.0;
  if (per_col < 0.0) per_col = 0.0; else if (per_col > 1.0) per_col = 1.0;
  int bin = (int) std::floor(per_col * (PCC_BINS - 1) * L + 0.5);
  int hi = L * (PCC_BINS - 1);
  if (bin < 0) bin = 0;
  if (bin > hi) bin = hi;
  return bin;
}

/* ============================================================
 * IUPAC consensus (yamcmp.c:986-1016)
 * ============================================================ */

static inline char iupac_from_mask(int mask) {
  switch (mask) {
    case 0x1: return 'A'; case 0x2: return 'C';
    case 0x4: return 'G'; case 0x8: return 'T';
    case 0x3: return 'M'; case 0x5: return 'R';
    case 0x9: return 'W'; case 0x6: return 'S';
    case 0xA: return 'Y'; case 0xC: return 'K';
    case 0x7: return 'V'; case 0xB: return 'H';
    case 0xD: return 'D'; case 0xE: return 'B';
    default:  return 'N';
  }
}

static char consensus_char(const double *probs) {
  int order[4] = {0, 1, 2, 3};
  for (int i = 1; i < 4; ++i) {
    int j = i;
    while (j > 0 && probs[order[j]] > probs[order[j-1]]) {
      int tmp = order[j]; order[j] = order[j-1]; order[j-1] = tmp;
      --j;
    }
  }
  double cum = 0.0;
  int mask = 0;
  for (int k = 0; k < 4; ++k) {
    cum += probs[order[k]];
    mask |= 1 << order[k];
    if (cum >= CONSENSUS_THRESHOLD) break;
  }
  return iupac_from_mask(mask);
}

/* Build IUPAC string for a column range of a motif. */
static std::string consensus_range(const list_num_t &mot,
                                   int start, int L) {
  std::string s;
  s.resize((std::size_t) L);
  for (int j = 0; j < L; ++j) s[j] = consensus_char(&mot[start + j][0]);
  return s;
}

/* ============================================================
 * RC of a PPM (DNA/RNA; A<->T, C<->G, reverse columns)
 * ============================================================ */

static list_num_t rc_motif(const list_num_t &mot) {
  std::size_t w = mot.size();
  list_num_t out(w, vec_num_t(4, 0.0));
  for (std::size_t i = 0; i < w; ++i) {
    /* Reverse complement: position w-1-i, and swap A<->T, C<->G. */
    out[w - 1 - i][0] = mot[i][3];   /* A <- T */
    out[w - 1 - i][1] = mot[i][2];   /* C <- G */
    out[w - 1 - i][2] = mot[i][1];   /* G <- C */
    out[w - 1 - i][3] = mot[i][0];   /* T <- A */
  }
  return out;
}

/* ============================================================
 * Parametric grid (yamcmp.c:1092-1133)
 * Build once per call, given user-supplied background.
 * ============================================================ */

struct enum_col_t {
  double col[4];
  double weight;   /* sum over all enum cols == 1 after normalisation */
};

static std::vector<enum_col_t> build_enum_cols(const vec_num_t &bkg) {
  std::vector<enum_col_t> enum_cols;
  enum_cols.reserve(ENUM_N_COLS);

  double alpha[4];
  for (int i = 0; i < 4; ++i) alpha[i] = (double) ENUM_DIRICHLET_N * bkg[i];
  for (int i = 0; i < 4; ++i) if (alpha[i] < 0.1) alpha[i] = 0.1;
  double alpha_sum = alpha[0] + alpha[1] + alpha[2] + alpha[3];

  for (int a = 0; a <= ENUM_GRID_K; ++a)
    for (int b = 0; b <= ENUM_GRID_K - a; ++b)
      for (int c = 0; c <= ENUM_GRID_K - a - b; ++c) {
        int d = ENUM_GRID_K - a - b - c;
        int counts[4] = {a, b, c, d};
        enum_col_t e;
        e.col[0] = (double) a / ENUM_GRID_K;
        e.col[1] = (double) b / ENUM_GRID_K;
        e.col[2] = (double) c / ENUM_GRID_K;
        e.col[3] = (double) d / ENUM_GRID_K;
        double lw = std::lgamma((double) ENUM_GRID_K + 1.0);
        for (int i = 0; i < 4; ++i) {
          lw -= std::lgamma((double) counts[i] + 1.0);
          lw += std::lgamma((double) counts[i] + alpha[i]) - std::lgamma(alpha[i]);
        }
        lw -= std::lgamma((double) ENUM_GRID_K + alpha_sum) - std::lgamma(alpha_sum);
        e.weight = lw;   /* store log temporarily */
        enum_cols.push_back(e);
      }

  /* log-sum-exp normalisation (yamcmp.c:1124-1132) */
  double max_lw = enum_cols[0].weight;
  for (std::size_t i = 1; i < enum_cols.size(); ++i)
    if (enum_cols[i].weight > max_lw) max_lw = enum_cols[i].weight;
  double zsum = 0.0;
  for (auto &e : enum_cols) { e.weight = std::exp(e.weight - max_lw); zsum += e.weight; }
  for (auto &e : enum_cols) e.weight /= zsum;
  return enum_cols;
}

/* ============================================================
 * Per-column null PMF builders (yamcmp.c:1137-1166)
 * ============================================================ */

/* Empirical null: histogram of PCC(q_col, t_col) over all target columns.
 * `target_cols_flat` is a flat list of target columns
 * (each is a vec_num_t of length 4). */
static vec_num_t build_col_pmf_empirical(const double *q_col,
                                         const list_num_t &target_cols_flat) {
  vec_num_t pmf(PCC_BINS, 0.0);
  if (target_cols_flat.empty()) return pmf;
  for (const auto &tc : target_cols_flat) {
    double v = pcc_column(q_col, &tc[0]);
    pmf[score_to_bin_pcc(v)] += 1.0;
  }
  double inv = 1.0 / (double) target_cols_flat.size();
  for (auto &p : pmf) p *= inv;
  return pmf;
}

/* Parametric null: weighted sum over the K=5 simplex grid. */
static vec_num_t build_col_pmf_parametric(const double *q_col,
                                          const std::vector<enum_col_t> &enum_cols) {
  vec_num_t pmf(PCC_BINS, 0.0);
  for (const auto &e : enum_cols) {
    double v = pcc_column(q_col, e.col);
    pmf[score_to_bin_pcc(v)] += e.weight;
  }
  return pmf;
}

/* ============================================================
 * Incremental convolution (yamcmp.c:1170-1179)
 * out_len = a_len + PCC_BINS - 1
 * ============================================================ */

static vec_num_t convolve_with_col(const vec_num_t &a, const vec_num_t &b) {
  std::size_t a_len = a.size();
  std::size_t b_len = b.size();
  vec_num_t out(a_len + b_len - 1, 0.0);
  for (std::size_t i = 0; i < a_len; ++i) {
    double ai = a[i];
    if (ai == 0.0) continue;
    for (std::size_t j = 0; j < b_len; ++j) out[i + j] += ai * b[j];
  }
  return out;
}

/* Tail probability P(S' >= observed_bin) (yamcmp.c:1230-1240) */
static double tail_prob(const vec_num_t &pmf, int observed_bin) {
  int len = (int) pmf.size();
  if (observed_bin < 0) observed_bin = 0;
  if (observed_bin >= len) return 0.0;
  double s = 0.0;
  for (int i = observed_bin; i < len; ++i) s += pmf[i];
  if (s < 0.0) s = 0.0;
  if (s > 1.0) s = 1.0;
  return s;
}

/* ============================================================
 * Alignment scoring (yamcmp.c:1247-1329)
 * ============================================================ */

struct align_result_t {
  double score;
  int    L;
  int    offset;             /* d_in_q: query column corresponding to target col 0 */
  int    strand;             /* 0 = '+', 1 = '-' (RC query) */
  int    q_start_oriented;   /* start col in the oriented query (for PMF lookup) */
  int    q_start_in_orig;    /* start col in the original (forward) query (for output) */
  int    t_start;
  int    n_tested;
};

static double align_score_pcc(const list_num_t &q, int wq,
                              const list_num_t &t, int wt,
                              int d_in_q,
                              int *out_L, int *out_q_start, int *out_t_start) {
  int t_start = (d_in_q >= 0) ? 0 : -d_in_q;
  int t_end   = (wq - d_in_q < wt) ? (wq - d_in_q) : wt;
  int L = t_end - t_start;
  if (L <= 0) { *out_L = 0; return 0.0; }
  int q_start = d_in_q + t_start;
  double s = 0.0;
  for (int j = 0; j < L; ++j) {
    s += pcc_column(&q[q_start + j][0], &t[t_start + j][0]);
  }
  *out_L = L;
  *out_q_start = q_start;
  *out_t_start = t_start;
  return s;
}

static align_result_t best_alignment_pcc(const list_num_t &q,
                                         const list_num_t &q_rc,
                                         bool have_rc,
                                         const list_num_t &t,
                                         int min_overlap) {
  int wq = (int) q.size();
  int wt = (int) t.size();
  int mo = min_overlap;
  if (mo > wq) mo = wq;
  if (mo > wt) mo = wt;
  int d_lo = mo - wt;
  int d_hi = wq - mo;

  align_result_t best;
  best.score             = -std::numeric_limits<double>::infinity();
  best.L                 = 0;
  best.offset            = 0;
  best.strand            = 0;
  best.q_start_oriented  = 0;
  best.q_start_in_orig   = 0;
  best.t_start           = 0;
  best.n_tested          = 0;

  int orientations = have_rc ? 2 : 1;
  for (int orient = 0; orient < orientations; ++orient) {
    const list_num_t &qq = (orient == 0) ? q : q_rc;
    for (int d = d_lo; d <= d_hi; ++d) {
      int L, q_start, t_start;
      double s = align_score_pcc(qq, wq, t, wt, d, &L, &q_start, &t_start);
      if (L < mo) continue;
      best.n_tested++;
      bool better = (s > best.score) || best.n_tested == 1;
      if (better) {
        best.score            = s;
        best.L                = L;
        best.offset           = d;
        best.strand           = orient;
        best.q_start_oriented = q_start;
        if (orient == 0) {
          best.q_start_in_orig = q_start;
        } else {
          best.q_start_in_orig = wq - 1 - (q_start + L - 1);
        }
        best.t_start = t_start;
      }
    }
  }
  return best;
}

/* ============================================================
 * Helpers to extract list<NumericMatrix> -> vector<list_num_t>
 * (position-major). Done serially before parallel regions.
 * ============================================================ */

static std::vector<list_num_t> extract_motif_list(const Rcpp::List &mots) {
  std::vector<list_num_t> out((std::size_t) mots.size());
  for (R_xlen_t m = 0; m < mots.size(); ++m) {
    Rcpp::NumericMatrix mat = mots[m];
    std::size_t w = (std::size_t) mat.ncol();
    out[(std::size_t) m].assign(w, vec_num_t(4, 0.0));
    for (std::size_t i = 0; i < w; ++i)
      for (std::size_t j = 0; j < 4; ++j)
        out[(std::size_t) m][i][j] = mat(j, i);
  }
  return out;
}

/* Flatten all columns of all targets into a single list_num_t for empirical
 * PMF construction. */
static list_num_t flatten_target_cols(const std::vector<list_num_t> &targets) {
  std::size_t total = 0;
  for (const auto &t : targets) total += t.size();
  list_num_t flat;
  flat.reserve(total);
  for (const auto &t : targets)
    for (const auto &col : t) flat.push_back(col);
  return flat;
}

/* ============================================================
 * EXPORT 1: per-pair best alignment, no p-value.
 * ============================================================ */

// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame compare_motifs2_align_cpp(
    const Rcpp::List &mots,
    const Rcpp::IntegerVector &qi,
    const Rcpp::IntegerVector &ti,
    const int min_overlap,
    const bool RC,
    const int nthreads = 1) {

  R_xlen_t n_pairs = qi.size();
  if (ti.size() != n_pairs)
    Rcpp::stop("compare_motifs2_align_cpp: qi and ti must have equal length");

  /* Serial pre-extract */
  std::vector<list_num_t> motifs = extract_motif_list(mots);
  std::size_t n_mot = motifs.size();
  std::vector<list_num_t> motifs_rc(n_mot);
  if (RC) {
    for (std::size_t i = 0; i < n_mot; ++i) motifs_rc[i] = rc_motif(motifs[i]);
  }

  /* Convert 1-based -> 0-based pair indices */
  std::vector<int> qi_0(n_pairs), ti_0(n_pairs);
  for (R_xlen_t k = 0; k < n_pairs; ++k) {
    qi_0[(std::size_t) k] = qi[k] - 1;
    ti_0[(std::size_t) k] = ti[k] - 1;
  }

  std::vector<double> out_score(n_pairs);
  std::vector<int>    out_L(n_pairs), out_offset(n_pairs),
                      out_strand(n_pairs), out_qstart(n_pairs),
                      out_ntested(n_pairs);

  RcppThread::parallelFor(0, (std::size_t) n_pairs,
    [&](std::size_t k) {
      int qix = qi_0[k];
      int tix = ti_0[k];
      const list_num_t &q  = motifs[(std::size_t) qix];
      const list_num_t &qr = RC ? motifs_rc[(std::size_t) qix]
                                : motifs[(std::size_t) qix];
      const list_num_t &t  = motifs[(std::size_t) tix];
      align_result_t a = best_alignment_pcc(q, qr, RC, t, min_overlap);
      out_score[k]   = a.score;
      out_L[k]       = a.L;
      out_offset[k]  = a.offset;
      out_strand[k]  = a.strand;
      out_qstart[k]  = a.q_start_oriented;
      out_ntested[k] = a.n_tested;
    }, nthreads);

  return Rcpp::DataFrame::create(
    Rcpp::_["score"]            = out_score,
    Rcpp::_["overlap"]          = out_L,
    Rcpp::_["offset"]           = out_offset,
    Rcpp::_["strand"]           = out_strand,
    Rcpp::_["q_start_oriented"] = out_qstart,
    Rcpp::_["n_tested"]         = out_ntested,
    Rcpp::_["stringsAsFactors"] = false
  );
}

/* ============================================================
 * EXPORT 2: per-pair Bonferroni p-value (yamtk-style).
 *
 * Inputs (per pair k):
 *   qi[k], ti[k], qstart[k], overlap[k], strand[k], score[k], n_tested[k]
 * Modes:
 *   null_mode = 0 -> empirical from target_mats
 *   null_mode = 1 -> parametric Dirichlet-Multinomial over K=5 simplex
 * Returns a NumericVector of Bonferroni-corrected p-values, one per pair.
 *
 * Parallelisation: over unique queries. Each thread:
 *   - builds per-column PMFs for its query (and its RC if any pair uses '-')
 *   - builds the per-(start, L) convolution cache
 *   - looks up tail probability for each of its pairs
 *   - frees the cache when done
 * ============================================================ */

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector compare_motifs2_pvalue_cpp(
    const Rcpp::List &query_mats,
    const Rcpp::List &target_mats,
    const Rcpp::IntegerVector &qi,
    const Rcpp::IntegerVector &qstart,
    const Rcpp::IntegerVector &overlap,
    const Rcpp::IntegerVector &strand,
    const Rcpp::NumericVector &score,
    const Rcpp::IntegerVector &n_tested,
    const Rcpp::NumericVector &bkg,
    const int null_mode,
    const int nthreads = 1) {

  R_xlen_t n_pairs = qi.size();
  if (qstart.size()   != n_pairs ||
      overlap.size()  != n_pairs ||
      strand.size()   != n_pairs ||
      score.size()    != n_pairs ||
      n_tested.size() != n_pairs)
    Rcpp::stop("compare_motifs2_pvalue_cpp: input vectors must be the same length");
  if (bkg.size() != 4)
    Rcpp::stop("compare_motifs2_pvalue_cpp: bkg must have length 4");
  if (null_mode != 0 && null_mode != 1)
    Rcpp::stop("compare_motifs2_pvalue_cpp: null_mode must be 0 (empirical) or 1 (parametric)");

  /* Serial pre-extract */
  std::vector<list_num_t> queries = extract_motif_list(query_mats);
  std::vector<list_num_t> targets = extract_motif_list(target_mats);
  std::size_t n_q = queries.size();

  std::vector<list_num_t> queries_rc(n_q);
  for (std::size_t i = 0; i < n_q; ++i) queries_rc[i] = rc_motif(queries[i]);

  vec_num_t bkg_v(4);
  for (int i = 0; i < 4; ++i) bkg_v[i] = bkg[i];

  /* For empirical mode, flatten all target columns once. */
  list_num_t target_cols_flat;
  if (null_mode == 0) target_cols_flat = flatten_target_cols(targets);

  /* For parametric mode, build the enum grid once. */
  std::vector<enum_col_t> enum_cols;
  if (null_mode == 1) enum_cols = build_enum_cols(bkg_v);

  /* 1-based -> 0-based; group pair indices by query index. */
  std::vector<std::vector<std::size_t>> pairs_by_query(n_q);
  for (R_xlen_t k = 0; k < n_pairs; ++k) {
    int qix = qi[k] - 1;
    if (qix < 0 || (std::size_t) qix >= n_q)
      Rcpp::stop("compare_motifs2_pvalue_cpp: qi out of range");
    pairs_by_query[(std::size_t) qix].push_back((std::size_t) k);
  }

  std::vector<double> pvals(n_pairs, 1.0);

  /* Lambda: build per-column PMF for a given query column. */
  auto build_col_pmf = [&](const double *q_col) -> vec_num_t {
    if (null_mode == 0) return build_col_pmf_empirical(q_col, target_cols_flat);
    else                return build_col_pmf_parametric(q_col, enum_cols);
  };

  /* Lambda: build full per-(start, L) convolution cache for an oriented query. */
  auto build_cache = [&](const list_num_t &q_oriented)
      -> std::vector<std::vector<vec_num_t>> {
    std::size_t wq = q_oriented.size();
    std::vector<vec_num_t> col_pmfs(wq);
    for (std::size_t i = 0; i < wq; ++i) col_pmfs[i] = build_col_pmf(&q_oriented[i][0]);

    std::vector<std::vector<vec_num_t>> cache(wq);
    for (std::size_t start = 0; start < wq; ++start) {
      std::size_t max_L = wq - start;
      cache[start].resize(max_L + 1);   /* [0] unused; L in [1, max_L] */
      cache[start][1] = col_pmfs[start];
      for (std::size_t L = 2; L <= max_L; ++L) {
        cache[start][L] = convolve_with_col(cache[start][L - 1],
                                            col_pmfs[start + L - 1]);
      }
    }
    return cache;
  };

  RcppThread::parallelFor(0, n_q,
    [&](std::size_t q_idx) {
      // TODO: per-query cache build is O(wq^2) -- candidate site for RcppThread::isInterrupted() early-return at the natural boundaries.
      const auto &pairs = pairs_by_query[q_idx];
      if (pairs.empty()) return;

      /* Are any pairs on RC strand? Only build the RC cache if needed. */
      bool need_rc = false;
      for (auto k : pairs) if (strand[(R_xlen_t) k] == 1) { need_rc = true; break; }

      auto cache_fwd = build_cache(queries[q_idx]);
      std::vector<std::vector<vec_num_t>> cache_rc;
      if (need_rc) cache_rc = build_cache(queries_rc[q_idx]);

      std::size_t wq = queries[q_idx].size();

      for (auto k : pairs) {
        int L  = overlap[(R_xlen_t) k];
        int qs = qstart[(R_xlen_t) k];
        int st = strand[(R_xlen_t) k];
        double s = score[(R_xlen_t) k];
        int nt = n_tested[(R_xlen_t) k];

        double p = 1.0;
        if (L > 0 && qs >= 0 && (std::size_t) qs < wq && nt > 0) {
          const auto &cache = (st == 0) ? cache_fwd : cache_rc;
          if (!cache.empty() && (std::size_t) qs < cache.size() &&
              (std::size_t) L < cache[(std::size_t) qs].size()) {
            const vec_num_t &pmf_L = cache[(std::size_t) qs][(std::size_t) L];
            int obs_bin = sum_score_to_bin_pcc(s, L);
            double p_align = tail_prob(pmf_L, obs_bin);
            p = p_align * (double) nt;
            if (p > 1.0) p = 1.0;
            if (p < 0.0) p = 0.0;
          }
        }
        pvals[k] = p;
      }
    }, nthreads);

  return Rcpp::wrap(pvals);
}

/* ============================================================
 * EXPORT 3 (helper): IUPAC consensus of a column range of a motif.
 * Used by R wrapper to attach query/target consensus strings to
 * long-format output. Trivial, but avoids re-implementing in R.
 * ============================================================ */

// [[Rcpp::export(rng = false)]]
Rcpp::CharacterVector compare_motifs2_consensus_cpp(
    const Rcpp::List &mots,
    const Rcpp::IntegerVector &mot_i,    /* 1-based motif index per row */
    const Rcpp::IntegerVector &start,    /* 0-based start column */
    const Rcpp::IntegerVector &len) {

  R_xlen_t n = mot_i.size();
  std::vector<list_num_t> motifs = extract_motif_list(mots);

  Rcpp::CharacterVector out(n);
  for (R_xlen_t k = 0; k < n; ++k) {
    int mix = mot_i[k] - 1;
    int s   = start[k];
    int L   = len[k];
    if (mix < 0 || (std::size_t) mix >= motifs.size() || L <= 0) {
      out[k] = "";
      continue;
    }
    const list_num_t &mot = motifs[(std::size_t) mix];
    if (s < 0) s = 0;
    if ((std::size_t)(s + L) > mot.size()) L = (int) mot.size() - s;
    if (L <= 0) { out[k] = ""; continue; }
    out[k] = consensus_range(mot, s, L);
  }
  return out;
}
