#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <string>
#include <vector>
#include "types.h"

/* compare_motifs_lite.cpp
 * --------------------
 * yamtk cmp-aligned motif comparison primitives.
 *
 * This file is a port of the algorithm in yamtk/src/yamcmp.c (Tremblay 2026)
 * adapted to operate on Rcpp inputs in a thread-safe way (pre-extract
 * Rcpp objects in serial, parallelise over queries with RcppThread, wrap
 * in serial).
 *
 * Exports:
 *   compare_motifs_lite_medians_cpp:
 *     the median score of each query column against the null columns,
 *     which complete (TOMTOM-style) scores subtract.
 *   compare_motifs_lite_align_cpp:
 *     given a flat list of motif PPMs and per-pair query/target index
 *     vectors, find the best alignment per pair and return
 *     score/chosen/offset/strand/overlap/n_tested/q_start_oriented.
 *   compare_motifs_lite_pvalue_cpp:
 *     given the same motif list plus per-pair alignment metadata, build
 *     the yamtk-style null PMF (empirical or parametric Dirichlet-
 *     Multinomial), convolve per-(start, L) for each unique query, and
 *     return a p-value per pair.
 *
 * The conventions:
 *   - Motifs are PPMs, position-major (mat[i][j] = probability of base j
 *     at position i). 4-base alphabet only (DNA/RNA, treated as ACGT).
 *   - Columns are scored by their Pearson correlation, centred on a uniform
 *     0.25 background (yamcmp.c:1028-1039), or by their Euclidean distance,
 *     negated so that higher scores are always better.
 *   - Score bins: 201, spanning [-1, 1] for correlations and [-sqrt(2), 0]
 *     for negated distances.
 *   - Parametric grid: K=5 simplex, 56 columns (ENUM_N_COLS); Dirichlet
 *     concentration alpha = 4 * bkg (ENUM_DIRICHLET_N).
 *   - Consensus strings follow the package's motif consensus rule.
 */

static const int    PCC_BINS              = 201;
static const int    ENUM_GRID_K           = 5;
static const int    ENUM_DIRICHLET_N      = 4;
static const int    ENUM_N_COLS           = 56;   /* C(K+3, 3) for K=5 */
static const double TIE_TOLERANCE         = 1e-10;

/* Column scoring methods. */
static const int METHOD_PCC  = 0;
static const int METHOD_EUCL = 1;

/* How alignments are chosen beyond the sum over the overlap. */
static const int RULE_SUM        = 0;
static const int RULE_MEAN       = 1;
static const int RULE_COMPLETE   = 2;
static const int RULE_BACKGROUND = 3;

std::string get_consensusC(std::vector<double> pos, const std::string &alphabet,
    const std::string &type, const double pseudocount);

/* ============================================================
 * Column scores and discretization (yamcmp.c:1028-1079)
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

/* Euclidean distance between two columns, negated. */
static inline double eucl_column(const double *q, const double *t) {
  double s = 0.0;
  for (int i = 0; i < 4; ++i) {
    double d = q[i] - t[i];
    s += d * d;
  }
  return -std::sqrt(s);
}

static inline double col_score(int method, const double *q, const double *t) {
  return method == METHOD_EUCL ? eucl_column(q, t) : pcc_column(q, t);
}

static inline double score_lo(int method) {
  return method == METHOD_EUCL ? -std::sqrt(2.0) : -1.0;
}

static inline double score_hi(int method) {
  return method == METHOD_EUCL ? 0.0 : 1.0;
}

/* Sum of L column scores -> bin in [0, L*(PCC_BINS-1)]. yamcmp.c:1072-1079.
 * With L = 1 this bins a single column score. */
static inline int sum_score_to_bin(double s, int L, int method) {
  double lo = score_lo(method), hi = score_hi(method);
  double per_col = (s / (double)L - lo) / (hi - lo);
  if (per_col < 0.0) per_col = 0.0; else if (per_col > 1.0) per_col = 1.0;
  int bin = (int) std::floor(per_col * (PCC_BINS - 1) * L + 0.5);
  int top = L * (PCC_BINS - 1);
  if (bin < 0) bin = 0;
  if (bin > top) bin = top;
  return bin;
}

/* Consensus string for a column range of a motif, by the same rule as the
 * motif's own consensus. */
static std::string consensus_range(const list_num_t &mot,
                                   int start, int L, const std::string &alph) {
  std::string s;
  for (int j = 0; j < L; ++j) s += get_consensusC(mot[start + j], alph, "PPM", 0.0);
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

/* Empirical null: histogram of the scores of q_col against all target
 * columns. `target_cols_flat` is a flat list of target columns
 * (each is a vec_num_t of length 4). */
static vec_num_t build_col_pmf_empirical(const double *q_col,
                                         const list_num_t &target_cols_flat,
                                         int method) {
  vec_num_t pmf(PCC_BINS, 0.0);
  if (target_cols_flat.empty()) return pmf;
  for (const auto &tc : target_cols_flat) {
    double v = col_score(method, q_col, &tc[0]);
    pmf[sum_score_to_bin(v, 1, method)] += 1.0;
  }
  double inv = 1.0 / (double) target_cols_flat.size();
  for (auto &p : pmf) p *= inv;
  return pmf;
}

/* Parametric null: weighted sum over the K=5 simplex grid. */
static vec_num_t build_col_pmf_parametric(const double *q_col,
                                          const std::vector<enum_col_t> &enum_cols,
                                          int method) {
  vec_num_t pmf(PCC_BINS, 0.0);
  for (const auto &e : enum_cols) {
    double v = col_score(method, q_col, e.col);
    pmf[sum_score_to_bin(v, 1, method)] += e.weight;
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

/* Index of the last bin with a positive probability, or -1. */
static inline int top_bin(const vec_num_t &pmf) {
  int top = (int) pmf.size() - 1;
  while (top >= 0 && pmf[top] <= 0.0) --top;
  return top;
}

/* Tail probability P(S' >= observed_bin) (yamcmp.c:1230-1240) */
static double tail_prob(const vec_num_t &pmf, int observed_bin) {
  int len = (int) pmf.size();
  if (observed_bin < 0) observed_bin = 0;
  /* A score beyond what the null reaches, as near-perfect alignments can
   * be against the coarse parametric grid, gets the probability of the
   * null's highest score rather than 0. */
  int top = top_bin(pmf);
  if (top < 0) return 0.0;
  if (observed_bin > top) observed_bin = top;
  double s = 0.0;
  for (int i = observed_bin; i < len; ++i) s += pmf[i];
  if (s < 0.0) s = 0.0;
  if (s > 1.0) s = 1.0;
  return s;
}

/* Running sums of per-column values, so that the sum over any run of
 * columns is a difference of two entries. */
static vec_num_t running_sums(const vec_num_t &v) {
  vec_num_t out(v.size() + 1, 0.0);
  for (std::size_t i = 0; i < v.size(); ++i) out[i + 1] = out[i] + v[i];
  return out;
}

/* The sum of a motif's per-column values outside the `len` columns from
 * `start`, given their running sums. */
static inline double sum_outside(const vec_num_t &sums, int start, int len) {
  return sums[(std::size_t) start] +
    (sums[sums.size() - 1] - sums[(std::size_t) (start + len)]);
}

/* ============================================================
 * Alignment scoring (yamcmp.c:1247-1329)
 * ============================================================ */

struct align_result_t {
  double score;              /* sum of column scores over the overlap */
  double chosen;             /* the value the alignment was chosen by */
  int    L;
  int    offset;             /* d_in_q: query column corresponding to target col 0 */
  int    strand;             /* 0 = '+', 1 = '-' (RC query) */
  int    q_start_oriented;   /* start col in the oriented query (for PMF lookup) */
  int    q_start_in_orig;    /* start col in the original (forward) query (for output) */
  int    t_start;
  int    n_tested;
};

/* What an alignment is chosen by: the sum of column scores over the
 * overlap, its mean, the sum less the query's column medians over the
 * overlap (complete scores), or the sum over every column of both motifs,
 * each overlapping pair counted from both sides and overhanging columns
 * scored against the background. The vectors are running sums, of the
 * query's medians or background scores in each orientation and of the
 * target's background scores. */
struct rule_t {
  int kind;
  const vec_num_t *q_plus;
  const vec_num_t *q_minus;
  const vec_num_t *t_plus;
};

static inline double rule_value(const rule_t &r, double sum, int strand,
                                int qs, int ts, int L) {
  switch (r.kind) {
    case RULE_MEAN:
      return sum / (double) L;
    case RULE_COMPLETE: {
      const vec_num_t &b = strand == 0 ? *r.q_plus : *r.q_minus;
      return sum - (b[(std::size_t) (qs + L)] - b[(std::size_t) qs]);
    }
    case RULE_BACKGROUND: {
      const vec_num_t &b = strand == 0 ? *r.q_plus : *r.q_minus;
      return 2.0 * sum + sum_outside(b, qs, L) + sum_outside(*r.t_plus, ts, L);
    }
    default:
      return sum;
  }
}

static double align_score(const list_num_t &q, int wq,
                          const list_num_t &t, int wt,
                          int d_in_q, int method,
                          int *out_L, int *out_q_start, int *out_t_start) {
  int t_start = (d_in_q >= 0) ? 0 : -d_in_q;
  int t_end   = (wq - d_in_q < wt) ? (wq - d_in_q) : wt;
  int L = t_end - t_start;
  if (L <= 0) { *out_L = 0; return 0.0; }
  int q_start = d_in_q + t_start;
  double s = 0.0;
  for (int j = 0; j < L; ++j) {
    s += col_score(method, &q[q_start + j][0], &t[t_start + j][0]);
  }
  *out_L = L;
  *out_q_start = q_start;
  *out_t_start = t_start;
  return s;
}

static align_result_t best_alignment(const list_num_t &q,
                                     const list_num_t &q_rc,
                                     bool have_rc,
                                     const list_num_t &t,
                                     int min_overlap,
                                     int method,
                                     const rule_t &rule) {
  int wq = (int) q.size();
  int wt = (int) t.size();
  int mo = min_overlap;
  if (mo > wq) mo = wq;
  if (mo > wt) mo = wt;
  int d_lo = mo - wt;
  int d_hi = wq - mo;

  align_result_t best;
  best.score             = -std::numeric_limits<double>::infinity();
  best.chosen            = -std::numeric_limits<double>::infinity();
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
      double s = align_score(qq, wq, t, wt, d, method, &L, &q_start, &t_start);
      if (L < mo || L <= 0) continue;
      double v = rule_value(rule, s, orient, q_start, t_start, L);
      best.n_tested++;
      /* Ties, such as both strands of a palindrome, go to the longer
       * overlap, so that comparing two motifs either way round gives the
       * same mean score, and then to the first alignment, rather than to
       * whichever sum rounds higher. */
      bool tied_longer = std::fabs(v - best.chosen) <= TIE_TOLERANCE && L > best.L;
      bool better = v > best.chosen + TIE_TOLERANCE || tied_longer || best.n_tested == 1;
      if (better) {
        best.score            = s;
        best.chosen           = v;
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

static vec_num_t extract_bkg(const Rcpp::NumericVector &bkg) {
  if (bkg.size() != 4)
    Rcpp::stop("compare_motifs_lite: bkg must have length 4");
  vec_num_t out(4);
  for (int i = 0; i < 4; ++i) out[i] = bkg[i];
  return out;
}

/* Per-motif column medians from R, as running sums in each orientation:
 * the minus orientation reverses the medians of the plus one. Motifs
 * without medians (not used as queries) get empty vectors. */
static void median_sums(const Rcpp::List &medians, std::size_t n_mot,
                        std::vector<vec_num_t> &plus,
                        std::vector<vec_num_t> &minus) {
  plus.assign(n_mot, vec_num_t());
  minus.assign(n_mot, vec_num_t());
  if ((std::size_t) medians.size() != n_mot)
    Rcpp::stop("compare_motifs_lite: need one medians vector per motif");
  for (std::size_t i = 0; i < n_mot; ++i) {
    Rcpp::NumericVector m = medians[(R_xlen_t) i];
    vec_num_t v(m.begin(), m.end());
    plus[i] = running_sums(v);
    std::reverse(v.begin(), v.end());
    minus[i] = running_sums(v);
  }
}

/* Running sums of each motif's column scores against the background, in
 * each orientation. */
static void background_sums(const std::vector<list_num_t> &motifs,
                            const vec_num_t &bkg, int method, bool rc,
                            std::vector<vec_num_t> &plus,
                            std::vector<vec_num_t> &minus) {
  std::size_t n = motifs.size();
  plus.assign(n, vec_num_t());
  minus.assign(n, vec_num_t());
  for (std::size_t i = 0; i < n; ++i) {
    vec_num_t s(motifs[i].size());
    for (std::size_t j = 0; j < s.size(); ++j)
      s[j] = col_score(method, &motifs[i][j][0], &bkg[0]);
    plus[i] = running_sums(s);
    if (rc) {
      list_num_t r = rc_motif(motifs[i]);
      for (std::size_t j = 0; j < s.size(); ++j)
        s[j] = col_score(method, &r[j][0], &bkg[0]);
      minus[i] = running_sums(s);
    }
  }
}

/* ============================================================
 * EXPORT 0: column medians for complete scores.
 *
 * The median score of each column of each wanted motif against the null
 * columns: those of all the motifs, in both orientations with RC, or the
 * parametric grid with its weights.
 * ============================================================ */

// [[Rcpp::export(rng = false)]]
Rcpp::List compare_motifs_lite_medians_cpp(
    const Rcpp::List &mots,
    const Rcpp::LogicalVector &wanted,
    const bool RC,
    const int null_mode,
    const Rcpp::NumericVector &bkg,
    const int method,
    const int nthreads = 1) {

  std::vector<list_num_t> motifs = extract_motif_list(mots);
  std::size_t n_mot = motifs.size();
  if ((std::size_t) wanted.size() != n_mot)
    Rcpp::stop("compare_motifs_lite_medians_cpp: need one `wanted` per motif");
  std::vector<int> want(n_mot);
  for (std::size_t i = 0; i < n_mot; ++i) want[i] = wanted[(R_xlen_t) i] == TRUE;
  vec_num_t bkg_v = extract_bkg(bkg);

  list_num_t null_cols;
  vec_num_t null_w;
  if (null_mode == 1) {
    for (const auto &e : build_enum_cols(bkg_v)) {
      null_cols.push_back(vec_num_t(e.col, e.col + 4));
      null_w.push_back(e.weight);
    }
  } else {
    null_cols = flatten_target_cols(motifs);
    if (RC) {
      for (const auto &m : motifs)
        for (const auto &col : rc_motif(m)) null_cols.push_back(col);
    }
    null_w.assign(null_cols.size(), 1.0);
  }
  double half = 0.0;
  for (double w : null_w) half += w;
  half /= 2.0;

  std::vector<vec_num_t> out(n_mot);
  RcppThread::parallelFor(0, n_mot,
    [&](std::size_t i) {
      if (!want[i]) return;
      out[i].assign(motifs[i].size(), 0.0);
      std::vector<std::pair<double, double>> v(null_cols.size());
      for (std::size_t j = 0; j < motifs[i].size(); ++j) {
        for (std::size_t k = 0; k < null_cols.size(); ++k)
          v[k] = std::make_pair(col_score(method, &motifs[i][j][0], &null_cols[k][0]),
                                null_w[k]);
        std::sort(v.begin(), v.end(),
                  [](const std::pair<double, double> &a,
                     const std::pair<double, double> &b) { return a.first < b.first; });
        double cum = 0.0, med = 0.0;
        for (const auto &x : v) {
          cum += x.second;
          if (cum >= half) { med = x.first; break; }
        }
        out[i][j] = med;
      }
    }, nthreads);

  Rcpp::List res((R_xlen_t) n_mot);
  for (std::size_t i = 0; i < n_mot; ++i) res[(R_xlen_t) i] = Rcpp::wrap(out[i]);
  return res;
}

/* ============================================================
 * EXPORT 1: per-pair best alignment, no p-value.
 *
 * method: 0 = PCC, 1 = EUCL. overhang: 0 = ignore, 1 = background.
 * With `complete`, `medians` holds each query's column medians.
 * ============================================================ */

// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame compare_motifs_lite_align_cpp(
    const Rcpp::List &mots,
    const Rcpp::IntegerVector &qi,
    const Rcpp::IntegerVector &ti,
    const int min_overlap,
    const bool RC,
    const int nthreads,
    const int method,
    const int overhang,
    const bool complete,
    const Rcpp::List &medians,
    const Rcpp::NumericVector &bkg) {

  R_xlen_t n_pairs = qi.size();
  if (ti.size() != n_pairs)
    Rcpp::stop("compare_motifs_lite_align_cpp: qi and ti must have equal length");

  /* Serial pre-extract */
  std::vector<list_num_t> motifs = extract_motif_list(mots);
  std::size_t n_mot = motifs.size();
  std::vector<list_num_t> motifs_rc(n_mot);
  if (RC) {
    for (std::size_t i = 0; i < n_mot; ++i) motifs_rc[i] = rc_motif(motifs[i]);
  }

  std::vector<vec_num_t> extra_plus, extra_minus;
  int kind = method == METHOD_EUCL ? RULE_MEAN : RULE_SUM;
  if (complete) {
    kind = RULE_COMPLETE;
    median_sums(medians, n_mot, extra_plus, extra_minus);
  } else if (overhang == 1) {
    kind = RULE_BACKGROUND;
    background_sums(motifs, extract_bkg(bkg), method, RC, extra_plus, extra_minus);
  }

  /* Convert 1-based -> 0-based pair indices */
  std::vector<int> qi_0(n_pairs), ti_0(n_pairs);
  for (R_xlen_t k = 0; k < n_pairs; ++k) {
    qi_0[(std::size_t) k] = qi[k] - 1;
    ti_0[(std::size_t) k] = ti[k] - 1;
    if (qi_0[(std::size_t) k] < 0 || (std::size_t) qi_0[(std::size_t) k] >= n_mot ||
        ti_0[(std::size_t) k] < 0 || (std::size_t) ti_0[(std::size_t) k] >= n_mot)
      Rcpp::stop("compare_motifs_lite_align_cpp: motif index out of range");
    if (kind == RULE_COMPLETE &&
        extra_plus[(std::size_t) qi_0[(std::size_t) k]].size() !=
          motifs[(std::size_t) qi_0[(std::size_t) k]].size() + 1)
      Rcpp::stop("compare_motifs_lite_align_cpp: missing column medians for a query");
  }

  std::vector<double> out_score(n_pairs), out_chosen(n_pairs);
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
      rule_t rule;
      rule.kind = kind;
      rule.q_plus = rule.q_minus = rule.t_plus = nullptr;
      if (kind == RULE_COMPLETE || kind == RULE_BACKGROUND) {
        rule.q_plus  = &extra_plus[(std::size_t) qix];
        rule.q_minus = &extra_minus[(std::size_t) qix];
        rule.t_plus  = &extra_plus[(std::size_t) tix];
      }
      align_result_t a = best_alignment(q, qr, RC, t, min_overlap, method, rule);
      out_score[k]   = a.score;
      out_chosen[k]  = a.chosen;
      out_L[k]       = a.L;
      out_offset[k]  = a.offset;
      out_strand[k]  = a.strand;
      out_qstart[k]  = a.q_start_oriented;
      out_ntested[k] = a.n_tested;
    }, nthreads);

  return Rcpp::DataFrame::create(
    Rcpp::_["score"]            = out_score,
    Rcpp::_["chosen"]           = out_chosen,
    Rcpp::_["overlap"]          = out_L,
    Rcpp::_["offset"]           = out_offset,
    Rcpp::_["strand"]           = out_strand,
    Rcpp::_["q_start_oriented"] = out_qstart,
    Rcpp::_["n_tested"]         = out_ntested,
    Rcpp::_["stringsAsFactors"] = false
  );
}

/* ============================================================
 * P-value helpers for complete scores and background overhangs.
 * ============================================================ */

/* Null distributions of an oriented query's column sums, by first column
 * and number of columns. */
typedef std::vector<std::vector<vec_num_t>> null_cache_t;

/* The probability that the best of all the alignments of a target of width
 * `wt` to a query of width `wq` reaches a complete score of `score`,
 * treating the alignments as independent. Each orientation tried brings
 * its null distributions and the running sums of its column medians. */
static double best_of_all(double score, int wq, int wt, int min_overlap,
                          int method,
                          const std::vector<const null_cache_t *> &caches,
                          const std::vector<const vec_num_t *> &medians) {
  int mo = min_overlap;
  if (mo > wq) mo = wq;
  if (mo > wt) mo = wt;
  double hi = score_hi(method);
  /* The log of the probability that no alignment reaches the score. */
  double log_none = 0.0;
  for (std::size_t o = 0; o < caches.size(); ++o) {
    const null_cache_t &cache = *caches[o];
    const vec_num_t &med = *medians[o];
    for (int d = mo - wt; d <= wq - mo; ++d) {
      int t_start = d >= 0 ? 0 : -d;
      int len = std::min(wq - d, wt) - t_start;
      if (len < mo || len <= 0) continue;
      int start = d + t_start;
      /* The complete score is the overlap's sum less its medians, so
       * reaching it takes an overlap sum of the score plus the medians. */
      double raw = score + med[(std::size_t) (start + len)] - med[(std::size_t) start];
      /* A sum of column scores cannot exceed the best each column can
       * score. */
      if (raw > (double) len * hi + 1e-9) continue;
      const vec_num_t &pmf = cache[(std::size_t) start][(std::size_t) len];
      int bin = sum_score_to_bin(raw, len, method);
      double tail = 0.0;
      for (std::size_t i = (std::size_t) bin; i < pmf.size(); ++i) tail += pmf[i];
      if (tail > 1.0) tail = 1.0;
      if (tail >= 1.0) return 1.0;
      log_none += std::log1p(-tail);
    }
  }
  double p = -std::expm1(log_none);
  if (p < 0.0) p = 0.0;
  if (p > 1.0) p = 1.0;
  return p;
}

/* The probability that twice a sum drawn from `overlap`, the null of the
 * overlap's column scores, plus a sum drawn from `outside`, the null of
 * the overhanging target columns' scores against the background, reaches
 * `score`, a sum of `units` column scores. As in tail_prob(), a score
 * beyond what the null can reach gets the probability of the null's
 * highest. */
static double overhang_tail(const vec_num_t &overlap, const vec_num_t &outside,
                            double score, int units, int method) {
  int top_x = top_bin(overlap), top_y = top_bin(outside);
  if (top_x < 0 || top_y < 0) return 0.0;
  int bin = std::min(sum_score_to_bin(score, units, method), 2 * top_x + top_y);
  vec_num_t at_least(overlap.size() + 1, 0.0);
  for (std::size_t i = overlap.size(); i-- > 0; )
    at_least[i] = at_least[i + 1] + overlap[i];
  double p = 0.0;
  for (std::size_t j = 0; j < outside.size(); ++j) {
    double py = outside[j];
    if (py == 0.0) continue;
    /* Twice the overlap's bin must reach what the outside leaves. */
    int rest = bin - (int) j;
    if (rest < 0) rest = 0;
    std::size_t i = (std::size_t) ((rest + 1) / 2);
    if (i < at_least.size()) p += py * at_least[i];
  }
  if (p < 0.0) p = 0.0;
  if (p > 1.0) p = 1.0;
  return p;
}

/* ============================================================
 * EXPORT 2: per-pair p-value (yamtk-style).
 *
 * Inputs (per pair k):
 *   qi[k], ti[k], qstart[k], overlap[k], strand[k], offset[k], score[k],
 *   chosen[k], n_tested[k]
 * Modes:
 *   null_mode = 0 -> empirical from target_mats
 *   null_mode = 1 -> parametric Dirichlet-Multinomial over K=5 simplex
 * Returns a NumericVector of p-values, one per pair: Bonferroni-corrected
 * for the alignments tried, or with `complete` the probability that the
 * best of all of them reaches the score.
 *
 * Parallelisation: over unique queries. Each thread:
 *   - builds per-column PMFs for its query (and its RC if any pair uses '-')
 *   - builds the per-(start, L) convolution cache
 *   - looks up tail probability for each of its pairs
 *   - frees the cache when done
 * ============================================================ */

// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector compare_motifs_lite_pvalue_cpp(
    const Rcpp::List &query_mats,
    const Rcpp::List &target_mats,
    const Rcpp::IntegerVector &qi,
    const Rcpp::IntegerVector &ti,
    const Rcpp::IntegerVector &qstart,
    const Rcpp::IntegerVector &overlap,
    const Rcpp::IntegerVector &strand,
    const Rcpp::IntegerVector &offset,
    const Rcpp::NumericVector &score,
    const Rcpp::NumericVector &chosen,
    const Rcpp::IntegerVector &n_tested,
    const Rcpp::NumericVector &bkg,
    const int null_mode,
    const int nthreads,
    const int method,
    const int overhang,
    const bool complete,
    const Rcpp::List &medians,
    const bool RC,
    const int min_overlap) {

  R_xlen_t n_pairs = qi.size();
  if (ti.size()       != n_pairs ||
      qstart.size()   != n_pairs ||
      overlap.size()  != n_pairs ||
      strand.size()   != n_pairs ||
      offset.size()   != n_pairs ||
      score.size()    != n_pairs ||
      chosen.size()   != n_pairs ||
      n_tested.size() != n_pairs)
    Rcpp::stop("compare_motifs_lite_pvalue_cpp: input vectors must be the same length");
  if (null_mode != 0 && null_mode != 1)
    Rcpp::stop("compare_motifs_lite_pvalue_cpp: null_mode must be 0 (empirical) or 1 (parametric)");

  /* Serial pre-extract */
  std::vector<list_num_t> queries = extract_motif_list(query_mats);
  std::vector<list_num_t> targets = extract_motif_list(target_mats);
  std::size_t n_q = queries.size();

  std::vector<list_num_t> queries_rc(n_q);
  for (std::size_t i = 0; i < n_q; ++i) queries_rc[i] = rc_motif(queries[i]);

  vec_num_t bkg_v = extract_bkg(bkg);

  /* For empirical mode, flatten all target columns once. */
  list_num_t target_cols_flat;
  if (null_mode == 0) target_cols_flat = flatten_target_cols(targets);

  /* For parametric mode, build the enum grid once. */
  std::vector<enum_col_t> enum_cols;
  if (null_mode == 1) enum_cols = build_enum_cols(bkg_v);

  /* Lambda: build per-column PMF for a given query column. */
  auto build_col_pmf = [&](const double *q_col) -> vec_num_t {
    if (null_mode == 0) return build_col_pmf_empirical(q_col, target_cols_flat, method);
    else                return build_col_pmf_parametric(q_col, enum_cols, method);
  };

  std::vector<vec_num_t> med_plus, med_minus, bkg_plus, bkg_minus;
  if (complete) median_sums(medians, n_q, med_plus, med_minus);

  /* With background overhangs, the null distributions of the sums of
   * target columns' scores against the background, by number of columns. */
  std::vector<vec_num_t> overhang_pmfs;
  if (!complete && overhang == 1) {
    background_sums(targets, bkg_v, method, false, bkg_plus, bkg_minus);
    std::size_t widest = 0;
    for (const auto &t : targets) widest = std::max(widest, t.size());
    vec_num_t one = build_col_pmf(&bkg_v[0]);
    overhang_pmfs.push_back(vec_num_t(1, 1.0));
    for (std::size_t k = 1; k <= widest; ++k)
      overhang_pmfs.push_back(convolve_with_col(overhang_pmfs[k - 1], one));
  }

  /* 1-based -> 0-based; group pair indices by query index. */
  std::vector<std::vector<std::size_t>> pairs_by_query(n_q);
  for (R_xlen_t k = 0; k < n_pairs; ++k) {
    int qix = qi[k] - 1;
    if (qix < 0 || (std::size_t) qix >= n_q)
      Rcpp::stop("compare_motifs_lite_pvalue_cpp: qi out of range");
    int tix = ti[k] - 1;
    if (tix < 0 || (std::size_t) tix >= targets.size())
      Rcpp::stop("compare_motifs_lite_pvalue_cpp: ti out of range");
    if (complete && med_plus[(std::size_t) qix].size() != queries[(std::size_t) qix].size() + 1)
      Rcpp::stop("compare_motifs_lite_pvalue_cpp: missing column medians for a query");
    pairs_by_query[(std::size_t) qix].push_back((std::size_t) k);
  }

  std::vector<double> pvals(n_pairs, 1.0);

  /* Lambda: build full per-(start, L) convolution cache for an oriented query. */
  auto build_cache = [&](const list_num_t &q_oriented) -> null_cache_t {
    std::size_t wq = q_oriented.size();
    std::vector<vec_num_t> col_pmfs(wq);
    for (std::size_t i = 0; i < wq; ++i) col_pmfs[i] = build_col_pmf(&q_oriented[i][0]);

    null_cache_t cache(wq);
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

      /* Are any pairs on RC strand? Only build the RC cache if needed;
       * complete scores look at every alignment of both strands. */
      bool need_rc = complete && RC;
      for (auto k : pairs) if (strand[(R_xlen_t) k] == 1) { need_rc = true; break; }

      auto cache_fwd = build_cache(queries[q_idx]);
      null_cache_t cache_rc;
      if (need_rc) cache_rc = build_cache(queries_rc[q_idx]);

      std::size_t wq = queries[q_idx].size();

      std::vector<const null_cache_t *> orient_caches;
      std::vector<const vec_num_t *> orient_meds;
      if (complete) {
        orient_caches.push_back(&cache_fwd);
        orient_meds.push_back(&med_plus[q_idx]);
        if (RC && !cache_rc.empty()) {
          orient_caches.push_back(&cache_rc);
          orient_meds.push_back(&med_minus[q_idx]);
        }
      }

      for (auto k : pairs) {
        int L  = overlap[(R_xlen_t) k];
        int qs = qstart[(R_xlen_t) k];
        int st = strand[(R_xlen_t) k];
        double s = score[(R_xlen_t) k];
        int nt = n_tested[(R_xlen_t) k];
        std::size_t tix = (std::size_t) (ti[(R_xlen_t) k] - 1);

        double p = 1.0;
        if (L > 0 && qs >= 0 && (std::size_t) qs < wq && nt > 0) {
          const auto &cache = (st == 0) ? cache_fwd : cache_rc;
          const vec_num_t *pmf_L = nullptr;
          if (!cache.empty() && (std::size_t) qs < cache.size() &&
              (std::size_t) L < cache[(std::size_t) qs].size())
            pmf_L = &cache[(std::size_t) qs][(std::size_t) L];
          if (complete) {
            /* An observed score beyond what the null reaches keeps the
             * probability of the null's highest. */
            double floor_p = pmf_L == nullptr ? 0.0 :
              tail_prob(*pmf_L, sum_score_to_bin(s, L, method));
            p = best_of_all(chosen[(R_xlen_t) k], (int) wq,
                            (int) targets[tix].size(), min_overlap, method,
                            orient_caches, orient_meds);
            if (floor_p > p) p = floor_p;
          } else if (overhang == 1) {
            if (pmf_L != nullptr) {
              /* The query's overhanging columns score the same whatever
               * the target, so only the target's columns, facing the
               * query or the background, vary. */
              int wt = (int) targets[tix].size();
              int n_out = wt - L;
              int t_start = std::max(0, -offset[(R_xlen_t) k]);
              double t_out = sum_outside(bkg_plus[tix], t_start, L);
              double tail = overhang_tail(*pmf_L, overhang_pmfs[(std::size_t) n_out],
                                          2.0 * s + t_out, 2 * L + n_out, method);
              p = tail * (double) nt;
            }
          } else if (pmf_L != nullptr) {
            int obs_bin = sum_score_to_bin(s, L, method);
            double p_align = tail_prob(*pmf_L, obs_bin);
            p = p_align * (double) nt;
          }
          if (p > 1.0) p = 1.0;
          if (p < 0.0) p = 0.0;
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
Rcpp::CharacterVector compare_motifs_lite_consensus_cpp(
    const Rcpp::List &mots,
    const Rcpp::IntegerVector &mot_i,    /* 1-based motif index per row */
    const Rcpp::IntegerVector &start,    /* 0-based start column */
    const Rcpp::IntegerVector &len,
    const std::string &alph = "DNA") {

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
    out[k] = consensus_range(mot, s, L, alph);
  }
  return out;
}
