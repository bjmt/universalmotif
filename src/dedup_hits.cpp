#include <Rcpp.h>
#include <algorithm>
#include <utility>
#include <vector>
#include "types.h"

/* dedup_hits_cpp
 * ---------------
 * yamtk-style greedy-by-priority deduplication of motif hits.
 *
 * Within each (seq_i, motif_i, strand_i) group -- with motif_i and/or
 * strand_i optionally ignored -- hits are clustered via a single linear
 * sweep (a new cluster starts when a hit's start strictly exceeds the
 * running max-end of the current cluster). Inside each cluster, hits are
 * ordered by priority (descending, ties broken by the original input
 * order: lower index wins). The first hit is always kept; each subsequent
 * hit is kept iff it does not overlap any already-kept hit (inclusive
 * coordinates).
 *
 * This matches the semantics of `yamtk dedup` (yamtk/src/yamdedup.c).
 *
 * Complexity:
 *   - O(N log N) overall: dominated by the sort step
 *   - O(K log K) per cluster (K = cluster size); K = 1 for most clusters
 *   - O(1) extra memory beyond the index permutation
 *
 * Caller responsibilities:
 *   - `priority` must already be direction-normalised so HIGHER = BETTER.
 *     (e.g. for p-values, the R wrapper passes `-pvalue`; for scores it
 *     passes the score as-is.)
 *   - If `ignore_strand`, the values in `strand_i` are not consulted.
 *   - If `ignore_motif`,  the values in `motif_i`  are not consulted.
 *
 * Returns a LogicalVector of length N: TRUE for kept hits, FALSE for
 * dropped.
 */

// [[Rcpp::export(rng = false)]]
Rcpp::LogicalVector dedup_hits_cpp(const Rcpp::IntegerVector &seq_i,
                                   const Rcpp::IntegerVector &motif_i,
                                   const Rcpp::IntegerVector &start,
                                   const Rcpp::IntegerVector &end,
                                   const Rcpp::IntegerVector &strand_i,
                                   const Rcpp::NumericVector &priority,
                                   const bool ignore_strand = false,
                                   const bool ignore_motif  = false) {

  const R_xlen_t N = seq_i.size();
  Rcpp::LogicalVector keep(N);  // initialised to FALSE
  if (N == 0) return keep;

  if (motif_i.size() != N || start.size() != N || end.size() != N ||
      strand_i.size() != N || priority.size() != N) {
    Rcpp::stop("dedup_hits_cpp: all input vectors must have the same length");
  }

  // Effective key components (collapsed when ignored).
  auto eff_motif  = [&](R_xlen_t i) -> int {
    return ignore_motif  ? 0 : motif_i[i];
  };
  auto eff_strand = [&](R_xlen_t i) -> int {
    return ignore_strand ? 0 : strand_i[i];
  };

  // 1) Sort indices by (seq_i, eff_motif, eff_strand, start).
  std::vector<R_xlen_t> ord(N);
  for (R_xlen_t i = 0; i < N; ++i) ord[i] = i;
  std::sort(ord.begin(), ord.end(),
    [&](R_xlen_t a, R_xlen_t b) {
      if (seq_i[a] != seq_i[b])       return seq_i[a] < seq_i[b];
      int ma = eff_motif(a), mb = eff_motif(b);
      if (ma != mb)                   return ma < mb;
      int sa = eff_strand(a), sb = eff_strand(b);
      if (sa != sb)                   return sa < sb;
      if (start[a] != start[b])       return start[a] < start[b];
      return a < b;  // stable on (seq, motif, strand, start) -> original order
    });

  // Flush a cluster [c_lo, c_hi) of `ord` to the keep mask.
  auto flush_cluster = [&](R_xlen_t c_lo, R_xlen_t c_hi) {
    R_xlen_t n = c_hi - c_lo;
    if (n == 1) {
      keep[ord[c_lo]] = true;
      return;
    }
    // Sort cluster indices by priority desc, ties broken by original index asc.
    std::vector<R_xlen_t> cluster(ord.begin() + c_lo, ord.begin() + c_hi);
    std::sort(cluster.begin(), cluster.end(),
      [&](R_xlen_t a, R_xlen_t b) {
        if (priority[a] != priority[b]) return priority[a] > priority[b];
        return a < b;
      });
    // Greedy: keep each hit iff it does not overlap any already-kept hit.
    std::vector<std::pair<int,int>> kept;
    kept.reserve(n);
    for (R_xlen_t idx : cluster) {
      int s = start[idx], e = end[idx];
      bool collides = false;
      for (const auto &r : kept) {
        if (s <= r.second && r.first <= e) {  // inclusive overlap
          collides = true;
          break;
        }
      }
      if (!collides) {
        keep[idx] = true;
        kept.emplace_back(s, e);
      }
    }
  };

  // 2) Sweep over groups; within each group, sweep clusters by max-end.
  R_xlen_t i = 0;
  while (i < N) {
    const R_xlen_t gi = i;
    const int g_seq    = seq_i[ord[gi]];
    const int g_motif  = eff_motif(ord[gi]);
    const int g_strand = eff_strand(ord[gi]);

    // Locate end of group.
    R_xlen_t gj = gi + 1;
    while (gj < N &&
           seq_i[ord[gj]]    == g_seq &&
           eff_motif(ord[gj])  == g_motif &&
           eff_strand(ord[gj]) == g_strand) {
      ++gj;
    }

    // Within [gi, gj), find clusters by running max-end.
    R_xlen_t cluster_lo = gi;
    int cluster_max_end = end[ord[gi]];
    for (R_xlen_t k = gi + 1; k < gj; ++k) {
      if (start[ord[k]] > cluster_max_end) {
        flush_cluster(cluster_lo, k);
        cluster_lo = k;
        cluster_max_end = end[ord[k]];
      } else if (end[ord[k]] > cluster_max_end) {
        cluster_max_end = end[ord[k]];
      }
    }
    flush_cluster(cluster_lo, gj);

    i = gj;
  }

  return keep;
}
