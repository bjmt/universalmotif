#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
NumericMatrix trim_motif_internal(const NumericMatrix &motif,
    const NumericVector &ic_scores, double min_ic) {

  R_xlen_t n_row = motif.nrow();
  R_xlen_t n = motif.ncol();

  R_xlen_t new_cols = n;
  R_xlen_t cut_left = 0;

  for (R_xlen_t i = 0; i < n; ++i) {
    if (ic_scores[i] < min_ic) {
      --new_cols;
      ++cut_left;
    } else break;
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    if (ic_scores[n - i - 1] < min_ic) {
      --new_cols;
    } else break;
  }

  if (new_cols <= 0) {
    NumericMatrix out;
    return out;
  }

  NumericMatrix out(n_row, new_cols);
  for (R_xlen_t i = 0; i < new_cols; ++i) {
    out(_, i) = motif(_, i + cut_left);
  }

  CharacterVector rown = rownames(motif);
  rownames(out) = rown;

  return out;

}
