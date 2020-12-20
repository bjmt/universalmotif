#include <Rcpp.h>

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix trim_motif_internal(const Rcpp::NumericMatrix &motif,
    const Rcpp::NumericVector &ic_scores, double min_ic, const int trim_from) {

  R_xlen_t n_row = motif.nrow();
  R_xlen_t n = motif.ncol();

  R_xlen_t new_cols = n;
  R_xlen_t cut_left = 0;

  if (trim_from == 0 || trim_from == 1) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (ic_scores[i] < min_ic) {
        --new_cols;
        ++cut_left;
      } else break;
    }
  }

  if (trim_from == 0 || trim_from == 2) {
    for (R_xlen_t i = 0; i < n; ++i) {
      if (ic_scores[n - i - 1] < min_ic) {
        --new_cols;
      } else break;
    }
  }

  if (new_cols <= 0) {
    Rcpp::NumericMatrix out;
    return out;
  }

  Rcpp::NumericMatrix out(n_row, new_cols);
  for (R_xlen_t i = 0; i < new_cols; ++i) {
    out(Rcpp::_, i) = motif(Rcpp::_, i + cut_left);
  }

  Rcpp::CharacterVector rown = Rcpp::rownames(motif);
  Rcpp::rownames(out) = rown;

  return out;

}
