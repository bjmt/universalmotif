#include <Rcpp.h>
#include "types.h"

extern const vec_str_t AMINOACIDS2 {
  "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
  "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
};

extern const Rcpp::StringVector AMINOACIDS {
  "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
  "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
};

extern const Rcpp::StringVector DNA { "A", "C", "G", "T" };

extern const Rcpp::StringVector RNA { "A", "C", "G", "U" };

list_int_t R_to_cpp_motif_no_inf(const Rcpp::IntegerMatrix &motif) {

  list_int_t mat(motif.ncol(), vec_int_t(motif.nrow()));
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    for (R_xlen_t j = 0; j < motif.nrow(); ++j) {
      if (motif(j, i) <= -std::numeric_limits<int>::max()) {
        mat[i][j] = -std::numeric_limits<int>::max();
      } else {
        // mat[i][j] = int(motif(j, i) * 1000);
        mat[i][j] = int(motif(j, i));
      }
    }
  }

  return mat;

}

list_int_t R_to_cpp_motif(const Rcpp::NumericMatrix &motif) {

  list_int_t mat(motif.ncol(), vec_int_t(motif.nrow()));
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    for (R_xlen_t j = 0; j < motif.nrow(); ++j) {
      mat[i][j] = int(motif(j, i) * 1000);
    }
  }

  return mat;

}

list_num_t R_to_cpp_motif_num(const Rcpp::NumericMatrix &motif) {

  list_num_t mat(motif.ncol(), vec_num_t(motif.nrow()));
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    for (R_xlen_t j = 0; j < motif.nrow(); ++j) {
      mat[i][j] = motif(j, i);
    }
  }

  return mat;

}

list_int_t R_to_cpp_motif(const Rcpp::IntegerMatrix &motif) {

  list_int_t mat(motif.ncol(), vec_int_t(motif.nrow()));
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    for (R_xlen_t j = 0; j < motif.nrow(); ++j) {
      mat[i][j] = motif(j, i);
    }
  }

  return mat;

}

Rcpp::NumericMatrix cpp_to_R_motif(const list_int_t &motif) {

  Rcpp::NumericMatrix out(motif[0].size(), motif.size());
  for (std::size_t i = 0; i < motif.size(); ++i) {
    Rcpp::NumericVector tmp = Rcpp::wrap(motif[i]);
    out(Rcpp::_, i) = tmp;
  }

  return out;

}

Rcpp::NumericMatrix cpp_to_R_motif(const list_num_t &motif) {

  Rcpp::NumericMatrix out(motif[0].size(), motif.size());
  for (std::size_t i = 0; i < motif.size(); ++i) {
    Rcpp::NumericVector tmp = Rcpp::wrap(motif[i]);
    out(Rcpp::_, i) = tmp;
  }

  return out;

}
