#include <Rcpp.h>
#include "types.h"

list_int_t R_to_cpp_motif(const Rcpp::NumericMatrix &motif) {

  list_int_t mat(motif.ncol());
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    mat[i].reserve(motif.nrow());
    for (R_xlen_t j = 0; j < motif.nrow(); ++j) {
      mat[i].push_back(int(motif(j, i) * 1000));
    }
  }

  return mat;

}

list_int_t R_to_cpp_motif(const Rcpp::IntegerMatrix &motif) {

  list_int_t mat(motif.ncol());
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    mat[i].reserve(motif.nrow());
    for (R_xlen_t j = 0; j < motif.nrow(); ++j) {
      mat[i].push_back(motif(j, i));
    }
  }

  return mat;

}
