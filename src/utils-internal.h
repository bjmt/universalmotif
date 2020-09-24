#ifndef _UTILS_INTERNAL_
#define _UTILS_INTERNAL_

#include <Rcpp.h>
#include "types.h"

/* each entry is multiplied by 1000 and converted to an integer */
list_int_t R_to_cpp_motif(const Rcpp::NumericMatrix &motif);

list_int_t R_to_cpp_motif(const Rcpp::IntegerMatrix &motif);

list_int_t R_to_cpp_motif_allow_inf(const Rcpp::NumericMatrix &motif);

list_int_t R_to_cpp_motif_no_inf(const Rcpp::IntegerMatrix &motif);

list_num_t R_to_cpp_motif_num(const Rcpp::NumericMatrix &motif);

Rcpp::NumericMatrix cpp_to_R_motif(const list_int_t &motif);

Rcpp::NumericMatrix cpp_to_R_motif(const list_num_t &motif);

void print_motif(const list_int_t &motif);

void print_motif(const list_num_t &motif);

extern const Rcpp::StringVector AMINOACIDS;

extern const vec_str_t AMINOACIDS2;

extern const Rcpp::StringVector DNA;

extern const Rcpp::StringVector RNA;

#endif
