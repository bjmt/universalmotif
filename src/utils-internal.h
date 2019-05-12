#ifndef _UTILS_INTERNAL_
#define _UTILS_INTERNAL_

#include <Rcpp.h>
#include "types.h"

/* each entry is multiplied by 1000 and converted to an integer */
list_int_t R_to_cpp_motif(const Rcpp::NumericMatrix &motif);

list_int_t R_to_cpp_motif(const Rcpp::IntegerMatrix &motif);

#endif
