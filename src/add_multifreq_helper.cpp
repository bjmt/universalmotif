#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
StringVector single_to_k(StringVector seq1, int k) {

  // NOTE: switching the assignment of rows and columns makes the function
  //       slightly faster for smaller sequences (<1000), but slower for bigger
  //       sequences.

  R_xlen_t n = seq1.size();
  R_xlen_t n2 = n - k + 1;
  StringMatrix seq_mat(k, n2);
  StringVector out(n2);

  seq_mat(0, _) = seq1[seq(0, n - k)];

  for (int i = 1; i < k; ++i) {
    seq_mat(i, _) = seq1[seq(i, i + n2 - 1)];
  }

  for (R_xlen_t i = 0; i < n2; ++i) {
    out[i] = collapse(seq_mat(_, i));
  }

  return out;

}
