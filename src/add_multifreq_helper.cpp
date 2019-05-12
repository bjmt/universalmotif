#include <Rcpp.h>
#include "types.h"

// [[Rcpp::export(rng = false)]]
std::vector<std::string> single_to_k(const std::vector<std::string> &seq1,
    const int &k) {

  vec_str_t out(seq1.size() - k + 1, "");
  std::size_t counter;

  for (int i = 0; i < k; ++i) {
    counter = 0;
    while (counter < out.size()) {
      out[counter] += seq1[counter + i];
      ++counter;
    }
  }

  return out;

}
