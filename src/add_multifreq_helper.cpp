#include <Rcpp.h>
#include <unordered_map>
#include "types.h"
#include "shuffle_sequences_helper.h"

std::vector<std::string> single_to_k(const std::string &seq1,
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

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix add_multi_cpp(const std::vector<std::string> &seqs,
    const int k, const std::vector<std::string> &alph) {

  std::size_t seqlen = seqs[0].size(), seqnum = seqs.size();
  if (int(seqlen) < k - 1)
    Rcpp::stop("motif is not long enough");

  list_str_t seqs_k(seqnum);
  for (std::size_t i = 0; i < seqnum; ++i) {
    seqs_k[i] = single_to_k(seqs[i], k);
  }

  vec_str_t klets = get_klet_strings(alph, k);

  std::unordered_map<std::string, int> mklets;
  mklets.reserve(klets.size());
  for (std::size_t i = 0; i < klets.size(); ++i) {
    mklets[klets[i]] = int(i);
  }

  Rcpp::NumericMatrix out(klets.size(), int(seqlen) - k + 1);
  Rcpp::rownames(out) = Rcpp::wrap(klets);

  vec_int_t rowsums(out.ncol(), 0);

  for (std::size_t i = 0; i < seqs_k.size(); ++i) {
    for (std::size_t j = 0; j < seqs_k[0].size(); ++j) {
      out(mklets[seqs_k[i][j]], j) += 1;
      ++rowsums[j];
    }
  }

  for (R_xlen_t i = 0; i < out.ncol(); ++i) {
    for (R_xlen_t j = 0; j < out.nrow(); ++j) {
      out(j, i) /= rowsums[i];
    }
  }

  return out;

}
