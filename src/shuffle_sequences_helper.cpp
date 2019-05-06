#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
String shuffle_markov_loop(R_xlen_t seq_i_l, R_xlen_t seq_i_r, int k,
    StringVector seqout, const StringVector &lets, const NumericMatrix &trans,
    const StringVector &trans_cols) {

  StringVector prev_k_split;
  R_xlen_t trans_i;
  String prev_k, seq_i;
  IntegerVector prev_i;
  NumericVector curr_prob;

  for (R_xlen_t i = seq_i_l; i < seq_i_r; ++i) {

    prev_i = seq(i - k + 1, i - 1);
    prev_k_split = seqout[prev_i];
    prev_k = collapse(prev_k_split);

    trans_i = trans_cols.findName(prev_k);
    curr_prob = trans(_, trans_i);

    seqout[i] = sample(lets, 1, false, curr_prob)[0];

  }

  return collapse(seqout);

}

// [[Rcpp::export(rng = false)]]
StringVector eulerian_walk_cpp(const StringVector &edgelist,
    const StringVector &firstl, R_xlen_t seqlen, int k, const StringVector &last,
    IntegerVector indices) {

  R_xlen_t next_i;
  StringVector nextl;  // if I declare nextl as String, it's converted to int
  String currentl;     // for some reason during
                       // nextl = edgelist[currentl][indices[currentl]]

  StringVector out(seqlen);
  for (int i = 0; i < k - 1; ++i) {
    out[i] = firstl[i];
  }

  for (R_xlen_t i = k - 2; i < seqlen - 2; ++i) {

    currentl = "";
    for (R_xlen_t j = i - k + 2; j <= i; ++j) {  // TODO: this can probably be
      currentl += out[j];                   //       replaced by int indexing
    }

    nextl = edgelist[currentl][indices[currentl]];

    // indices[currentl] = indices[currentl] + 1;  // fails on windows
    next_i = indices[currentl];
    ++next_i;
    indices[currentl] = next_i;

    out[i + 1] = nextl[0];

  }

  out[out.size() - 1] = last[last.size() - 1];

  return out;

}
