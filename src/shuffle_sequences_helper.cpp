#include <Rcpp.h>
using namespace Rcpp;

// StringVector sample_string_cpp(StringVector x, int size, bool replace = false,
    // sugar::probs_t prob = R_NilValue) {
  // // about 1.5 times faster than base::sample
  // return sample(x, size, replace, prob);
// }

// [[Rcpp::export(rng = false)]]
CharacterMatrix shuffle_random_loop(int seqs_k_n, int k,
    IntegerVector seqs_k_new_i, CharacterMatrix new_seq,
    CharacterMatrix seqs_k) {

  int j, del1, del2;
  IntegerVector to_del, to_keep;

  for (int i = 0; i < seqs_k_n; ++i) {

    if (seqs_k_new_i.length() == 0) break;

    j = seqs_k_new_i[0];
    new_seq(_, i) = seqs_k(j, _);

    del1 = j - k + 1;
    del2 = j + k - 1;

    to_del = seq(del1, del2);
    to_keep = match(seqs_k_new_i, to_del);

    seqs_k_new_i = seqs_k_new_i[is_na(to_keep)];

  }

  return new_seq;

}

// [[Rcpp::export]]
String shuffle_markov_loop(int seq_i_l, int seq_i_r, int k,
    StringVector seqout, StringVector lets, NumericMatrix trans,
    StringVector trans_cols) {

  StringVector prev_k_split;
  int trans_i;
  String prev_k, seq_i;
  IntegerVector prev_i;
  NumericVector curr_prob;

  for (int i = seq_i_l; i < seq_i_r; ++i) {

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
StringVector collapse_rows_mat(CharacterMatrix seqs_k) {
  // ~10 times as fast as apply(seqs_k_matrix, 1, collapse_cpp)
  // ~15 faster than collapse_rows_df(seqs_k_df)

  StringVector out(seqs_k.nrow());

  for (int i = 0; i < seqs_k.nrow(); ++i) {
    out[i] = collapse(seqs_k(i, _));
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
StringVector collapse_rows_df(DataFrame seqs_k) {
  // >2 times as fast as apply(seqs_k_df, 1, collapse_cpp)

  CharacterMatrix seqs_k_mat(seqs_k.nrow(), seqs_k.size());

  for (int i = 0; i < seqs_k.size(); ++i)
    seqs_k_mat(_, i) = StringVector(seqs_k[i]);

  return collapse_rows_mat(seqs_k_mat);

}

// [[Rcpp::export(rng = false)]]
String collapse_cpp(StringVector x) {
  // collapse_cpp(x) is about 3 times faster than base::paste(x, collapse = "")
  // collapse_cpp(c(x, y)) about 2 times faster than base::paste0(x, y)
  return collapse(x);
}
