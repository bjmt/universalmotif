#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double score_seq(IntegerVector tmp_seq, NumericMatrix score_mat) {
  double score = 0;
  for (int i = 0; i < tmp_seq.length(); ++i) {
    score += score_mat(tmp_seq(i), i);
  }
  return score;
}

// [[Rcpp::export]]
NumericVector scan_seq_internal(IntegerVector sequence, NumericMatrix score_mat,
    double min_score) {

  NumericVector to_keep(sequence.length());

  double tmp_score;
  int max_step = sequence.size() - score_mat.ncol() + 1;

  for (int i = 0; i < max_step; ++i) {

    tmp_score = 0;
    for (int j = 0; j < score_mat.ncol(); ++j) {
      tmp_score += score_mat(sequence(i + j), j);
    }
    if (tmp_score >= min_score) to_keep[i] = 1;

  }

  return to_keep;

}

// [[Rcpp::export]]
IntegerVector LETTER_to_int(IntegerVector seqs, int k, IntegerVector letters) {

  IntegerVector out(seqs.length() / k);
  double out_i;
  double l_;
  double let_length = letters.length();

  for (int i = 0; i < seqs.length(); ++i) {
    if (i % k == 0) {

      for (int l = 0; l < k; ++l) {

        l_ = k - l;
        l_ = pow(let_length, l_);
        out_i = seqs(i + l);

        l_ *= out_i / let_length;

        if (l_ == 0) out[i / k] += out_i;
        else out[i / k] += l_;

      }

    }
  }

  return out;

}

// [[Rcpp::export]]
IntegerVector string_to_factor(StringVector x, StringVector y) {

  StringVector lvls = sort_unique(y);
  IntegerVector out = match(x, lvls);

  out.attr("levels") = as<CharacterVector>(lvls);
  out.attr("class") = "factor";

  return out;

}
