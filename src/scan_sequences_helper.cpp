#include <Rcpp.h>
// // [[Rcpp::plugins(cpp11)]]
// // [[Rcpp::depends(BH, bigmemory)]]
// #include <bigmemory/MatrixAccessor.hpp>
// if using bigmemory: add BH and bigmemory to 'LinkingTo'
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
  for (int i = 0; i < sequence.length(); ++i) {
    to_keep[i] = 0;
  }

  double tmp_score;
  int max_step = sequence.size() - score_mat.ncol() + 1;
  IntegerVector tmp_seq(score_mat.ncol());

  for (int i = 0; i < max_step; ++i) {

    for (int j = 0; j < score_mat.ncol(); ++j) {
      tmp_seq[j] = sequence(i + j);
    }

    tmp_score = 0;
    for (int j = 0; j < tmp_seq.length(); ++j) {
      tmp_score += score_mat(tmp_seq(j), j);
    }
    if (tmp_score >= min_score) to_keep[i] = 1;

  }

  return to_keep;

}

// // [[Rcpp::export]]
// NumericVector scan_seq_internal_bigmem(IntegerVector sequence, SEXP score_mat,
    // double min_score) {
//
  // XPtr<BigMatrix> pMat(score_mat);
  // MatrixAccessor<double> mat(*pMat);
  // int score_mat_ncol = pMat->ncol();
//
  // NumericVector to_keep(sequence.length());
  // for (int i = 0; i < sequence.length(); ++i) {
    // to_keep[i] = 0;
  // }
//
  // double tmp_score;
  // int max_step = sequence.size() - score_mat_ncol + 1;
  // IntegerVector tmp_seq(score_mat_ncol);
//
  // for (int i = 0; i < max_step; ++i) {
    // for (int j = 0; j < score_mat_ncol; ++j) {
      // tmp_seq[j] = sequence(i + j);
    // }
    // tmp_score = 0;
    // for (int j = 0; j < tmp_seq.length(); ++j) {
      // tmp_score += mat[tmp_seq(j)][j];
    // }
    // if (tmp_score >= min_score) to_keep[i] = 1;
  // }
//
  // return to_keep;
//
// }

// [[Rcpp::export]]
IntegerVector DNA_to_int_k(StringVector seqs, int k) {

  IntegerVector out(seqs.length() / k);
  for (int i = 0; i < out.length(); ++i) {
    out[i] = 0;
  }

  IntegerVector out_i;
  int l_;

  for (int i = 0; i < seqs.length(); ++i) {
    if (i % k == 0) {

      for (int j = 0; j < k; ++j) {
             if (as<std::string>(seqs(i + j)) == "A") out_i[j] = 0;
        else if (as<std::string>(seqs(i + j)) == "C") out_i[j] = 1;
        else if (as<std::string>(seqs(i + j)) == "G") out_i[j] = 2;
        else if (as<std::string>(seqs(i + j)) == "T") out_i[j] = 3;
      }

      for (int l = 0; l < k; ++l) {
        l_ = k - l;
        l_ = pow(4, l_);
             if (out_i[l] == 0) l_ = 0;
        else if (out_i[l] == 1) l_ = l_ * 0.25;
        else if (out_i[l] == 2) l_ = l_ * 0.5;
        else if (out_i[l] == 3) l_ = l_ * 0.75;
        if (l_ == 0) out[i / k] += out_i[l];
        else out[i / k] += l_;
      }

    }
  }

  return out;

}
