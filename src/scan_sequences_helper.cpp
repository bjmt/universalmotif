#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double score_seq(IntegerVector tmp_seq, NumericMatrix score_mat) {
  double score = 0;
  for (int i = 0; i < tmp_seq.size(); ++i) {
    score += score_mat(tmp_seq(i), i);
  }
  return score;
}

// [[Rcpp::export]]
NumericVector scan_1st_order(IntegerVector sequence, NumericMatrix score_mat,
    double min_score) {

  NumericVector to_keep(sequence.size());
  for (int i = 0; i < sequence.size(); ++i) {
    to_keep[i] = 0;
  }

  double tmp_score;
  int max_step = sequence.size() - score_mat.ncol() + 1;
  IntegerVector tmp_seq(score_mat.ncol());

  for (int i = 0; i < max_step; ++i) {

    for (int j = 0; j < score_mat.ncol(); ++j) {
      tmp_seq[j] = sequence(i + j);
    }

    tmp_score = score_seq(tmp_seq, score_mat);
    if (tmp_score >= min_score) to_keep[i] = 1;

  }

  return to_keep;

}

// [[Rcpp::export]]
IntegerVector DNA_to_int_di(StringVector seqs) {

  IntegerVector out(seqs.length() / 2);

  for (int i = 0; i < seqs.length(); ++i) {
    if (i % 2 == 0) {

         if (as<std::string>(seqs(i)) == "A" &&
        as<std::string>(seqs(i + 1)) == "A") out[i / 2] = 0;
    else if (as<std::string>(seqs(i)) == "A" &&
        as<std::string>(seqs(i + 1)) == "C") out[i / 2] = 1;
    else if (as<std::string>(seqs(i)) == "A" &&
        as<std::string>(seqs(i + 1)) == "G") out[i / 2] = 2;
    else if (as<std::string>(seqs(i)) == "A" &&
        as<std::string>(seqs(i + 1)) == "T") out[i / 2] = 3;
    else if (as<std::string>(seqs(i)) == "C" &&
        as<std::string>(seqs(i + 1)) == "A") out[i / 2] = 4;
    else if (as<std::string>(seqs(i)) == "C" &&
        as<std::string>(seqs(i + 1)) == "C") out[i / 2] = 5;
    else if (as<std::string>(seqs(i)) == "C" &&
        as<std::string>(seqs(i + 1)) == "G") out[i / 2] = 6;
    else if (as<std::string>(seqs(i)) == "C" &&
        as<std::string>(seqs(i + 1)) == "T") out[i / 2] = 7;
    else if (as<std::string>(seqs(i)) == "G" &&
        as<std::string>(seqs(i + 1)) == "A") out[i / 2] = 8;
    else if (as<std::string>(seqs(i)) == "G" &&
        as<std::string>(seqs(i + 1)) == "C") out[i / 2] = 9;
    else if (as<std::string>(seqs(i)) == "G" &&
        as<std::string>(seqs(i + 1)) == "G") out[i / 2] = 10;
    else if (as<std::string>(seqs(i)) == "G" &&
        as<std::string>(seqs(i + 1)) == "T") out[i / 2] = 11;
    else if (as<std::string>(seqs(i)) == "T" &&
        as<std::string>(seqs(i + 1)) == "A") out[i / 2] = 12;
    else if (as<std::string>(seqs(i)) == "T" &&
        as<std::string>(seqs(i + 1)) == "C") out[i / 2] = 13;
    else if (as<std::string>(seqs(i)) == "T" &&
        as<std::string>(seqs(i + 1)) == "G") out[i / 2] = 14;
    else if (as<std::string>(seqs(i)) == "T" &&
        as<std::string>(seqs(i + 1)) == "T") out[i / 2] = 15;

    }
  }

  return out;

}

// [[Rcpp::export]]
IntegerVector DNA_to_int_tri(StringVector seqs) {

  IntegerVector out(seqs.length() / 3);
  StringVector score_di(2);
  IntegerVector score_di_i;

  for (int i = 0; i < seqs.length(); ++i) {
    if (i % 3 == 0) {

      score_di[0] = seqs[i + 1];
      score_di[1] = seqs[i + 2];
      score_di_i = DNA_to_int_di(score_di);
    
      if (as<std::string>(seqs(i)) == "A") {
        out[i / 3] = score_di_i[0]; 
      } else if (as<std::string>(seqs(i)) == "C") {
        out[i / 3] = score_di_i[0] + 16;
      } else if (as<std::string>(seqs(i)) == "G") {
        out[i / 3] = score_di_i[0] + 32;
      } else if (as<std::string>(seqs(i)) == "T") {
        out[i / 3] = score_di_i[0] + 48;
      }
    
    }
  }

  return out;

}

// [[Rcpp::export]]
NumericVector scan_2nd_order(IntegerVector sequence, NumericMatrix score_mat,
    double min_score) {

  NumericVector to_keep(sequence.size());
  for (int i = 0; i < sequence.size(); ++i) {
    to_keep[i] = 0;
  }

  double tmp_score;
  int max_step = sequence.size() - score_mat.ncol() + 1;
  IntegerVector tmp_seq(score_mat.ncol());

  for (int i = 0; i < max_step; ++i) {
    for (int j = 0; j < score_mat.ncol(); ++j) {
      tmp_seq[j] = sequence(i + j);
    }
    tmp_score = score_seq(tmp_seq, score_mat);
    if (tmp_score >= min_score) to_keep[i] = 1;
  }

  return to_keep;

}
