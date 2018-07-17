#include <Rcpp.h>
using namespace Rcpp;

double motif_euclidean(NumericMatrix mot1, NumericMatrix mot2) {

  // Adapted from TFBSTools:::PWMEuclidean

  int mat_size = mot1.nrow() * mot1.ncol();
  NumericMatrix diff_mat(mot1.nrow(), mot1.ncol());

  for (int i = 0; i < mat_size; ++i) {
    diff_mat[i] = mot1[i] - mot2[i];
    diff_mat[i] = pow(diff_mat[i], 2.0);
  }

  NumericVector diff_mat_colsums = colSums(diff_mat);
  for (int i = 0; i < diff_mat_colsums.length(); ++i) {
    diff_mat_colsums[i] = sqrt(diff_mat_colsums[i]);
  }

  double k = sum(diff_mat_colsums);
  return k / sqrt(2.0) / mot1.ncol();  // report smallest distance;
                                       // distance measure
}

double motif_pearson(NumericMatrix mot1, NumericMatrix mot2) {

  // Adapted from TFBSTools:::PWMPearson

  int mat_size = mot1.nrow() * mot1.ncol();
  
  for (int i = 0; i < mat_size; ++i) {
    mot1[i] -= 1.0 / mot1.nrow();
    mot2[i] -= 1.0 / mot1.nrow();
  }

  NumericMatrix mult_mat(mot1.nrow(), mot1.ncol());
  for (int i = 0; i < mat_size; ++i) {
    mult_mat[i] = mot1[i] * mot2[i];
  }

  NumericVector top = colSums(mult_mat);

  for (int i = 0; i < mat_size; ++i) {
    mot1[i] = pow(mot1[i], 2.0);
    mot2[i] = pow(mot2[i], 2.0);
  }

  NumericVector mot1_colsums = colSums(mot1);
  NumericVector mot2_colsums = colSums(mot2);
  NumericVector bottom(mot1.ncol());
  for (int i = 0; i < mot1_colsums.length(); ++i) {
    bottom[i] = mot1_colsums[i] * mot2_colsums[i];
    bottom[i] = sqrt(bottom[i]);
  }

  NumericVector top_bot(mot1.ncol());
  for (int i = 0; i < mot1.ncol(); ++i) {
    top_bot[i] = top[i] / bottom[i];
  }

  return 1.0 / mot1.ncol() * sum(top_bot);  // report largest correlation;
                                            // similarity measure
}

double motif_kl(NumericMatrix mot1, NumericMatrix mot2) {

  // Adapted from TFBSTools:::PWMKL

  int mat_size = mot1.nrow() * mot1.ncol();
  NumericMatrix new_mat(mot1.nrow(), mot1.ncol());

  for (int i = 0; i < mat_size; ++i) {
    new_mat[i] = mot1[i] * log(mot1[i] / mot2[i]);
    new_mat[i] += mot2[i] * log(mot2[i] / mot1[i]);
  }

  NumericVector mat_colsums = colSums(new_mat);

  return 0.5 / mot1.ncol() * sum(mat_colsums);  // report smallest divergence;
                                                // distance measure
}

// [[Rcpp::export]]
double motif_simil_internal(NumericMatrix mot1, NumericMatrix mot2,
    StringVector method) {

  // Adapted from TFBSTools::PWMSimilarity

  NumericVector widths(2);
  widths[0] = mot1.ncol();
  widths[1] = mot2.ncol();
  int min_width = min(widths);
  int for_i = 1 + mot1.ncol() - min_width;
  int for_j = 1 + mot2.ncol() - min_width;
  NumericMatrix mot1_tmp(mot1.nrow(), min_width);
  NumericMatrix mot2_tmp(mot2.nrow(), min_width);
  widths[1] = mot2.ncol() - min_width + 1;
  NumericVector ans(for_i + for_j - 1);

  int method_i = 0;
       if (method[0] == "Euclidean") method_i = 1;
  else if (method[0] == "Pearson")   method_i = 2;
  else if (method[0] == "KL")        method_i = 3;

  for (int i = 0; i < for_i; ++i) {
    for (int j = 0; j < for_j; ++j) {
      
      for (int k = 0; k < min_width; ++k) {
        mot1_tmp(_, k) = mot1(_, k + i);
        mot2_tmp(_, k) = mot2(_, k + j);
      }

      switch(method_i) {

        case 1: ans(i + j) = motif_euclidean(mot1_tmp, mot2_tmp);
                break;
        case 2: ans(i + j) = motif_pearson(mot1_tmp, mot2_tmp);
                break;
        case 3: ans(i + j) = motif_kl(mot1_tmp, mot2_tmp);
                break;

      }

    }
  }

  switch(method_i) {

    case 1: return min(ans);
    case 2: return max(ans);
    case 3: return min(ans);

  }

  return 999.0;

}
