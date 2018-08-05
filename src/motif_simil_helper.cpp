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

  LogicalVector to_keep = !is_na(diff_mat_colsums);
  diff_mat_colsums = diff_mat_colsums[to_keep];

  double k = sum(diff_mat_colsums);
  return k / sqrt(2.0) / diff_mat_colsums.length();
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

  LogicalVector to_keep = !is_na(top_bot);
  top_bot = top_bot[to_keep];

  return 1.0 / top_bot.length() * sum(top_bot);
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

  LogicalVector to_keep = !is_na(mat_colsums);
  mat_colsums = mat_colsums[to_keep];

  return 0.5 / mat_colsums.length() * sum(mat_colsums);
}

// [[Rcpp::export]]
List add_cols(NumericMatrix mot1, NumericMatrix mot2,
    NumericVector ic1, NumericVector ic2, int overlap) {

  List out;

  int ncol_1 = mot1.ncol();
  int nrow_1 = mot1.nrow();
  int total_1 = ncol_1 * nrow_1;

  int ncol_2 = mot2.ncol();
  int nrow_2 = mot2.nrow();
  int total_2 = ncol_2 * nrow_2;

  int ncol_1_toadd;
  if (overlap > ncol_2) ncol_1_toadd = 0;
  else ncol_1_toadd = ncol_2 - overlap;

  int ncol_2_toadd;
  if (overlap > ncol_1) ncol_2_toadd = 0;
  else ncol_2_toadd = ncol_1 - overlap;

  if (ncol_1_toadd == 0 || ncol_2_toadd == 0) {
    out = List::create(mot1, mot2, ic1, ic2);
    return out;
  }

  if (mot2.ncol() >= mot1.ncol()) {

    int total_1_toadd = ncol_1_toadd * nrow_1;
    int total_1_new = total_1_toadd * 2 + total_1;
    NumericMatrix mot1_new(nrow_1, ncol_1_toadd * 2 + ncol_1);
    NumericVector ic1_new(ic1.length() + ncol_1_toadd * 2, NA_REAL);

    for (int i = 0; i < total_1_new; ++i) {
      mot1_new[i] = NA_REAL;
    }
    for (int i = ncol_1_toadd; i < ncol_1_toadd + ncol_1; ++i) {
      mot1_new(_, i) = mot1(_, i - ncol_1_toadd);
    }

    for (int i = ncol_1_toadd; i < ncol_1_toadd + ic1.length(); ++i) {
      ic1_new[i] = ic1[i - ncol_1_toadd];
    }

  out = List::create(mot1_new, mot2, ic1_new, ic2);

  } else {
  
    int total_2_toadd = ncol_2_toadd * nrow_2;
    int total_2_new = total_2_toadd * 2 + total_2;
    NumericMatrix mot2_new(nrow_2, ncol_2_toadd * 2 + ncol_2);
    NumericVector ic2_new(ic2.length() + ncol_2_toadd * 2, NA_REAL);

    for (int i = 0; i < total_2_new; ++i) {
      mot2_new[i] = NA_REAL;
    }
    for (int i = ncol_2_toadd; i < ncol_2_toadd + ncol_2; ++i) {
      mot2_new(_, i) = mot2(_, i - ncol_2_toadd);
    }

    for (int i = ncol_2_toadd; i < ncol_2_toadd + ic2.length(); ++i) {
      ic2_new[i] = ic2[i - ncol_2_toadd];
    }

  out = List::create(mot1, mot2_new, ic1, ic2_new);
  
  }

  return out;

}

NumericMatrix rev_motif(NumericMatrix motif) {

  NumericMatrix out(motif.nrow(), motif.ncol());
  for (int i = 0; i < motif.nrow(); ++i) {
    out(i, _) = motif(motif.nrow() - i - 1, _);
  }
  NumericMatrix out2 = Rcpp::clone(out);
  for (int i = 0; i < motif.ncol(); ++i) {
    out(_, i) = out2(_, motif.ncol() - i - 1);
  }

  return out;

}

double motif_simil_internal(NumericMatrix mot1, NumericMatrix mot2,
    String method, int min_overlap, bool tryRC, NumericVector ic1,
    NumericVector ic2, double min_ic);

double motif_simil_rc(NumericMatrix mot1, NumericMatrix mot2,
    String method, int min_overlap, NumericVector ic1,
    NumericVector ic2, double min_ic) {

  NumericMatrix mot1_rc = rev_motif(mot1);
  double ans = motif_simil_internal(mot1_rc, mot2, method, min_overlap, false,
      ic1, ic2, min_ic);

  return ans;

}

// [[Rcpp::export]]
double motif_simil_internal(NumericMatrix mot1, NumericMatrix mot2,
    String method, int min_overlap, bool tryRC, NumericVector ic1,
    NumericVector ic2, double min_ic) {

  // Adapted from TFBSTools::PWMSimilarity

  if (method == "KL") {
    mot1 = mot1 + 0.0001;
    mot2 = mot2 + 0.0001;
  }

  NumericVector ans_rc;
  if (tryRC) ans_rc = motif_simil_rc(mot1, mot2, method, min_overlap,
      ic1, ic2, min_ic);
  
  List new_mots = add_cols(mot1, mot2, ic1, ic2, min_overlap);
  NumericMatrix mot1_new = new_mots(0);
  NumericMatrix mot2_new = new_mots(1);
  NumericVector ic1_new = new_mots(2);
  NumericVector ic2_new = new_mots(3);

  NumericVector widths(2);
  widths[0] = mot1_new.ncol();
  widths[1] = mot2_new.ncol();
  int min_width = min(widths);
  int for_i = 1 + mot1_new.ncol() - min_width;
  int for_j = 1 + mot2_new.ncol() - min_width;
  NumericMatrix mot1_tmp(mot1_new.nrow(), min_width);
  NumericMatrix mot2_tmp(mot2_new.nrow(), min_width);
  NumericVector ic1_tmp(min_width);
  NumericVector ic2_tmp(min_width);
  NumericVector ans(for_i + for_j - 1);

  int method_i = 0;
  double ic1_mean;
  double ic2_mean;
  bool low_ic = false;

       if (method == "Euclidean") method_i = 1;
  else if (method == "Pearson")   method_i = 2;
  else if (method == "KL")        method_i = 3;

  for (int i = 0; i < for_i; ++i) {
    for (int j = 0; j < for_j; ++j) {
      
      for (int k = 0; k < min_width; ++k) {
        mot1_tmp(_, k) = mot1_new(_, k + i);
        mot2_tmp(_, k) = mot2_new(_, k + j);
        ic1_tmp[k] = ic1_new[k + i];
        ic2_tmp[k] = ic2_new[k + j];
      }

      ic1_mean = mean(na_omit(ic1_tmp));
      ic2_mean = mean(na_omit(ic2_tmp));
      if (ic1_mean < min_ic || ic2_mean < min_ic) low_ic = true;

      switch(method_i) {

        case 1: if (low_ic) ans(i + j) = 1;
                else ans(i + j) = motif_euclidean(mot1_tmp, mot2_tmp);
                break;
        case 2: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_pearson(mot1_tmp, mot2_tmp);
                break;
        case 3: if (low_ic) ans(i + j) = 10;
                else ans(i + j) = motif_kl(mot1_tmp, mot2_tmp);
                break;

      }

      low_ic = false;

    }
  }

  if (tryRC) {

    NumericVector final_ans(2);

    switch(method_i) {

      case 1: final_ans(0) = min(ans);
              final_ans(1) = min(ans_rc);
              return min(na_omit(final_ans));
      case 2: final_ans(0) = max(ans);
              final_ans(1) = max(ans_rc);
              return max(final_ans);
      case 3: final_ans(0) = min(ans);
              final_ans(1) = min(ans_rc);
              return min(final_ans);

    }
  
  }

  switch(method_i) {

    case 1: return min(ans);
    case 2: return max(ans);
    case 3: return min(ans);

  }

  return 999.0;

}

// [[Rcpp::export]]
NumericMatrix list_to_matrix_simil(List comparisons, StringVector mot_names,
    String method) {

  NumericMatrix out(mot_names.length(), mot_names.length());

  if (method == "Euclidean" || method == "KL") out.fill_diag(0.0);
  else if (method == "Pearson") out.fill_diag(1.0);

  int n = comparisons.length();
  for (int i = 0; i < n; ++i) {
    
    NumericVector toadd = comparisons(i);
    int n_ = toadd.length();
    for (int j = 0; j < n_; ++j) {
      out(j + i + 1, i) = toadd[j];
      out(i, j+ i + 1) = toadd[j];
    }
  
  }

  rownames(out) = mot_names;
  colnames(out) = mot_names;

  return out;

}
