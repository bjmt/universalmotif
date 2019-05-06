#include <Rcpp.h>
#include <unordered_map>
using namespace Rcpp;

std::unordered_map<String, int> metrics_e = {
  {"EUCL", 1}, {"MEUCL", 2},
  {"PCC", 3}, {"MPCC", 4},
  {"KL", 5}, {"MKL", 8},
  {"MSW", 6}, {"SW", 7}
};

double motif_euclidean(const NumericMatrix &mot1, const NumericMatrix &mot2) {

  // Initially adapted from TFBSTools:::PWMEuclidean

  R_xlen_t mat_size = mot1.nrow() * mot1.ncol();
  NumericMatrix diff_mat(mot1.nrow(), mot1.ncol());

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    diff_mat[i] = mot1[i] - mot2[i];
    diff_mat[i] = pow(diff_mat[i], 2.0);
  }

  NumericVector diff_mat_colsums = colSums(diff_mat);
  for (R_xlen_t i = 0; i < diff_mat_colsums.size(); ++i) {
    diff_mat_colsums[i] = sqrt(diff_mat_colsums[i]);
  }

  LogicalVector to_keep = !is_na(diff_mat_colsums);
  diff_mat_colsums = diff_mat_colsums[to_keep];

  double k = sum(diff_mat_colsums);
  return k / sqrt(2.0);

}

double motif_euclidean_norm(const NumericMatrix &mot1, const NumericMatrix &mot2) {

  // Initially adapted from TFBSTools:::PWMEuclidean

  R_xlen_t mat_size = mot1.nrow() * mot1.ncol();
  NumericMatrix diff_mat(mot1.nrow(), mot1.ncol());

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    diff_mat[i] = mot1[i] - mot2[i];
    diff_mat[i] = pow(diff_mat[i], 2.0);
  }

  NumericVector diff_mat_colsums = colSums(diff_mat);
  for (R_xlen_t i = 0; i < diff_mat_colsums.size(); ++i) {
    diff_mat_colsums[i] = sqrt(diff_mat_colsums[i]);
  }

  LogicalVector to_keep = !is_na(diff_mat_colsums);
  diff_mat_colsums = diff_mat_colsums[to_keep];

  double k = sum(diff_mat_colsums);
  return k / diff_mat_colsums.length() / sqrt(2.0);

}

double motif_pearson(NumericMatrix mot1, NumericMatrix mot2) {

  // Initially adapted from TFBSTools:::PWMPearson

  R_xlen_t mat_size = mot1.nrow() * mot1.ncol();

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    mot1[i] -= 1.0 / mot1.nrow();
    mot2[i] -= 1.0 / mot1.nrow();
  }

  NumericMatrix mult_mat(mot1.nrow(), mot1.ncol());
  for (R_xlen_t i = 0; i < mat_size; ++i) {
    mult_mat[i] = mot1[i] * mot2[i];
  }

  NumericVector top = colSums(mult_mat);

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    mot1[i] = pow(mot1[i], 2.0);
    mot2[i] = pow(mot2[i], 2.0);
  }

  NumericVector mot1_colsums = colSums(mot1);
  NumericVector mot2_colsums = colSums(mot2);
  NumericVector bottom(mot1.ncol());
  for (R_xlen_t i = 0; i < mot1_colsums.size(); ++i) {
    bottom[i] = mot1_colsums[i] * mot2_colsums[i];
    bottom[i] = sqrt(bottom[i]);
  }

  NumericVector top_bot(mot1.ncol());
  for (R_xlen_t i = 0; i < mot1.ncol(); ++i) {
    top_bot[i] = top[i] / bottom[i];
  }

  LogicalVector to_keep = !is_na(top_bot);
  top_bot = top_bot[to_keep];

  return top_bot.size() * sum(top_bot);
}

double motif_pearson_norm(NumericMatrix mot1, NumericMatrix mot2) {

  // Initially adapted from TFBSTools:::PWMPearson

  R_xlen_t mat_size = mot1.nrow() * mot1.ncol();

  // should this be turned off for ICM?
  for (R_xlen_t i = 0; i < mat_size; ++i) {
    mot1[i] -= 1.0 / mot1.nrow();
    mot2[i] -= 1.0 / mot1.nrow();
  }

  NumericMatrix mult_mat(mot1.nrow(), mot1.ncol());
  for (R_xlen_t i = 0; i < mat_size; ++i) {
    mult_mat[i] = mot1[i] * mot2[i];
  }

  NumericVector top = colSums(mult_mat);

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    mot1[i] = pow(mot1[i], 2.0);
    mot2[i] = pow(mot2[i], 2.0);
  }

  NumericVector mot1_colsums = colSums(mot1);
  NumericVector mot2_colsums = colSums(mot2);
  NumericVector bottom(mot1.ncol());
  for (R_xlen_t i = 0; i < mot1_colsums.size(); ++i) {
    bottom[i] = mot1_colsums[i] * mot2_colsums[i];
    bottom[i] = sqrt(bottom[i]);
  }

  NumericVector top_bot(mot1.ncol());
  for (R_xlen_t i = 0; i < mot1.ncol(); ++i) {
    top_bot[i] = top[i] / bottom[i];
  }

  LogicalVector to_keep = !is_na(top_bot);
  top_bot = top_bot[to_keep];

  return 1.0 / top_bot.size() * sum(top_bot);
}

double motif_kl(const NumericMatrix &mot1, const NumericMatrix &mot2) {

  // Initially adapted from TFBSTools:::PWMKL

  R_xlen_t mat_size = mot1.nrow() * mot1.ncol();
  NumericMatrix new_mat(mot1.nrow(), mot1.ncol());

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    new_mat[i] = mot1[i] * log(mot1[i] / mot2[i]);
    new_mat[i] += mot2[i] * log(mot2[i] / mot1[i]);
  }

  NumericVector mat_colsums = colSums(new_mat);

  LogicalVector to_keep = !is_na(mat_colsums);
  mat_colsums = mat_colsums[to_keep];

  return 0.5 * sum(mat_colsums);

}

double motif_kl_norm(const NumericMatrix &mot1, const NumericMatrix &mot2) {

  // Initially adapted from TFBSTools:::PWMKL

  R_xlen_t mat_size = mot1.nrow() * mot1.ncol();
  NumericMatrix new_mat(mot1.nrow(), mot1.ncol());

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    new_mat[i] = mot1[i] * log(mot1[i] / mot2[i]);
    new_mat[i] += mot2[i] * log(mot2[i] / mot1[i]);
  }

  NumericVector mat_colsums = colSums(new_mat);

  LogicalVector to_keep = !is_na(mat_colsums);
  mat_colsums = mat_colsums[to_keep];

  return 0.5 / mat_colsums.size() * sum(mat_colsums);

}

double motif_sw(const NumericMatrix &mot1, const NumericMatrix &mot2) {

  // Adapted from RSAT compare-matrices

  R_xlen_t mat_size = mot1.nrow() * mot1.ncol();
  NumericMatrix res_mat(mot1.nrow(), mot1.ncol());

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    res_mat[i] = pow(mot1[i] - mot2[i], 2);
  }

  NumericVector mat_colsums = colSums(res_mat);
  mat_colsums = 2 - mat_colsums;

  LogicalVector to_keep = !is_na(mat_colsums);
  mat_colsums = mat_colsums[to_keep];

  return sum(mat_colsums);

}

double motif_sw_norm(const NumericMatrix &mot1, const NumericMatrix &mot2) {

  // Adapted from RSAT compare-matrices

  R_xlen_t mat_size = mot1.nrow() * mot1.ncol();
  NumericMatrix res_mat(mot1.nrow(), mot1.ncol());

  for (R_xlen_t i = 0; i < mat_size; ++i) {
    res_mat[i] = pow(mot1[i] - mot2[i], 2);
  }

  NumericVector mat_colsums = colSums(res_mat);
  mat_colsums = 2 - mat_colsums;

  LogicalVector to_keep = !is_na(mat_colsums);
  mat_colsums = mat_colsums[to_keep];

  return sum(mat_colsums) / mat_colsums.length();

}

// [[Rcpp::export(rng = false)]]
List add_cols(const NumericMatrix &mot1, const NumericMatrix &mot2,
    const NumericVector &ic1, const NumericVector &ic2, double overlap) {

  List out;

  R_xlen_t ncol_1 = mot1.ncol(), nrow_1 = mot1.nrow();
  R_xlen_t total_1 = ncol_1 * nrow_1;

  R_xlen_t ncol_2 = mot2.ncol(), nrow_2 = mot2.nrow();
  R_xlen_t total_2 = ncol_2 * nrow_2;

  R_xlen_t overlap1 = overlap, overlap2 = overlap;
  if (overlap < 1) {
    overlap1 = overlap * ncol_1;
    overlap2 = overlap * ncol_2;
  } else {
    overlap1 = overlap;
    overlap2 = overlap;
  }

  R_xlen_t ncol_1_toadd;
  if (overlap1 > ncol_2) ncol_1_toadd = 0;
  else ncol_1_toadd = ncol_2 - overlap1;

  R_xlen_t ncol_2_toadd;
  if (overlap2 > ncol_1) ncol_2_toadd = 0;
  else ncol_2_toadd = ncol_1 - overlap2;

  if (ncol_1_toadd == 0 || ncol_2_toadd == 0) {
    out = List::create(mot1, mot2, ic1, ic2);
    return out;
  }

  if (mot2.ncol() >= mot1.ncol()) {

    R_xlen_t total_1_toadd = ncol_1_toadd * nrow_1;
    R_xlen_t total_1_new = total_1_toadd * 2 + total_1;
    NumericMatrix mot1_new(nrow_1, ncol_1_toadd * 2 + ncol_1);
    NumericVector ic1_new(ic1.size() + ncol_1_toadd * 2, NA_REAL);

    for (R_xlen_t i = 0; i < total_1_new; ++i) {
      mot1_new[i] = NA_REAL;
    }
    for (R_xlen_t i = ncol_1_toadd; i < ncol_1_toadd + ncol_1; ++i) {
      mot1_new(_, i) = mot1(_, i - ncol_1_toadd);
    }

    for (R_xlen_t i = ncol_1_toadd; i < ncol_1_toadd + ic1.size(); ++i) {
      ic1_new[i] = ic1[i - ncol_1_toadd];
    }

  out = List::create(mot1_new, mot2, ic1_new, ic2);

  } else {

    R_xlen_t total_2_toadd = ncol_2_toadd * nrow_2;
    R_xlen_t total_2_new = total_2_toadd * 2 + total_2;
    NumericMatrix mot2_new(nrow_2, ncol_2_toadd * 2 + ncol_2);
    NumericVector ic2_new(ic2.size() + ncol_2_toadd * 2, NA_REAL);

    for (R_xlen_t i = 0; i < total_2_new; ++i) {
      mot2_new[i] = NA_REAL;
    }
    for (R_xlen_t i = ncol_2_toadd; i < ncol_2_toadd + ncol_2; ++i) {
      mot2_new(_, i) = mot2(_, i - ncol_2_toadd);
    }

    for (R_xlen_t i = ncol_2_toadd; i < ncol_2_toadd + ic2.size(); ++i) {
      ic2_new[i] = ic2[i - ncol_2_toadd];
    }

  out = List::create(mot1, mot2_new, ic1, ic2_new);

  }

  return out;

}

NumericMatrix rev_motif(const NumericMatrix &motif) {

  NumericMatrix out(motif.nrow(), motif.ncol());
  for (R_xlen_t i = 0; i < motif.nrow(); ++i) {
    out(i, _) = motif(motif.nrow() - i - 1, _);
  }
  NumericMatrix out2 = Rcpp::clone(out);
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    out(_, i) = out2(_, motif.ncol() - i - 1);
  }

  return out;

}

R_xlen_t get_align_len(const NumericMatrix &mot1, const NumericMatrix &mot2) {

  NumericVector mot1_ = mot1(1, _);
  NumericVector mot2_ = mot2(1, _);

  R_xlen_t out = 0;

  for (R_xlen_t i = 0; i < mot1.ncol(); ++i) {
    if (!NumericVector::is_na(mot1_[i]) && !NumericVector::is_na(mot2_[i]))
      out += 1;
  }

  return out;

}

double motif_simil_internal(NumericMatrix mot1, NumericMatrix mot2,
    String method, double min_overlap, bool tryRC, NumericVector ic1,
    NumericVector ic2, double min_ic, bool norm);

double motif_simil_rc(NumericMatrix mot1, NumericMatrix mot2,
    String method, double min_overlap, NumericVector ic1,
    NumericVector ic2, double min_ic, bool norm) {

  NumericMatrix mot2_rc = rev_motif(mot2);
  NumericVector ic2_rc = rev(ic2);
  double ans = motif_simil_internal(mot1, mot2_rc, method, min_overlap, false,
      ic1, ic2_rc, min_ic, norm);

  return ans;

}

// [[Rcpp::export(rng = false)]]
double motif_simil_internal(NumericMatrix mot1, NumericMatrix mot2,
    String method, double min_overlap, bool tryRC, NumericVector ic1,
    NumericVector ic2, double min_ic, bool norm) {

  // Initially adapted from TFBSTools::PWMSimilarity

  if (method == "KL" || method == "MKL") {
    mot1 = mot1 + 0.01;
    mot2 = mot2 + 0.01;
  }

  double len_total;
  if (mot1.ncol() >= mot2.ncol()) len_total = mot1.ncol();
  else len_total = mot2.ncol();

  NumericVector ans_rc;
  if (tryRC) ans_rc = motif_simil_rc(mot1, mot2, method, min_overlap,
      ic1, ic2, min_ic, norm);

  List new_mots = add_cols(mot1, mot2, ic1, ic2, min_overlap);
  NumericMatrix mot1_new = new_mots(0);
  NumericMatrix mot2_new = new_mots(1);
  NumericVector ic1_new = new_mots(2);
  NumericVector ic2_new = new_mots(3);

  NumericVector widths(2);
  widths[0] = mot1_new.ncol();
  widths[1] = mot2_new.ncol();
  R_xlen_t min_width = min(widths);
  R_xlen_t for_i = 1 + mot1_new.ncol() - min_width;
  R_xlen_t for_j = 1 + mot2_new.ncol() - min_width;
  NumericMatrix mot1_tmp(mot1_new.nrow(), min_width);
  NumericMatrix mot2_tmp(mot2_new.nrow(), min_width);
  NumericVector ic1_tmp(min_width);
  NumericVector ic2_tmp(min_width);
  NumericVector ans(for_i + for_j - 1);
  double align_len;

  int method_i = ::metrics_e[method];
  double ic1_mean, ic2_mean;
  bool low_ic = false;

  for (R_xlen_t i = 0; i < for_i; ++i) {
    for (R_xlen_t j = 0; j < for_j; ++j) {

      for (R_xlen_t k = 0; k < min_width; ++k) {
        mot1_tmp(_, k) = mot1_new(_, k + i);
        mot2_tmp(_, k) = mot2_new(_, k + j);
        ic1_tmp[k] = ic1_new[k + i];
        ic2_tmp[k] = ic2_new[k + j];
      }

      if (norm) align_len = get_align_len(mot1_tmp, mot2_tmp);
      else align_len = len_total;

      ic1_mean = mean(na_omit(ic1_tmp));
      ic2_mean = mean(na_omit(ic2_tmp));
      if (ic1_mean < min_ic || ic2_mean < min_ic) low_ic = true;

      switch(method_i) {

        case 1: if (low_ic) ans(i + j) = sqrt(2.0);
                else ans(i + j) = motif_euclidean(mot1_tmp, mot2_tmp) *
                                  (len_total / align_len);
                break;
        case 2: if (low_ic) ans(i + j) = sqrt(2.0) * align_len;
                else ans(i + j) = motif_euclidean_norm(mot1_tmp, mot2_tmp) *
                                  (len_total / align_len);
                break;
        case 3: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_pearson(mot1_tmp, mot2_tmp) *
                                  (align_len / len_total);
                break;
        case 4: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_pearson_norm(mot1_tmp, mot2_tmp) *
                                  (align_len / len_total);
                break;
        case 5: if (low_ic) ans(i + j) = 5 * align_len;
                else ans(i + j) = motif_kl(mot1_tmp, mot2_tmp) *
                                  (len_total / align_len);
                break;
        case 6: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_sw_norm(mot1_tmp, mot2_tmp) *
                                  (align_len / len_total);
                break;
        case 7: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_sw(mot1_tmp, mot2_tmp) *
                                  (align_len / len_total);
                break;
        case 8: if (low_ic) ans(i + j) = 5;
                else ans(i + j) = motif_kl_norm(mot1_tmp, mot2_tmp) *
                                  (len_total / align_len);
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
      case 2: final_ans(0) = min(ans);
              final_ans(1) = min(ans_rc);
              return min(na_omit(final_ans));
      case 3: final_ans(0) = max(ans);
              final_ans(1) = max(ans_rc);
              return max(final_ans);
      case 4: final_ans(0) = max(ans);
              final_ans(1) = max(ans_rc);
              return max(final_ans);
      case 5: final_ans(0) = min(ans);
              final_ans(1) = min(ans_rc);
              return min(final_ans);
      case 6: final_ans(0) = max(ans);
              final_ans(1) = max(ans_rc);
              return max(final_ans);
      case 7: final_ans(0) = max(ans);
              final_ans(1) = max(ans_rc);
              return max(final_ans);
      case 8: final_ans(0) = min(ans);
              final_ans(1) = min(ans_rc);
              return min(final_ans);

    }

  }

  switch(method_i) {

    case 1: return min(ans);
    case 2: return min(ans);
    case 3: return max(ans);
    case 4: return max(ans);
    case 5: return min(ans);
    case 6: return max(ans);
    case 7: return max(ans);
    case 8: return min(ans);

  }

  return NA_REAL;

}

// [[Rcpp::export(rng = false)]]
NumericMatrix list_to_matrix_simil(List comparisons, StringVector mot_names,
    String method) {

  NumericMatrix out(mot_names.size(), mot_names.size());

  switch (::metrics_e[method]) {
    case 4: out.fill_diag(1.0);
            break;
    case 6: out.fill_diag(2.0);
            break;
    default: out.fill_diag(0.0);
             break;
  }

  R_xlen_t n = comparisons.size();
  for (R_xlen_t i = 0; i < n; ++i) {

    NumericVector toadd = comparisons(i);
    R_xlen_t n_ = toadd.size();
    for (R_xlen_t j = 0; j < n_; ++j) {
      out(j + i + 1, i) = toadd[j];
      out(i, j+ i + 1) = toadd[j];
    }

  }

  rownames(out) = mot_names;
  colnames(out) = mot_names;

  return out;

}

NumericMatrix merge_mats(const NumericMatrix &mat1, const NumericMatrix &mat2,
    double weight1, double weight2) {

  NumericMatrix new_mat(mat1.nrow(), mat1.ncol());

  for (R_xlen_t i = 0; i < mat1.nrow(); ++i) {
    for (R_xlen_t j = 0; j < mat1.ncol(); ++j) {

      if (NumericVector::is_na(mat1(i, j)) && NumericVector::is_na(mat2(i, j)))
        new_mat(i, j) = NA_REAL;
      else if (NumericVector::is_na(mat1(i, j)))
        new_mat(i, j) = mat2(i, j);
      else if (NumericVector::is_na(mat2(i, j)))
        new_mat(i, j) = mat1(i, j);
      else
        new_mat(i, j) = (mat1(i, j) * weight1 + mat2(i, j) * weight2) / (weight1 + weight2);

    }
  }

  return new_mat;

}

// [[Rcpp::export(rng = false)]]
List merge_add_cols(List out) {

  NumericMatrix mot1 = out(0);
  NumericMatrix mot2 = out(1);
  R_xlen_t offset = out(2);

  R_xlen_t ncol1 = mot1.ncol();
  R_xlen_t nrow1 = mot1.nrow();
  R_xlen_t ncol2 = mot2.ncol();

  R_xlen_t inner_mat, ncol3;
  bool dontskip = true;
  NumericMatrix tmp_mat;

  if (ncol1 > ncol2) {
    ncol3 = ncol1;
    tmp_mat = mot2;
    inner_mat = ncol2;
  } else if (ncol2 > ncol1) {
    ncol3 = ncol2;
    tmp_mat = mot1;
    inner_mat = ncol1;
  } else dontskip = false;

  if (dontskip) {

    NumericMatrix new_mat(nrow1, ncol3); 

    for (R_xlen_t i = 0; i < offset; ++i) {
      for (R_xlen_t j = 0; j < nrow1; ++j) {
        new_mat(j, i) = NA_REAL;
      }
    }

    for (R_xlen_t i = offset; i < offset + inner_mat; ++i) {
      for (R_xlen_t j = 0; j < nrow1; ++j) {
        new_mat(j, i) = tmp_mat(j, i - offset);
      }
    }

    for (R_xlen_t i = offset + inner_mat; i < ncol3; ++i) {
      for (R_xlen_t j = 0; j < nrow1; ++j) {
        new_mat(j, i) = NA_REAL;
      }
    }

    if (ncol1 > ncol2) out(1) = new_mat;
    else if (ncol2 > ncol1) out(0) = new_mat;

  }

  return out;

}

// [[Rcpp::export(rng = false)]]
List merge_motifs_get_offset(NumericMatrix mot1, NumericMatrix mot2,
    String method, double min_overlap, NumericVector ic1,
    NumericVector ic2, double min_ic, bool norm) {

  if (method == "KL" || method == "MKL") {
    mot1 = mot1 + 0.01;
    mot2 = mot2 + 0.01;
  }

  double len_total;
  if (mot1.ncol() >= mot2.ncol()) len_total = mot1.ncol();
  else len_total = mot2.ncol();

  List new_mots = add_cols(mot1, mot2, ic1, ic2, min_overlap);
  NumericMatrix mot1_new = new_mots(0);
  NumericMatrix mot2_new = new_mots(1);
  NumericVector ic1_new = new_mots(2);
  NumericVector ic2_new = new_mots(3);

  NumericVector widths(2);
  widths[0] = mot1_new.ncol();
  widths[1] = mot2_new.ncol();
  R_xlen_t min_width = min(widths);
  R_xlen_t for_i = 1 + mot1_new.ncol() - min_width;
  R_xlen_t for_j = 1 + mot2_new.ncol() - min_width;
  NumericMatrix mot1_tmp(mot1_new.nrow(), min_width);
  NumericMatrix mot2_tmp(mot2_new.nrow(), min_width);
  NumericVector ic1_tmp(min_width);
  NumericVector ic2_tmp(min_width);
  NumericVector ans(for_i + for_j - 1);
  double align_len;

  int method_i = ::metrics_e[method];
  double ic1_mean;
  double ic2_mean;
  bool low_ic = false;

  for (R_xlen_t i = 0; i < for_i; ++i) {
    for (R_xlen_t j = 0; j < for_j; ++j) {

      for (R_xlen_t k = 0; k < min_width; ++k) {
        mot1_tmp(_, k) = mot1_new(_, k + i);
        mot2_tmp(_, k) = mot2_new(_, k + j);
        ic1_tmp[k] = ic1_new[k + i];
        ic2_tmp[k] = ic2_new[k + j];
      }

      if (norm) align_len = get_align_len(mot1_tmp, mot2_tmp);
      else align_len = len_total;

      ic1_mean = mean(na_omit(ic1_tmp));
      ic2_mean = mean(na_omit(ic2_tmp));
      if (ic1_mean < min_ic || ic2_mean < min_ic) low_ic = true;

      switch(method_i) {

        case 1: if (low_ic) ans(i + j) = sqrt(2.0);
                else ans(i + j) = motif_euclidean(mot1_tmp, mot2_tmp) *
                                  (len_total / align_len);
                break;
        case 2: if (low_ic) ans(i + j) = sqrt(2.0) * align_len;
                else ans(i + j) = motif_euclidean_norm(mot1_tmp, mot2_tmp) *
                                  (len_total / align_len);
                break;
        case 3: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_pearson(mot1_tmp, mot2_tmp) *
                                  (align_len / len_total);
                break;
        case 4: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_pearson_norm(mot1_tmp, mot2_tmp) *
                                  (align_len / len_total);
                break;
        case 5: if (low_ic) ans(i + j) = 5 * align_len;
                else ans(i + j) = motif_kl(mot1_tmp, mot2_tmp) *
                                  (len_total / align_len);
                break;
        case 6: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_sw_norm(mot1_tmp, mot2_tmp) *
                                  (align_len / len_total);
                break;
        case 7: if (low_ic) ans(i + j) = 0;
                else ans(i + j) = motif_sw(mot1_tmp, mot2_tmp) *
                                  (align_len / len_total);
                break;
        case 8: if (low_ic) ans(i + j) = 5;
                else ans(i + j) = motif_kl_norm(mot1_tmp, mot2_tmp) *
                                  (len_total / align_len);
                break;

      }

      low_ic = false;

    }
  }

  R_xlen_t offset;
  double score;
  switch(method_i) {

    case 1: offset = which_min(ans);
            score = min(ans);
            break;
    case 2: offset = which_min(ans);
            score = min(ans);
            break;
    case 3: offset = which_max(ans);
            score = max(ans);
            break;
    case 4: offset = which_max(ans);
            score = max(ans);
            break;
    case 5: offset = which_min(ans);
            score = min(ans);
            break;
    case 6: offset = which_max(ans);
            score = max(ans);
            break;
    case 7: offset = which_max(ans);
            score = max(ans);
            break;
    case 8: offset = which_min(ans);
            score = min(ans);
            break;

  }

  List out = List::create(_["mot1_new"] = mot1_new, _["mot2_new"] = mot2_new,
                          _["offset"] = offset, _["score"] = score);

  out = merge_add_cols(out);

  return out;

}

// [[Rcpp::export(rng = false)]]
NumericMatrix merge_motifs_internal(NumericMatrix mot1, NumericMatrix mot2,
    String method, double min_overlap, bool tryRC, NumericVector ic1,
    NumericVector ic2, double min_ic, double weight1, double weight2,
    bool norm) {

  List out = merge_motifs_get_offset(mot1, mot2, method, min_overlap,
      ic1, ic2, min_ic, norm);

  NumericMatrix mot1_;
  NumericMatrix mot2_;

  if (tryRC) {

    double score;
    double score_rc;
    NumericMatrix mot2_rc = rev_motif(mot2);
    NumericVector ic2_rc = rev(ic2);
    List out_rc = merge_motifs_get_offset(mot1, mot2_rc, method, min_overlap,
        ic1, ic2_rc, min_ic, norm);
    score = out(3);
    score_rc = out_rc(3);

    if (score >= score_rc) {

      out = merge_add_cols(out);
      NumericMatrix mot1_ = out(0);
      NumericMatrix mot2_ = out(1);
      NumericMatrix final_mot = merge_mats(mot1_, mot2_, weight1, weight2);
      return final_mot;

    } else {

      out_rc = merge_add_cols(out_rc);
      NumericMatrix mot1_ = out_rc(0);
      NumericMatrix mot2_ = out_rc(1);
      NumericMatrix final_mot = merge_mats(mot1_, mot2_, weight1, weight2);
      return final_mot;

    }

  } else {

    out = merge_add_cols(out);
    NumericMatrix mot1_ = out(0);
    NumericMatrix mot2_ = out(1);
    NumericMatrix final_mot = merge_mats(mot1_, mot2_, weight1, weight2);
    return final_mot;

  }

  NumericMatrix final_mot;
  return final_mot;

}
