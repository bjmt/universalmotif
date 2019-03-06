#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector calc_scores_cpp(IntegerMatrix paths, IntegerMatrix score_mat) {

  IntegerVector final_scores(paths.nrow());
  int tmp_score;

  for (int i = 0; i < paths.nrow(); ++i) {
    tmp_score = 0;
    for (int j = 0; j < paths.ncol(); ++j) {
      tmp_score += score_mat(paths(i, j) - 1, j);
    } 
    final_scores[i] = tmp_score;
  }

  return final_scores;

}

// [[Rcpp::export]]
NumericVector kmer_mat_to_probs_k1_cpp(IntegerMatrix bb_mat, NumericVector bkg,
    IntegerMatrix alph_sort) {

  // NumericMatrix prob_mat(bb_mat.nrow(), bb_mat.ncol());
  NumericVector probs_out(bb_mat.nrow(), 1.0);

  for (int i = 0; i < bb_mat.nrow(); ++i) {

    for (int j = 0; j < bb_mat.ncol(); ++j) {

      probs_out[i] *= bkg[alph_sort(bb_mat(i, j) - 1, j) - 1];

    }

  }

  return probs_out;

}

// // [[Rcpp::export]]
// NumericVector kmer_mat_to_probs_k2_cpp(IntegerMatrix bb_mat, NumericVector bkg,
    // IntegerMatrix alph_sort) {
//
  // NumericVector probs_out(bb_mat.nrow(), 1.0);
  // int let1, let2, let3;
//
  // for (int i = 0; i < bb_mat.nrow(); ++i) {
//
    // for (int j = 0; j < bb_mat.ncol() - 1; ++j) {
//
      // let1 = alph_sort(bb_mat(i, j) - 1, j) - 1;
      // let2 = alph_sort(bb_mat(i, j + 1) - 1, j + 1) - 1;
//
      // let3 = let1 * alph_sort.nrow() + let2;
//
      // probs_out[i] *= bkg[let3];
//
    // }
//
  // }
//
  // return probs_out;
//
// }

// // [[Rcpp::export]]
// NumericVector kmer_mat_to_probs_k3_cpp(IntegerMatrix bb_mat, NumericVector bkg,
    // IntegerMatrix alph_sort) {
//
  // NumericVector probs_out(bb_mat.nrow(), 1.0);
  // int let1, let2, let3, let4;
  // int alph_len = alph_sort.nrow();
//
  // for (int i = 0; i < bb_mat.nrow(); ++i) {
//
    // for (int j = 0; j < bb_mat.ncol() - 2; ++j) {
//
      // let1 = alph_sort(bb_mat(i, j) - 1, j) - 1;
      // let2 = alph_sort(bb_mat(i, j + 1) - 1, j + 1) - 1;
      // let3 = alph_sort(bb_mat(i, j + 2) - 1, j + 2) - 1;
//
      // let4 = let1 * alph_len * alph_len + let2 * alph_len + let3;
//
      // probs_out[i] *= bkg[let4];
//
    // }
//
  // }
//
  // return probs_out;
//
// }

// [[Rcpp::export]]
IntegerMatrix init_paths_cpp(IntegerMatrix score_mat, int score, 
    int max_score) {

  int alph_len = score_mat.nrow();

  IntegerVector path(alph_len);
  for (int i = 0; i < alph_len; ++i) {
    path[i] = i + 1;
  }

  LogicalVector tokeep(alph_len, false);

  int tmp_score;

  for (int i = 0; i < alph_len; ++i) {
    tmp_score = score_mat(i, 0) + max_score;
    if (tmp_score >= score) tokeep[i] = true;
  }

  path = path[tokeep];

  IntegerMatrix out(path.length(), 1);
  for (int i = 0; i < path.length(); ++i) {
    out(i, 0) = path[i];
  }

  return out;

}

IntegerMatrix calc_next_subworker_cpp(IntegerMatrix paths_totry,
    IntegerVector scores_tmp, int score) {

  int numrows = 0;

  for (int i = 0; i < scores_tmp.length(); ++i) {
    if (scores_tmp[i] >= score) numrows += 1;
  }

  IntegerMatrix new_mat(numrows, paths_totry.ncol());

  int currentrow = 0;
  for (int i = 0; i < paths_totry.nrow(); ++i) {
    if (scores_tmp[i] >= score) {
      new_mat(currentrow, _) = paths_totry(i, _);
      currentrow += 1;
    }
  }

  return new_mat;

}

// [[Rcpp::export]]
IntegerMatrix list_to_matrix(List paths) {

  IntegerMatrix tmp = paths[0];

  int num_mat(paths.length());
  int mat_ncols(tmp.ncol());

  IntegerVector path_nrows(num_mat);

  for (int i = 0; i < num_mat; ++i) {
    IntegerMatrix tmp = paths[i];
    path_nrows[i] = tmp.nrow();
  }

  int total_rows = sum(path_nrows);

  IntegerVector out_pre(total_rows * mat_ncols);

  int pre_index = 0;

  for (int i = 0; i < num_mat; ++i) {
    IntegerMatrix tmp = paths[i];
    for (int j = 0; j < tmp.nrow(); ++j) {
      for (int b = 0; b < tmp.ncol(); ++b) {
        out_pre[pre_index] = tmp(j, b);
        pre_index += 1;
      }
    }
  }

  pre_index = 0;
  IntegerMatrix out(total_rows, mat_ncols);
  for (int i = 0; i < out.nrow(); ++i) {
    for (int j = 0; j < out.ncol(); ++j) {
      out(i, j) = out_pre[pre_index];
      pre_index += 1;
    }
  }

  return out;

}

// [[Rcpp::export]]
IntegerMatrix calc_next_path_cpp(IntegerMatrix score_mat, IntegerMatrix paths,
    int score, int max_score) {

  int alph_len = score_mat.nrow();
  IntegerMatrix next_paths(alph_len, 1);
  for (int i = 0; i < alph_len; ++i) {
    next_paths(i, 0) = i + 1;
  }

  List final_paths(paths.nrow());

  IntegerMatrix paths_totry(score_mat.nrow(), paths.ncol() + 1);

  IntegerVector scores_tmp;

  for (int i = 0; i < paths.nrow(); ++i) {

    for (int j = 0; j < paths_totry.nrow(); ++j) {
      for (int b = 0; b < paths_totry.ncol() - 1; ++b) {
        paths_totry(j, b) = paths(i, b);
      }
    }

    paths_totry(_, paths_totry.ncol() - 1) = next_paths;

    scores_tmp = calc_scores_cpp(paths_totry, score_mat) + max_score;

    final_paths(i) = calc_next_subworker_cpp(paths_totry, scores_tmp, score);
  
  }

  return list_to_matrix(final_paths);
  // return final_paths;

}

// // [[Rcpp::export]]
IntegerMatrix branch_and_bound_cpp(IntegerMatrix score_mat, int min_score) {

  IntegerVector max_scores(score_mat.ncol() + 1);
  for (int i = 0; i < max_scores.length() - 1; ++i) {
    max_scores[i] = max(score_mat(_, i));
  }

  for (int i = 0; i < max_scores.length(); ++i) {
    for (int j = i + 1; j < max_scores.length(); ++j) {
      max_scores[i] += max_scores[j];
    }
  }

  if (min_score > max_scores(0)) min_score = max_scores(0) - 0.0001;

  int mot_len = score_mat.ncol();

  IntegerMatrix paths = init_paths_cpp(score_mat, min_score, max_scores[1]);
  if (mot_len == 1) return paths;

  for (int i = 1; i < mot_len; ++i) {
    paths = calc_next_path_cpp(score_mat, paths, min_score, max_scores[i + 1]);
  }

  return paths;

}

// [[Rcpp::export]]
NumericVector calc_final_probs_cpp(List all_probs, List all_scores,
    int score) {

  IntegerVector scores0 = all_scores[0];
  IntegerVector scores1 = all_scores[1];
  NumericVector probs0 = all_probs[0];
  NumericVector probs1 = all_probs[1];
  NumericVector final_probs(scores0.length());

  double tmp_prob;
  for (int i = 0; i < final_probs.length(); ++i) {
    tmp_prob = probs0[i];
    NumericVector tmp_probs_2 = probs1;
    LogicalVector which_probs_2 = scores1 > score - scores0(i);
    tmp_probs_2 = tmp_probs_2[which_probs_2];
    final_probs[i] = tmp_prob * sum(tmp_probs_2);
  }

  return final_probs;

}
