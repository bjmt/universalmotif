#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
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

// [[Rcpp::export(rng = false)]]
NumericVector kmer_mat_to_probs_k1_cpp(IntegerMatrix bb_mat, NumericVector bkg,
    IntegerMatrix alph_sort) {

  // This function calculates the probability of finding a string of letters.
  // The only thing I'm unsure of is taking into account the length of the
  // sequence of origin; in which case the P-value calculation becomes:
  // ( 1 - P(A)^k1 * P(C)^k2 * P(G)^k3 * P(T)^k4 ) ^ N
  // Where N = total number of positions in a positions which could contain
  // the sequence.

  NumericVector probs_out(bb_mat.nrow(), 1.0);

  for (int i = 0; i < bb_mat.nrow(); ++i) {

    for (int j = 0; j < bb_mat.ncol(); ++j) {

      probs_out[i] *= bkg[alph_sort(bb_mat(i, j) - 1, j) - 1];

    }

  }

  return probs_out;

}

// // [[Rcpp::export(rng = false)]]
// NumericVector kmer_mat_to_probs_k2_cpp(IntegerMatrix bb_mat, NumericVector bkg1,
    // IntegerMatrix alph_sort, NumericVector bkg2) {
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
      // probs_out[i] *= bkg2[let3];
      // probs_out[i] *= bkg1[let1];
//
    // }
//
  // }
//
  // return probs_out;
//
// }

// // [[Rcpp::export(rng = false)]]
// NumericVector kmer_mat_to_probs_k3_cpp(IntegerMatrix bb_mat, NumericVector bkg,
    // IntegerMatrix alph_sort, NumericVector bkg2, NumericVector bkg3) {
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

// [[Rcpp::export(rng = false)]]
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

  path.attr("dim") = Dimension(path.length(), 1);

  return as<IntegerMatrix>(path);

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

// [[Rcpp::export(rng = false)]]
IntegerMatrix list_to_matrix(List paths) {

  IntegerMatrix tmp = paths[0];
  int n_rows = 0, n_cols = tmp.ncol();

  for (int i = 0; i < paths.length(); ++i) {
    IntegerMatrix tmp = paths[i];
    n_rows += tmp.nrow();
  }

  IntegerMatrix out(n_rows, n_cols);

  int pos = 0;
  for (int i = 0; i < paths.length(); ++i) {
    IntegerMatrix tmp = paths[i];
    for (int j = pos; j < pos + tmp.nrow(); ++j) {
      out(j, _) = tmp(j - pos, _);
    }
    pos += tmp.nrow();
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
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

    final_paths[i] = calc_next_subworker_cpp(paths_totry, scores_tmp, score);

  }

  return list_to_matrix(final_paths);

}

// [[Rcpp::export(rng = false)]]
IntegerVector expand_scores(IntegerMatrix scores) {

  int n_row = scores.nrow(), n_col = scores.ncol();
  IntegerMatrix expanded(pow(n_row, n_col), n_col);

  for (int i = 0; i < n_col; ++i) {
    expanded(_, i) = rep(rep_each(scores(_, i), pow(n_row, n_col - i - 1)),
                         pow(n_row, i + 1));
  }

  return rowSums(expanded);

}
