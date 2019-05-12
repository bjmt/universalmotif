#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <vector>
#include <deque>
#include <numeric>
#include "types.h"

bool alph_sort(int i, int j) {
  return i > j;
}

bool position_sort(vec_int_t i, vec_int_t j) {
  return *std::max_element(i.begin(), i.end()) > *std::max_element(j.begin(), j.end());
}

int get_split_max_sum(const list_int_t mot) {

  int answer = 0;
  for (std::size_t i = 0; i < mot.size(); ++i) {
    answer += *std::max_element(mot[i].begin(), mot[i].end());
  }

  return answer;

}

vec_int_t get_split_mins(const vec_int_t &split_maxsums, const int &score) {

  vec_int_t split_mins(split_maxsums.size());
  int msum = std::accumulate(split_maxsums.begin(), split_maxsums.end(), 0);

  for (std::size_t i = 0; i < split_maxsums.size(); ++i) {
    split_mins[i] = score - (msum - split_maxsums[i]);
  }

  return split_mins;

}

list_int_t bb_init_paths(const list_int_t &mot, const int &score,
    const int &max_score, const std::size_t &alphlen) {

  vec_int_t ipaths(alphlen);
  std::iota(ipaths.begin(), ipaths.end(), 0);

  list_int_t paths(1);
  paths[0].reserve(alphlen);

  for (std::size_t i = 0; i < alphlen; ++i) {
    if (mot[0][i] + max_score >= score) paths[0].push_back(ipaths[i]);
  }

  return paths;

}

vec_int_t bb_calc_scores(const list_int_t &paths, const list_int_t &mot) {

  std::size_t pnrow = paths[0].size();
  std::size_t pncol = paths.size();
  vec_int_t final_scores(pnrow);
  int tmp_score;

  for (std::size_t i = 0; i < pnrow; ++i) {
    tmp_score = 0;
    for (std::size_t j = 0; j < pncol; ++j) {
      tmp_score += mot[j][paths[j][i]];
    }
    final_scores[i] = tmp_score;
  }

  return final_scores;

}

int calc_row_score(const vec_int_t &tpath, const list_int_t &mot) {

  int score = 0;
  for (std::size_t i = 0; i < tpath.size(); ++i) {
    score += mot[i][tpath[i]];
  }

  return score;

}


void fill_tpath(vec_int_t &tpath, const list_int_t &paths, const std::size_t &j) {

  for (std::size_t i = 0; i < paths.size(); ++i) {
    tpath[i] = paths[i][j];
  }

  return;

}

list_int_t bb_path_get_next(const list_int_t &mot, const list_int_t &paths,
    const int &score, const int &max_score, const std::size_t &alphlen) {

  std::size_t npaths = paths[0].size();

  list_int_t newpaths(paths.size() + 1);
  for (std::size_t i = 0; i < newpaths.size(); ++i) {
    /* might be a better solution here, ratio of reserved/filled can be huge */
    newpaths[i].reserve(alphlen * npaths);
  }

  vec_int_t tpath(paths.size() + 1);

  for (std::size_t i = 0; i < alphlen; ++i) {

    for (std::size_t j = 0; j < npaths; ++j) {

      fill_tpath(tpath, paths, j);
      tpath[tpath.size() - 1] = int(i);

      if (calc_row_score(tpath, mot) + max_score >= score) {
        for (std::size_t b = 0; b < tpath.size(); ++b) {
          newpaths[b].push_back(tpath[b]);
        }
      }

    }

  }

  return newpaths;

}

list_int_t branch_and_bound_kmers_cpp(const list_int_t &mot, const int &minscore) {

  std::size_t alphlen = mot[0].size();

  vec_int_t pos_maxes(mot.size() + 1, 0);
  for (std::size_t i = 0; i < mot.size(); ++i) {
    pos_maxes[i] = *std::max_element(mot[i].begin(), mot[i].end());
  }

  for (std::size_t i = 0; i < mot.size() - 1; ++i) {
    for (std::size_t j = i + 1; j < mot.size(); ++j) {
      pos_maxes[i] += pos_maxes[j];
    }
  }

  list_int_t paths = bb_init_paths(mot, minscore, pos_maxes[1], alphlen);
  if (mot.size() == 1) return paths;

  for (std::size_t i = 1; i < mot.size(); ++i) {
    paths = bb_path_get_next(mot, paths, minscore, pos_maxes[i + 1], alphlen);
  }

  return paths;

}

vec_lnum_t bb_paths_to_probs(const list_int_t &paths, const list_int_t &alph_indices,
    const vec_num_t &bkg) {

  std::size_t pnrow = paths[0].size();
  std::size_t pncol = paths.size();
  vec_lnum_t probs(pnrow, 1.0);

  for (std::size_t i = 0; i < pnrow; ++i) {
    for (std::size_t j = 0; j < pncol; ++j) {
      probs[i] *= bkg[alph_indices[j][paths[j][i]]];
    }
  }

  return probs;

}

bool sort_motpos(std::size_t j, std::size_t b, const vec_int_t mot) {

  return mot[j] > mot[b];

}

long double motif_pvalue_single(list_int_t mot, const double score,
    const int k, const vec_num_t bkg) {

  int iscore = score * 1000;
  std::size_t alphlen = bkg.size();
  std::size_t motlen = mot.size();

  if (k < 1)
    Rcpp::stop("k must be greater than 0");

  if (mot.size() == 0 || mot[0].size() == 0)
    Rcpp::stop("empty motif");
  if (bkg.size() == 0)
    Rcpp::stop("empty bkg vector");

  if (mot[0].size() != bkg.size())
    Rcpp::stop("bkg vector length does not match motif row number");

  int mmin = 0, mmax = 0;
  for (std::size_t i = 0; i < motlen; ++i) {
    mmin += *std::min_element(mot[i].begin(), mot[i].end());
    mmax += *std::max_element(mot[i].begin(), mot[i].end());
  }
  if (iscore > mmax || iscore < mmin) {
    Rcpp::stop("input score is outside of min/max motif score range");
  }

  std::sort(mot.begin(), mot.end(), position_sort);

  list_int_t sorted_alph_indices(motlen);
  for (std::size_t i = 0; i < motlen; ++i) {
    vec_int_t indices(alphlen);
    std::iota(indices.begin(), indices.end(), 0);
    sort(std::begin(indices), std::end(indices),
         std::bind(sort_motpos,std::placeholders::_1, std::placeholders::_2, mot[i]));
    sorted_alph_indices[i] = indices;
    std::sort(mot[i].begin(), mot[i].end(), alph_sort);
  }

  list_mat_t mot_split;
  list_mat_t alph_indices_split;
  if (int(motlen) <= k) {

    mot_split.push_back(mot);
    alph_indices_split.push_back(sorted_alph_indices);

  } else {

    int nsplit = motlen / k;
    if (motlen % k > 0) ++nsplit;
    int counter = 0;

    mot_split = list_mat_t(nsplit);
    alph_indices_split = list_mat_t(nsplit);
    for (int i = 0; i < nsplit; ++i) {
      mot_split[i].reserve(k);
      alph_indices_split[i].reserve(k);
      for (int j = 0; j < k; ++j) {
        if (counter == int(motlen)) break;
        mot_split[i].push_back(mot[counter]);
        alph_indices_split[i].push_back(sorted_alph_indices[counter]);
        ++counter;
      }
    }

  }

  vec_int_t split_maxsums(mot_split.size());
  for (std::size_t i = 0; i < mot_split.size(); ++i) {
    split_maxsums[i] = get_split_max_sum(mot_split[i]);
  }

  vec_int_t split_minsums = get_split_mins(split_maxsums, score);

  list_mat_t all_paths(mot_split.size());
  for (std::size_t i = 0; i < all_paths.size(); ++i) {
    all_paths[i] = branch_and_bound_kmers_cpp(mot_split[i], split_minsums[i]);
  }

  list_lnum_t all_probs(mot_split.size());
  for (std::size_t i = 0; i < mot_split.size(); ++i) {
    all_probs[i] = bb_paths_to_probs(all_paths[i], alph_indices_split[i], bkg);
  }

  long double pvalue = 0;

  if (mot_split.size() == 1) {
    pvalue = std::accumulate(all_probs[0].begin(), all_probs[0].end(), 0.0);
    return pvalue;
  } /* --> if (k >= motlen), then function end here */

  list_int_t all_scores(mot_split.size());
  for (std::size_t i = 0; i < mot_split.size(); ++i) {
    all_scores[i] = bb_calc_scores(all_paths[i], mot_split[i]);
  }

  vec_int_t split_maxes(all_scores.size());
  for (std::size_t i = 0; i < all_scores.size(); ++i) {
    split_maxes[i] = *std::max_element(all_scores[i].begin(), all_scores[i].end());
  }

  if (mot_split.size() > 2) {

    for (std::size_t i = 0; i < mot_split.size() - 2; ++i) {

      for (std::size_t j = 0; j < all_probs[i + 1].size(); ++j) {

        int tscore = 0;
        for (std::size_t b = 0; b <= i; ++b) {
          tscore += split_maxes[b];
        }

        vec_bool_t gscores(all_scores[i + 1].size(), false);
        for (std::size_t b = 0; b < all_scores[i + 2].size(); ++b) {
          if (all_scores[i + 2][b] > iscore - all_scores[i + 1][i] - tscore)
            gscores[b] = true;
        }

        vec_lnum_t tprobs;
        tprobs.reserve(all_probs[i + 2].size());
        for (std::size_t b = 0; b < gscores.size(); ++b) {
          if (gscores[b]) tprobs.push_back(all_probs[i + 2][b]);
        }

        all_probs[i + 1][j] *= std::accumulate(tprobs.begin(), tprobs.end(), 0.0);

      }

    }

  }

  vec_lnum_t final_probs(all_scores[0].size());

  for (std::size_t i = 0; i < final_probs.size(); ++i) {

    vec_bool_t gscores(all_scores[1].size(), false);
    for (std::size_t j = 0; j < gscores.size(); ++j) {
      if (all_scores[1][j] > iscore - all_scores[0][i])
        gscores[j] = true;
    }

    vec_lnum_t tprobs;
    tprobs.reserve(all_probs[1].size());
    for (std::size_t j = 0; j < gscores.size(); ++j) {
      if (gscores[j]) tprobs.push_back(all_probs[1][j]);
    }

    final_probs[i] = all_probs[0][i] * std::accumulate(tprobs.begin(), tprobs.end(), 0.0);

  }

  pvalue = std::accumulate(final_probs.begin(), final_probs.end(), 0.0);

  return pvalue;

}

list_int_t R_to_cpp_motif(const Rcpp::NumericMatrix &motif) {

  list_int_t mat(motif.ncol());
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    mat[i].reserve(motif.nrow());
    for (R_xlen_t j = 0; j < motif.nrow(); ++j) {
      mat[i].push_back(int(motif(j, i) * 1000));
    }
  }

  return mat;

}

list_int_t R_to_cpp_motif(const Rcpp::IntegerMatrix &motif) {

  list_int_t mat(motif.ncol());
  for (R_xlen_t i = 0; i < motif.ncol(); ++i) {
    mat[i].reserve(motif.nrow());
    for (R_xlen_t j = 0; j < motif.nrow(); ++j) {
      mat[i].push_back(motif(j, i));
    }
  }

  return mat;

}
/* C++ ENTRY ---------------------------------------------------------------- */

// [[Rcpp::export(rng = false)]]
std::vector<long double> motif_pvalue_cpp(const Rcpp::List &motifs,
    const Rcpp::List &bkg, const std::vector<double> &scores, const int &k = 6,
    const int &nthreads = 1) {

  list_mat_t vmotifs(motifs.size());
  for (R_xlen_t i = 0; i < motifs.size(); ++i) {
    Rcpp::NumericMatrix single = motifs[i];
    vmotifs[i] = R_to_cpp_motif(single);
  }
  std::vector<vec_num_t> vbkg(bkg.size());
  for (R_xlen_t i = 0; i < bkg.size(); ++i) {
    Rcpp::NumericVector single = bkg[i];
    vbkg[i] = Rcpp::as<std::vector<double>>(single);
  }

  std::vector<long double> pvalues(scores.size());
  RcppThread::parallelFor(0, pvalues.size(),
      [&vmotifs, &pvalues, &scores, &k, &vbkg] (std::size_t i) {
        pvalues[i] = motif_pvalue_single(vmotifs[i], scores[i], k, vbkg[i]);
      }, nthreads);

  return pvalues;

}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix branch_and_bound_cpp_exposed(Rcpp::IntegerMatrix mat,
    const int score) {

  list_int_t cppmat = R_to_cpp_motif(mat);

  list_int_t paths = branch_and_bound_kmers_cpp(cppmat, score);

  Rcpp::IntegerMatrix out(paths[0].size(), paths.size());
  for (std::size_t i = 0; i < paths.size(); ++i) {
    Rcpp::IntegerVector tmp = Rcpp::wrap(paths[i]);
    out(Rcpp::_, i) = tmp;
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector expand_scores(const Rcpp::IntegerMatrix &scores) {

  R_xlen_t n_row = scores.nrow(), n_col = scores.ncol();
  Rcpp::IntegerMatrix expanded(pow(n_row, n_col), n_col);

  for (R_xlen_t i = 0; i < n_col; ++i) {
    expanded(Rcpp::_, i) = Rcpp::rep(Rcpp::rep_each(scores(Rcpp::_, i),
          pow(n_row, n_col - i - 1)), pow(n_row, i + 1));
  }

  return Rcpp::rowSums(expanded);

}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerMatrix paths_alph_unsort(const Rcpp::IntegerMatrix &paths,
    const Rcpp::IntegerMatrix &alph) {

  Rcpp::IntegerMatrix paths2 = Rcpp::clone(paths);

  for (R_xlen_t i = 0; i < paths2.ncol(); ++i) {
    for (R_xlen_t j = 0; j < paths2.nrow(); ++j) {
      paths2(j, i) = alph(paths(j, i), i);
    }
  }

  return paths2;

}

// [[Rcpp::export(rng = false)]]
Rcpp::StringVector paths_to_alph(const Rcpp::IntegerMatrix &paths,
    const Rcpp::StringVector &alph) {

  Rcpp::StringMatrix out(paths.nrow(), paths.ncol());

  for (R_xlen_t i = 0; i < paths.nrow(); ++i) {
    for (R_xlen_t j = 0; j < paths.ncol(); ++j) {
      out(i, j) = alph[paths(i, j)];
    }
  }

  Rcpp::StringVector outstr(out.nrow());
  for (R_xlen_t i = 0; i < out.nrow(); ++i) {
    outstr[i] = Rcpp::collapse(out(i, Rcpp::_));
  }

  return outstr;

}
