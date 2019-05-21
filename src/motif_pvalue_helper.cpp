#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include "types.h"
#include "utils-internal.h"

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

  long double pvalue = 0;

  if (int(motlen) <= k) {

    list_int_t paths = branch_and_bound_kmers_cpp(mot, iscore);
    vec_lnum_t probs = bb_paths_to_probs(paths, sorted_alph_indices, bkg);
    pvalue = std::accumulate(probs.begin(), probs.end(), 0.0);
    return pvalue;

  }

  int nsplit = int(motlen) / k;
  if (int(motlen) % k > 0 || nsplit == 0) ++nsplit;

  list_mat_t mot_split;
  list_mat_t alph_indices_split;

  mot_split.reserve(nsplit);
  alph_indices_split.reserve(nsplit);

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

  vec_int_t split_maxsums(mot_split.size());
  for (std::size_t i = 0; i < mot_split.size(); ++i) {
    split_maxsums[i] = get_split_max_sum(mot_split[i]);
  }

  vec_int_t split_minsums = get_split_mins(split_maxsums, iscore);

  list_mat_t all_paths(mot_split.size());
  for (std::size_t i = 0; i < all_paths.size(); ++i) {
    all_paths[i] = branch_and_bound_kmers_cpp(mot_split[i], split_minsums[i]);
  }

  list_lnum_t all_probs(mot_split.size());
  for (std::size_t i = 0; i < mot_split.size(); ++i) {
    all_probs[i] = bb_paths_to_probs(all_paths[i], alph_indices_split[i], bkg);
  }

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

list_int_t expand_scores_cpp(const list_int_t &mot) {

  std::size_t nrow = mot[0].size();
  int npaths = pow(nrow, mot.size());

  list_int_t expanded(mot.size(), vec_int_t(npaths));

  std::size_t step;
  int counter;
  for (std::size_t i = 0; i < expanded.size(); ++i) {

    step = pow(nrow, i);
    counter = 0;

    while (counter < npaths) {
      for (std::size_t j = 0; j < nrow; ++j) {
        for (std::size_t b = 0; b < step; ++b) {
          expanded[i][counter] = mot[i][j];
          ++counter;
        }
      }
    }

  }

  return expanded;

}

vec_int_t rowsums_cpp(const list_int_t &mot) {

  std::size_t nrow = mot[0].size();
  vec_int_t rowsums(nrow, 0);

  for (std::size_t i = 0; i < nrow; ++i) {
    for (std::size_t j = 0; j < mot.size(); ++j) {
      rowsums[i] += mot[j][i];
    }
  }

  return rowsums;

}

double motif_score_single(const list_int_t &mot, const int k, const int randtries,
    std::default_random_engine gen, const double pval) {

  std::size_t motlen = mot.size();

  if (k < 1)
    Rcpp::stop("k must be greater than 0");

  if (mot.size() == 0 || mot[0].size() == 0)
    Rcpp::stop("empty motif");

  if (randtries < 1)
    Rcpp::stop("rand.tries must be greater than zero");

  if (int(motlen) > k) {

    int nsplit = int(motlen) / k;
    if (int(motlen) % k > 0 || nsplit == 0) ++nsplit;

    list_mat_t mot_split = list_mat_t(nsplit);

    int counter = 0;

    for (int i = 0; i < nsplit; ++i) {
      mot_split[i].reserve(k);
      for (int j = 0; j < k; ++j) {
        if (counter == int(motlen)) break;
        mot_split[i].push_back(mot[counter]);
        ++counter;
      }
    }

    list_mat_t paths(nsplit);
    list_int_t scores(nsplit);

    for (int i = 0; i < nsplit; ++i) {
      paths[i] = expand_scores_cpp(mot_split[i]);
    }

    for (int i = 0; i < nsplit; ++i) {
      scores[i] = rowsums_cpp(paths[i]);
    }

    vec_int_t scores0 = scores[0];
    vec_int_t tquants(randtries);

    vec_int_t scorelens(scores.size());
    for (std::size_t i = 0; i < scores.size(); ++i) {
      scorelens[i] = scores[i].size();
    }

    for (int i = 0; i < randtries; ++i) {

      for (int j = 1; j < nsplit; ++j) {
        for (std::size_t b = 0; b < scores[0].size(); ++b) {
          scores[0][b] += scores[j][gen() % scorelens[j]];
        }
      }

      int which = scores[0].size() * (1.0 - pval);
      if (which == int(scores[0].size())) --which;

      std::nth_element(scores[0].begin(), scores[0].begin() + which, scores[0].end());

      tquants[i] = scores[0][which];

      scores[0].assign(scores0.begin(), scores0.end());

    }

    double final_score = std::accumulate(tquants.begin(), tquants.end(), 0.0);

    final_score /= double(tquants.size());
    final_score /= 1000.0;

    return final_score;

  }

  list_int_t paths = expand_scores_cpp(mot);
  vec_int_t scores = rowsums_cpp(paths);

  int which = scores.size() * (1.0 - pval);
  if (which == int(scores.size())) --which;

  std::nth_element(scores.begin(), scores.begin() + which, scores.end());

  double final_score = scores[which];
  final_score /= 1000.0;

  return final_score;

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
std::vector<double> motif_score_cpp(const Rcpp::List &motifs,
    const std::vector<double> &pvals, const int seed = 1, const int k = 6,
    const int nthreads = 1, const int randtries = 100) {

  list_mat_t vmots(motifs.size());
  for (R_xlen_t i = 0; i < motifs.size(); ++i) {
    Rcpp::NumericMatrix single = motifs[i];
    vmots[i] = R_to_cpp_motif(single);
  }

  vec_num_t scores(pvals.size());
  RcppThread::parallelFor(0, scores.size(),
      [&vmots, &scores, &seed, &k, &randtries, &pvals] (std::size_t i) {
      std::default_random_engine gen(seed * (int(i) + 1));
        scores[i] = motif_score_single(vmots[i], k, randtries, gen, pvals[i]);
      }, nthreads);

  return scores;

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

  list_int_t scores2 = R_to_cpp_motif(scores);
  list_int_t scores3 = expand_scores_cpp(scores2);
  vec_int_t scores4 = rowsums_cpp(scores3);

  return Rcpp::wrap(scores4);

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
