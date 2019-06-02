#include <Rcpp.h>
#include <RcppThread.h>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <limits>
#include "types.h"
#include "utils-internal.h"

std::unordered_map<std::string, int> METRICS_enum = {

  /**
   * Instructions for adding a comparison metric:
   * - Add entry to METRICS_enum
   * - Write compare_*() function
   * - Add case entry to get_compare_ans()
   * - Add case entry to return_best_ans()
   * - Add case entry to return_best_ans_which()
   * - If zero values are not allowed, add case entry to fix_mot_bkg_zeros()
   * - Add case entry to compare_columns_cpp()
   */

  /* distance */      // timings from last make_DBscores() run:

  // Euclidean distance
  {"EUCL",      1},   // 35.78 s
  // Kullback-Leibler divergence
  {"KL",        2},   // 21.99 s
  // Hellinger distance
  {"HELL",      3},   // 17.37 s
  // Itakura-Saito distance
  {"IS",        4},   // 18.41 s
  // Squared Euclidean distance
  {"SEUCL",     5},   // 14.93 s
  // Manhattan distance
  {"MAN",       6},   // 14.92 s

  /* similarity */

  // Pearson correlation coefficient
  {"PCC",       7},   // 1.231 m
  // Sandelin-Wasserman similarity (aka SSD)
  {"SW",        8},   // 15.42 s
  // Average log-likelihood ratio
  {"ALLR",      9},   // 59.88 s
  // Bhattacharyya coefficient
  {"BHAT",     10},   // 15.67 s
  // Lower limit average log-likelihood ratio
  {"ALLR_LL",  11},

};

enum COMPARE_METRICS {

  EUCL    =  1,
  KL      =  2,
  HELL    =  3,
  IS      =  4,
  SEUCL   =  5,
  MAN     =  6,
  PCC     =  7,
  SW      =  8,
  ALLR    =  9,
  BHAT    = 10,
  ALLR_LL = 11

};

std::unordered_map<std::string, int> SCORESTRAT_enum = {

  /* Possible means to add: harmonic mean, weighted means */

  /* Weighted means would have to be handled differently for similarity vs
   * distance metrics
   */

  {"sum",     1},
  {"a.mean",  2},
  {"g.mean",  3},
  {"median",  4},
  {"wa.mean", 5},
  {"wg.mean", 6}

};

enum ADDSCORES_STRAT {

  SUM    = 1,
  AMEAN  = 2,
  GMEAN  = 3,
  MEDIAN = 4,
  WAMEAN = 5,
  WGMEAN = 6

};

std::unordered_map<std::string, int> STRATS_enum = {

  {"normal",   1},
  {"logistic", 2},
  {"weibull",  3}

};

enum DIST_FUNS {

  NORMAL   = 1,
  LOGISTIC = 2,
  WEIBULL  = 3

};

double pval_calculator(const double score, const double paramA,
    const double paramB, const int ltail, const str_t &dist) {

  switch (::STRATS_enum[dist]) {

    case NORMAL  : return R::pnorm(score, paramA, paramB, ltail, 1);
    case LOGISTIC: return R::plogis(score, paramA, paramB, ltail, 1);
    case WEIBULL : return R::pweibull(score, paramA, paramB, ltail, 1);

    default      : return -1.0;

  }

}

double score_sum (const vec_num_t &scores) {

  return std::accumulate(scores.begin(), scores.end(), 0.0);

}

double score_amean(const vec_num_t &scores, const int n) {

  return std::accumulate(scores.begin(), scores.end(), 0.0) / double(n);

}

double score_gmean(const vec_num_t &scores) {

  double sum = 0.0;
  for (std::size_t i = 0; i < scores.size(); ++i) {
    if (scores[i] > 0) sum += log(scores[i]);
  }

  if (sum == 0.0)
    return 0.0;
  else
    return exp(sum / double(scores.size()));

}

double score_median(vec_num_t scores) {

  if (scores.size() == 1) return scores[0];

  std::sort(scores.begin(), scores.end());

  if (scores.size() % 2 == 0)
    return (scores[scores.size() / 2 - 1] + scores[scores.size() / 2]) / 2.0;

  else
    return scores[scores.size() / 2];

}

double score_wamean(vec_num_t scores, const vec_num_t &ic) {

  double total_ic = std::accumulate(ic.begin(), ic.end(), 0.0);

  for (std::size_t i = 0; i < scores.size(); ++i) {
    scores[i] *= ic[i] / total_ic;
  }

  return std::accumulate(scores.begin(), scores.end(), 0.0);

}

double score_wgmean(vec_num_t scores, const vec_num_t &ic) {

  for (std::size_t i = 0; i < scores.size(); ++i) {
    if (scores[i] > 0.0) scores[i] = log(scores[i]);
  }

  if (std::accumulate(scores.begin(), scores.end(), 0.0) == 0.0)
    return 0.0;
  else
    return exp(score_wamean(scores, ic));

}

vec_num_t keep_good(const vec_num_t &scores, const vec_bool_t &good, const int n) {

  vec_num_t out;
  out.reserve(n);

  for (std::size_t i = 0; i < scores.size(); ++i) {
    if (good[i]) out.push_back(scores[i]);
  }

  return out;

}

vec_num_t combine_good_ic(const vec_num_t &ic1, const vec_num_t &ic2,
    const vec_bool_t &good, const int n) {

  vec_num_t final_ic;
  final_ic.reserve(n);

  for (std::size_t i = 0; i < good.size(); ++i) {
    if (good[i])
      final_ic.push_back(ic1[i] + ic2[i]);
  }

  return final_ic;

}

double calc_final_score(const vec_num_t &scores, const str_t &strat,
    const int n, const vec_bool_t &good, const vec_num_t &ic1,
    const vec_num_t &ic2) {

  switch (::SCORESTRAT_enum[strat]) {

    case SUM   : return score_sum(scores);
    case AMEAN : return score_amean(scores, n);
    case GMEAN : return score_gmean(keep_good(scores, good, n));
    case MEDIAN: return score_median(keep_good(scores, good, n));
    case WAMEAN: return score_wamean(keep_good(scores, good, n),
                                     combine_good_ic(ic1, ic2, good, n));
    case WGMEAN: return score_wgmean(keep_good(scores, good, n),
                                     combine_good_ic(ic1, ic2, good, n));

    default: return -333.333;

  }

}

double compare_hell(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();

  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  vec_num_t ans(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        ans[i] += pow(sqrt(mot1[i][j]) - sqrt(mot2[i][j]), 2.0);
      }
      ans[i] = sqrt(ans[i]) / sqrt(2.0);
    }
  }

  return calc_final_score(ans, strat, n, good, ic1, ic2);

}

double compare_is(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();

  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  vec_num_t ans(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        ans[i] += mot1[i][j] / mot2[i][j] - log(mot1[i][j] / mot2[i][j]) - 1.0;
      }
    }
  }

  return calc_final_score(ans, strat, n, good, ic1, ic2);

}

double compare_seucl(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();

  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  vec_num_t ans(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        ans[i] += pow(mot1[i][j] - mot2[i][j], 2.0);
      }
    }
  }

  return calc_final_score(ans, strat, n, good, ic1, ic2);

}

double compare_man(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();

  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  vec_num_t ans(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        ans[i] += abs(mot1[i][j] - mot2[i][j]);
      }
    }
  }

  return calc_final_score(ans, strat, n, good, ic1, ic2);

}

double compare_bhat(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();

  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  vec_num_t ans(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        ans[i] += sqrt(mot1[i][j] * mot2[i][j]);
      }
    }
  }

  return calc_final_score(ans, strat, n, good, ic1, ic2);

}

double compare_eucl(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  list_num_t diffmat(ncol, vec_num_t(nrow, 0.0));
  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  vec_num_t ans(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        ans[i] += pow(mot1[i][j] - mot2[i][j], 2.0);
      }
      ans[i] = sqrt(ans[i]);
    }
  }

  return calc_final_score(ans, strat, n, good, ic1, ic2);

}

double compare_pcc(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  // expensive!

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  list_num_t mmat(ncol, vec_num_t(nrow, 0.0));
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        mmat[i][j] = mot1[i][j] * mot2[i][j];
      }
    }
  }

  list_num_t pmat1(ncol, vec_num_t(nrow, 0.0));
  list_num_t pmat2(ncol, vec_num_t(nrow, 0.0));
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        pmat1[i][j] = pow(mot1[i][j], 2.0);
        pmat2[i][j] = pow(mot2[i][j], 2.0);
      }
    }
  }

  vec_num_t top(ncol);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      top[i] = double(nrow) * std::accumulate(mmat[i].begin(), mmat[i].end(), 0.0);
      top[i] -= std::accumulate(mot1[i].begin(), mot1[i].end(), 0.0)
                * std::accumulate(mot2[i].begin(), mot2[i].end(), 0.0);
    }
  }

  vec_num_t bot(ncol);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      bot[i] = double(nrow) * std::accumulate(pmat1[i].begin(), pmat1[i].end(), 0.0)
               - pow(std::accumulate(mot1[i].begin(), mot1[i].end(), 0.0), 2.0);
      bot[i] *= double(nrow) * std::accumulate(pmat2[i].begin(), pmat2[i].end(), 0.0)
                - pow(std::accumulate(mot2[i].begin(), mot2[i].end(), 0.0), 2.0);
      bot[i] = sqrt(bot[i]);
    }
  }

  vec_num_t topbot(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      /* Note:
       *   If both columns are uniform (i.e. {0.25, 0.25, 0.25, 0.25}) then
       *   bot = 0, leading to top / 0 = Inf. Since both being uniform also
       *   means they are identical, a perfect score is given (1.0).
       */
      topbot[i] = bot[i] == 0 ? 1.0 : top[i] / bot[i];
    }
  }
  return calc_final_score(topbot, strat, n, good, ic1, ic2);

}

double compare_kl(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  vec_num_t ans(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        ans[i] += mot1[i][j] * log(mot1[i][j] / mot2[i][j]);
        ans[i] += mot2[i][j] * log(mot2[i][j] / mot1[i][j]);
      }
      ans[i] *= 0.5;
    }
  }

  return calc_final_score(ans, strat, n, good, ic1, ic2);

}

double compare_allr(const list_num_t &mot1, const list_num_t &mot2,
    const vec_num_t &bkg1, const vec_num_t &bkg2, const double nsites1,
    const double nsites2, const str_t &strat, const vec_num_t &ic1,
    const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  vec_bool_t good(ncol, false);
  int counter = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++counter;
    }
  }

  list_num_t left(ncol, vec_num_t(nrow, 0.0)), right(ncol, vec_num_t(nrow, 0.0));
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        left[i][j] = (mot2[i][j] * nsites2) * log(mot1[i][j] / bkg1[j]);
        right[i][j] = (mot1[i][j] * nsites1) * log(mot2[i][j] / bkg2[j]);
      }
    }
  }

  vec_num_t answers(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      answers[i] += std::accumulate(left[i].begin(), left[i].end(), 0.0);
      answers[i] += std::accumulate(right[i].begin(), right[i].end(), 0.0);
      answers[i] /= nsites1 + nsites2;
    }
  }

  return calc_final_score(answers, strat, counter, good, ic1, ic2);

}

double compare_allr_ll(const list_num_t &mot1, const list_num_t &mot2,
    const vec_num_t &bkg1, const vec_num_t &bkg2, const double nsites1,
    const double nsites2, const str_t &strat, const vec_num_t &ic1,
    const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  vec_bool_t good(ncol, false);
  int counter = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++counter;
    }
  }

  list_num_t left(ncol, vec_num_t(nrow, 0.0)), right(ncol, vec_num_t(nrow, 0.0));
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        left[i][j] = (mot2[i][j] * nsites2) * log(mot1[i][j] / bkg1[j]);
        right[i][j] = (mot1[i][j] * nsites1) * log(mot2[i][j] / bkg2[j]);
      }
    }
  }

  vec_num_t answers(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      answers[i] += std::accumulate(left[i].begin(), left[i].end(), 0.0);
      answers[i] += std::accumulate(right[i].begin(), right[i].end(), 0.0);
      answers[i] /= nsites1 + nsites2;
      if (answers[i] < -2.0) answers[i] = -2.0;
    }
  }

  return calc_final_score(answers, strat, counter, good, ic1, ic2);

}

double compare_sw(const list_num_t &mot1, const list_num_t &mot2,
    const str_t &strat, const vec_num_t &ic1, const vec_num_t &ic2) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  vec_bool_t good(ncol, false);
  int n = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      good[i] = true;
      ++n;
    }
  }

  vec_num_t ans(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        ans[i] += pow(mot1[i][j] - mot2[i][j], 2.0);
      }
      ans[i] = 2.0 - ans[i];
    }
  }

  return calc_final_score(ans, strat, n, good, ic1, ic2);

}

void klfix(list_num_t &mot) {
  for (std::size_t i = 0; i < mot.size(); ++i) {
    for (std::size_t j = 0; j < mot[0].size(); ++j) {
      mot[i][j] += 0.01;
    }
  }
  return;
}

void bkgfix(vec_num_t &bkg) {
  bool fix = false;
  for (std::size_t i = 0; i < bkg.size(); ++i) {
    if (bkg[i] == 0) {
      fix = true;
      break;
    }
  }
  if (fix) {
    for (std::size_t i = 0; i < bkg.size(); ++i) {
      bkg[i] += 0.01 * (1.0 / double(bkg.size()));
    }
  }
  return;
}

void equalize_mot_cols(list_num_t &mot1, list_num_t &mot2,
    vec_num_t &ic1, vec_num_t &ic2, const double overlap) {

  std::size_t nrow = mot1[0].size();
  std::size_t ncol1 = mot1.size();
  std::size_t ncol2 = mot2.size();
  std::size_t overlap1 = overlap, overlap2 = overlap;

  if (overlap < 1) {
    overlap1 = overlap * ncol1;
    overlap2 = overlap * ncol2;
  }

  std::size_t ncol1_toadd = overlap1 > ncol2 ? 0 : ncol2 - overlap1;
  std::size_t ncol2_toadd = overlap2 > ncol1 ? 0 : ncol1 - overlap2;

  if (ncol1_toadd == 0 || ncol2_toadd == 0) return;

  if (ncol2 > ncol1) {

    list_num_t newmot(ncol1_toadd * 2 + ncol1, vec_num_t(nrow, -1.0));
    vec_num_t newic(ncol1_toadd * 2 + ncol1, 0.0);

    for (std::size_t i = ncol1_toadd; i < ncol1_toadd + ncol1; ++i) {
      newmot[i].assign(mot1[i - ncol1_toadd].begin(), mot1[i - ncol1_toadd].end());
    }

    for (std::size_t i = ncol1_toadd; i < ncol1_toadd + ncol1; ++i) {
      newic[i] = ic1[i - ncol1_toadd];
    }

    mot1.assign(newmot.begin(), newmot.end());
    ic1.assign(newic.begin(), newic.end());

  } else {

    list_num_t newmot(ncol2_toadd * 2 + ncol2, vec_num_t(nrow, -1.0));
    vec_num_t newic(ncol2_toadd * 2 + ncol2, 0.0);

    for (std::size_t i = ncol2_toadd; i < ncol2_toadd + ncol2; ++i) {
      newmot[i].assign(mot2[i - ncol2_toadd].begin(), mot2[i - ncol2_toadd].end());
    }

    for (std::size_t i = ncol2_toadd; i < ncol2_toadd + ncol2; ++i) {
      newic[i] = ic2[i - ncol2_toadd];
    }

    mot2.assign(newmot.begin(), newmot.end());
    ic2.assign(newic.begin(), newic.end());

  }

  return;

}

int get_alignlen(const list_num_t &mot1, const list_num_t &mot2) {

  int out = 0;
  for (std::size_t i = 0; i < mot1.size(); ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) ++out;
  }

  return out;

}

double calc_mic(const vec_num_t &tic) {

  double out = 0.0;
  int counter = 0;

  for (std::size_t i = 0; i < tic.size(); ++i) {
    if (tic[i] >= 0) {
      out += tic[i];
      ++counter;
    }
  }

  return out / double(counter);

}

list_num_t get_motif_rc(const list_num_t &mot) {

  list_num_t rcmot = mot;
  std::reverse(rcmot.begin(), rcmot.end());
  for (std::size_t i = 0; i < rcmot.size(); ++i) {
    std::reverse(rcmot[i].begin(), rcmot[i].end());
  }

  return rcmot;

}

void fix_lowic_pos(list_num_t &tmot1, list_num_t &tmot2,
    vec_num_t &tic1, vec_num_t &tic2, const double posic) {

  for (std::size_t i = 0; i < tmot1.size(); ++i) {
    if (tic1[i] < posic) {
      for (std::size_t j = 0; j < tmot1[0].size(); ++j) {
        tmot1[i][j] = -1.0;
      }
      tic1[i] = -1.0;
    }
    if (tic2[i] < posic) {
      for (std::size_t j = 0; j < tmot1[0].size(); ++j) {
        tmot2[i][j] = -1.0;
      }
      tic2[i] = -1.0;
    }
  }

  return;

}

void get_compare_ans(vec_num_t &ans, const std::size_t i,
    const list_num_t &tmot1, const list_num_t &tmot2, const bool lowic,
    const bool norm, const int alignlen, const std::size_t tlen,
    const std::string &method, const double nsites1, const double nsites2,
    const vec_num_t &bkg1, const vec_num_t &bkg2, const str_t &strat,
    const vec_num_t &ic1, const vec_num_t &ic2) {

  switch (::METRICS_enum[method]) {

    case EUCL   : ans[i] = lowic ? std::numeric_limits<double>::max()
                                 : compare_eucl(tmot1, tmot2, strat, ic1, ic2)
                                   * double(tlen) / double(alignlen);
                  break;
    case KL     : ans[i] = lowic ? std::numeric_limits<double>::max()
                                 : compare_kl(tmot1, tmot2, strat, ic1, ic2)
                                   * double(tlen) / double(alignlen);
                  break;
    case PCC    : ans[i] = lowic ? -std::numeric_limits<double>::max()
                                 : compare_pcc(tmot1, tmot2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;
    case SW     : ans[i] = lowic ? -std::numeric_limits<double>::max()
                                 : compare_sw(tmot1, tmot2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;
    case ALLR   : ans[i] = lowic ? -std::numeric_limits<double>::max()
                                 : compare_allr(tmot1, tmot2, bkg1,
                                   bkg2, nsites1, nsites2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;
    case BHAT   : ans[i] = lowic ? -std::numeric_limits<double>::max()
                                 : compare_bhat(tmot1, tmot2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;
    case HELL   : ans[i] = lowic ? std::numeric_limits<double>::max()
                                 : compare_hell(tmot1, tmot2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;
    case IS     : ans[i] = lowic ? std::numeric_limits<double>::max()
                                 : compare_is(tmot1, tmot2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;
    case SEUCL  : ans[i] = lowic ? std::numeric_limits<double>::max()
                                 : compare_seucl(tmot1, tmot2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;
    case MAN    : ans[i] = lowic ? std::numeric_limits<double>::max()
                                 : compare_man(tmot1, tmot2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;
    case ALLR_LL: ans[i] = lowic ? -std::numeric_limits<double>::max()
                                 : compare_allr_ll(tmot1, tmot2, bkg1,
                                     bkg2, nsites1, nsites2, strat, ic1, ic2)
                                   * double(alignlen) / double(tlen);
                  break;

  }

  return;

}

double return_best_ans(const vec_num_t &ans, const std::string &method) {

  switch (::METRICS_enum[method]) {

    case EUCL   : return *std::min_element(ans.begin(), ans.end());
    case KL     : return *std::min_element(ans.begin(), ans.end());
    case HELL   : return *std::min_element(ans.begin(), ans.end());
    case IS     : return *std::min_element(ans.begin(), ans.end());
    case SEUCL  : return *std::min_element(ans.begin(), ans.end());
    case MAN    : return *std::min_element(ans.begin(), ans.end());
    case PCC    : return *std::max_element(ans.begin(), ans.end());
    case SW     : return *std::max_element(ans.begin(), ans.end());
    case ALLR   : return *std::max_element(ans.begin(), ans.end());
    case BHAT   : return *std::max_element(ans.begin(), ans.end());
    case ALLR_LL: return *std::max_element(ans.begin(), ans.end());

  }

  return -1111.1111;

}

double compare_motif_pair(list_num_t mot1, list_num_t mot2,
    const std::string method, const double moverlap, const bool RC,
    vec_num_t ic1, vec_num_t ic2, const double minic, const bool norm,
    const double posic, const vec_num_t &bkg1, const vec_num_t &bkg2,
    const double nsites1, const double nsites2, const str_t &strat) {

  double ans_rc = 0.0;
  if (RC) {
    list_num_t rcmot2 = get_motif_rc(mot2);
    vec_num_t rcic2 = ic2;
    std::reverse(rcic2.begin(), rcic2.end());
    ans_rc = compare_motif_pair(mot1, rcmot2, method, moverlap, false,
        ic1, rcic2, minic, norm, posic, bkg1, bkg2, nsites1, nsites2, strat);
  }

  std::size_t tlen = mot1.size() >= mot2.size() ? mot1.size() : mot2.size();

  equalize_mot_cols(mot1, mot2, ic1, ic2, moverlap);

  std::size_t minw = mot1.size() <= mot2.size() ? mot1.size() : mot2.size();
  list_num_t tmot1(minw), tmot2(minw);
  vec_num_t tic1(minw), tic2(minw);
  std::size_t fori = 1 + mot1.size() - minw, forj = 1 + mot2.size() - minw;
  vec_num_t ans;
  ans = RC ? vec_num_t(fori * forj + 1) : vec_num_t(fori * forj);
  int alignlen;
  double mic1, mic2;
  bool lowic = false;
  std::size_t counter = 0;

  for (std::size_t i = 0; i < fori; ++i) {
    for (std::size_t j = 0; j < forj; ++j) {

      for (std::size_t k = 0; k < minw; ++k) {
        tmot1[k] = mot1[k + i];
        tmot2[k] = mot2[k + j];
        tic1[k] = ic1[k + i];
        tic2[k] = ic2[k + j];
      }

      if (posic > 0) fix_lowic_pos(tmot1, tmot2, tic1, tic2, posic);

      alignlen = norm ? get_alignlen(tmot1, tmot2) : tlen;

      mic1 = calc_mic(tic1);
      mic2 = calc_mic(tic2);
      if (mic1 < minic || mic2 < minic) lowic = true;

      get_compare_ans(ans, counter, tmot1, tmot2, lowic, norm, alignlen, tlen,
          method, nsites1, nsites2, bkg1, bkg2, strat, tic1, tic2);

      ++counter;
      lowic = false;

    }
  }

  if (RC) ans[counter] = ans_rc;

  return return_best_ans(ans, method);

}

double internal_posIC(vec_num_t pos, const vec_num_t &bkg,
    const int type, bool relative) {

  if (type == 2) return std::accumulate(pos.begin(), pos.end(), 0.0);

  if (relative) {
    for (std::size_t i = 0; i < pos.size(); ++i) {
      double tmp = pos[i] / bkg[i];
      pos[i] *= tmp >= 0 ? log2(tmp) : 0.0;
      if (pos[i] < 0) pos[i] = 0.0;
    }
    return std::accumulate(pos.begin(), pos.end(), 0.0);
  }

  vec_num_t heights(pos.size());
  for (std::size_t i = 0; i < pos.size(); ++i) {
    heights[i] = pos[i] > 0 ? -pos[i] * log2(pos[i]) : 0.0;
  }
  return log2(pos.size()) - std::accumulate(heights.begin(), heights.end(), 0.0);

}

list_num_t get_merged_motif(const list_num_t &mot1, const list_num_t &mot2,
    const int weight) {

  /* TODO: Is it possible for the motif to be internally trimmed as a result
   *       of min.position.ic > 0?
   */

  list_num_t out;
  out.reserve(mot1.size());
  for (std::size_t i = 0; i < mot1.size(); ++i) {
    if (mot1[i][0] < 0 && mot2[i][0] >= 0) {
      out.push_back(mot2[i]);
    } else if (mot2[i][0] < 0 && mot1[i][0] >= 0) {
      out.push_back(mot1[i]);
    } else if (mot1[i][0] >= 0 && mot2[i][0] >= 0) {
      vec_num_t tmp(mot1[0].size());
      for (std::size_t j = 0; j < mot1[0].size(); ++j) {
        tmp[j] = (mot1[i][j] * double(weight) + mot2[i][j]) / (double(weight) + 1);
      }
      out.push_back(tmp);
    }
  }

  return out;

}

int return_best_ans_which(const vec_num_t &ans, const std::string &method) {

  switch (::METRICS_enum[method]) {

    case EUCL   : return std::distance(ans.begin(), std::min_element(ans.begin(), ans.end()));
    case KL     : return std::distance(ans.begin(), std::min_element(ans.begin(), ans.end()));
    case HELL   : return std::distance(ans.begin(), std::min_element(ans.begin(), ans.end()));
    case IS     : return std::distance(ans.begin(), std::min_element(ans.begin(), ans.end()));
    case SEUCL  : return std::distance(ans.begin(), std::min_element(ans.begin(), ans.end()));
    case MAN    : return std::distance(ans.begin(), std::min_element(ans.begin(), ans.end()));
    case PCC    : return std::distance(ans.begin(), std::max_element(ans.begin(), ans.end()));
    case SW     : return std::distance(ans.begin(), std::max_element(ans.begin(), ans.end()));
    case ALLR   : return std::distance(ans.begin(), std::max_element(ans.begin(), ans.end()));
    case BHAT   : return std::distance(ans.begin(), std::max_element(ans.begin(), ans.end()));
    case ALLR_LL: return std::distance(ans.begin(), std::max_element(ans.begin(), ans.end()));

  }

  return -1;

}


void merge_motif_pair_subworker(list_num_t mot1, list_num_t mot2,
    const std::string &method, const double minoverlap, vec_num_t ic1,
    vec_num_t ic2, const bool norm, const double posic, const double minic,
    double &score, int &offset, const double nsites1, const double nsites2,
    const vec_num_t &bkg1, const vec_num_t &bkg2, const str_t &strat) {

  std::size_t tlen = mot1.size() >= mot2.size() ? mot1.size() : mot2.size();

  equalize_mot_cols(mot1, mot2, ic1, ic2, minoverlap);

  std::size_t minw = mot1.size() <= mot2.size() ? mot1.size() : mot2.size();
  list_num_t tmot1(minw), tmot2(minw);
  vec_num_t tic1(minw), tic2(minw);
  std::size_t fori = 1 + mot1.size() - minw, forj = 1 + mot2.size() - minw;
  int alignlen;
  double mic1, mic2;
  bool lowic = false;
  vec_num_t ans(fori * forj);
  std::size_t counter = 0;

  for (std::size_t i = 0; i < fori; ++i) {
    for (std::size_t j = 0; j < forj; ++j) {

      for (std::size_t k = 0; k < minw; ++k) {
        tmot1[k] = mot1[k + i];
        tmot2[k] = mot2[k + j];
        tic1[k] = ic1[k + i];
        tic2[k] = ic2[k + j];
      }

      if (posic > 0) fix_lowic_pos(tmot1, tmot2, tic1, tic2, posic);

      alignlen = norm ? get_alignlen(tmot1, tmot2) : tlen;

      mic1 = calc_mic(tic1);
      mic2 = calc_mic(tic2);
      if (mic1 < minic || mic2 < minic) lowic = true;

      get_compare_ans(ans, counter, tmot1, tmot2, lowic, norm, alignlen, tlen,
          method, nsites1, nsites2, bkg1, bkg2, strat, tic1, tic2);

      ++counter;
      lowic = false;

    }
  }

  score = return_best_ans(ans, method);
  offset = return_best_ans_which(ans, method);

}

list_num_t add_motif_columns(const list_num_t &mot, const int tlen,
    const int add) {

  list_num_t out(tlen, vec_num_t(mot[0].size(), -1.0));
  int lenadd = add + int(mot.size());
  int counter = 0;
  for (int i = add; i < lenadd; ++i) {
    out[i] = mot[counter];
    ++counter;
  }

  return out;

}

void trim_both_motifs(list_num_t &m1, list_num_t &m2) {

  int mlen = m1.size();
  int tleft = 0;

  for (int i = 0; i < mlen; ++i) {
    if (m1[i][0] < 0 && m2[i][0] < 0) {
      --mlen;
      ++tleft;
    } else {
      break;
    }
  }

  for (int i = mlen - 1; i >= 0; --i) {
    if (m1[i][0] < 0 && m2[i][0] < 0) {
      --mlen;
    } else {
      break;
    }
  }

  if (mlen > 0) {

    list_num_t tm1(mlen), tm2(mlen);
    for(int i = tleft; i < tleft + mlen; ++i) {
      tm1[i - tleft] = m1[i];
      tm2[i - tleft] = m2[i];
    }

    m1 = tm1;
    m2 = tm2;

  }

}

list_num_t merge_motif_pair(list_num_t mot1, list_num_t mot2,
    const std::string &method, const double minoverlap, const bool RC,
    vec_num_t ic1, vec_num_t ic2, const int weight,
    const bool norm, const double posic, const double minic,
    const double nsites1, const double nsites2, const vec_num_t &bkg1,
    const vec_num_t &bkg2, const str_t &strat) {

  double score;
  int offset;

  merge_motif_pair_subworker(mot1, mot2, method, minoverlap, ic1, ic2, norm,
      posic, minic, score, offset, nsites1, nsites2, bkg1, bkg2, strat);

  if (RC) {
    double score_rc;
    int offset_rc;
    list_num_t rcmot2 = get_motif_rc(mot2);
    vec_num_t rcic2 = ic2;
    std::reverse(rcic2.begin(), rcic2.end());
    merge_motif_pair_subworker(mot1, rcmot2, method, minoverlap, ic1, rcic2,
        norm, posic, minic, score_rc, offset_rc, nsites1, nsites2, bkg1, bkg2, strat);
    if (score_rc > score) {
      offset = offset_rc;
      mot2 = rcmot2;
    }
  }

  equalize_mot_cols(mot1, mot2, ic1, ic2, minoverlap);

  if (mot1.size() > mot2.size()) {
    mot2 = add_motif_columns(mot2, mot1.size(),
        offset % mot1.size() - offset / mot1.size());
  } else if (mot2.size() > mot1.size()) {
    mot1 = add_motif_columns(mot1, mot2.size(),
        offset % mot2.size() - offset / mot2.size());
  }

  trim_both_motifs(mot1, mot2);
  list_num_t merged = get_merged_motif(mot1, mot2, weight);

  return merged;

}

vec_num_t merge_bkg_pair(const vec_num_t &bkg1, const vec_num_t &bkg2,
    const int weight) {

  vec_num_t newbkg(bkg1.size());

  for (std::size_t i = 0; i < bkg1.size(); ++i) {
    newbkg[i] = (bkg1[i] * double(weight) + bkg2[i]) / (double(weight) + 1);
  }

  return newbkg;

}

vec_num_t calc_ic_motif(const list_num_t &motif, const vec_num_t &bkg,
    const bool relative) {

  vec_num_t out(motif.size());
  for (std::size_t i = 0; i < motif.size(); ++i) {
    out[i] = internal_posIC(motif[i], bkg, 1, relative);
  }

  return out;

}

void find_offsets(list_num_t mot1, list_num_t mot2,
    bool &use_rc, const std::string &method, const double minoverlap,
    vec_num_t ic1, vec_num_t ic2, const bool norm, const double posic,
    const double minic, const bool RC, const double nsites1, const double nsites2,
    const vec_num_t &bkg1, const vec_num_t &bkg2, int &offset, const str_t &strat) {

  double score;

  merge_motif_pair_subworker(mot1, mot2, method, minoverlap, ic1, ic2, norm,
      posic, minic, score, offset, nsites1, nsites2, bkg1, bkg2, strat);

  if (RC) {
    double score_rc;
    int offset_rc;
    list_num_t rcmot2 = get_motif_rc(mot2);
    vec_num_t rcic2 = ic2;
    std::reverse(rcic2.begin(), rcic2.end());
    merge_motif_pair_subworker(mot1, rcmot2, method, minoverlap, ic1, rcic2,
        norm, posic, minic, score_rc, offset_rc, nsites1, nsites2, bkg1, bkg2, strat);
    if (score_rc > score) {
      offset = offset_rc;
      mot2 = rcmot2;
      use_rc = true;
    }
  }

  return;

}

void neg_one_to_zero(list_num_t &mot) {

  for (std::size_t i = 0; i < mot.size(); ++i) {
    for (std::size_t j = 0; j < mot[0].size(); ++j) {
      if (mot[i][j] < 0) mot[i][j] = 0;
    }
  }

  return;

}

void fix_mot_bkg_zeros(list_num_t &mot, vec_num_t &bkg, const std::string &method) {

  switch (::METRICS_enum[method]) {
    case KL:
    case IS:
    case ALLR:
    case ALLR_LL: klfix(mot);
                  bkgfix(bkg);
  }

}

int count_left_empty(const list_num_t &m) {

  int left = 0;
  for (std::size_t i = 0; i < m.size(); ++i) {
    if (m[i][0] < 0)
      ++left;
    else
      break;
  }

  return left;

}

/* C++ ENTRY ---------------------------------------------------------------- */

// [[Rcpp::export(rng = false)]]
std::vector<double> compare_motifs_cpp(const Rcpp::List &mots,
    const std::vector<int> &index1, const std::vector<int> &index2,
    const std::string &method, double minoverlap, const bool RC,
    std::vector<std::vector<double>> &bkg, const int type, const bool relative,
    const double minic, const bool norm, const int nthreads, const double posic,
    const std::vector<double> &nsites, const std::string &strat) {

  /* compare motifs by indices, i.e. mots[index1[i]] vs mots[index2[i]] */

  if (minoverlap < 0) minoverlap = 1;

  if (type != 1 && type != 2)
    Rcpp::stop("type must be 1 or 2");
  if (minic < 0)
    Rcpp::stop("min.mean.ic must be positive");
  if (posic < 0)
    Rcpp::stop("min.position.ic must be positive");
  if (mots.size() == 0)
    Rcpp::stop("empty motif list");
  if (bkg.size() == 0)
    Rcpp::stop("empty bkg list");
  if (int(mots.size()) != int(bkg.size()))
    Rcpp::stop("different motif and bkg lengths");
  if (index1.size() != index2.size())
    Rcpp::stop("lengths of indices do not match [compare_motifs_cpp()]");

  list_nmat_t vmots(bkg.size());
  list_num_t icscores(vmots.size());

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    Rcpp::NumericMatrix tmp = mots(i);
    vmots[i] = R_to_cpp_motif_num(tmp);
    if (vmots[i].size() == 0)
      Rcpp::stop("encountered an empty motif [compare_motifs_cpp()]");
    fix_mot_bkg_zeros(vmots[i], bkg[i], method);
  }

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    icscores[i].reserve(vmots[i].size());
    for (std::size_t j = 0; j < vmots[i].size(); ++j) {
      icscores[i].push_back(internal_posIC(vmots[i][j], bkg[i], type, relative));
    }
  }

  vec_num_t answers(index1.size());
  RcppThread::parallelFor(0, answers.size(),
      [&answers, &vmots, &index1, &index2, &icscores, &method, minoverlap, RC,
        minic, norm, posic, &bkg, &nsites, &strat]
      (std::size_t i) {
        answers[i] = compare_motif_pair(vmots[index1[i]], vmots[index2[i]],
            method, minoverlap, RC, icscores[index1[i]], icscores[index2[i]],
            minic, norm, posic, bkg[index1[i]], bkg[index1[i]], nsites[index1[i]],
            nsites[index2[i]], strat);
      }, nthreads);

  return answers;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::vector<double>> compare_motifs_all_cpp(const Rcpp::List &mots,
    const std::string &method, double minoverlap, const bool RC,
    std::vector<std::vector<double>> &bkg, const int type, const bool relative,
    const double minic, const bool norm, const int nthreads, const double posic,
    const std::vector<double> &nsites, const std::string &strat) {

  /* compare all motifs to all motifs (without comparing the same motifs twice) */

  if (minoverlap < 0) minoverlap = 1;

  if (type != 1 && type != 2)
    Rcpp::stop("type must be 1 or 2");
  if (minic < 0)
    Rcpp::stop("min.mean.ic must be positive");
  if (posic < 0)
    Rcpp::stop("min.position.ic must be positive");
  if (mots.size() == 0)
    Rcpp::stop("empty motif list");
  if (bkg.size() == 0)
    Rcpp::stop("empty bkg list");
  if (int(mots.size()) != int(bkg.size()))
    Rcpp::stop("different motif and bkg lengths");

  list_nmat_t vmots(bkg.size());
  list_num_t icscores(vmots.size());

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    Rcpp::NumericMatrix tmp = mots(i);
    vmots[i] = R_to_cpp_motif_num(tmp);
    if (vmots[i].size() == 0)
      Rcpp::stop("encountered an empty motif [compare_motifs_all_cpp()]");
    fix_mot_bkg_zeros(vmots[i], bkg[i], method);
  }

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    icscores[i].reserve(vmots[i].size());
    for (std::size_t j = 0; j < vmots[i].size(); ++j) {
      icscores[i].push_back(internal_posIC(vmots[i][j], bkg[i], type, relative));
    }
  }

  list_num_t answers(vmots.size());
  RcppThread::parallelFor(0, answers.size(),
      [&answers, &vmots, &icscores, &method, minoverlap, RC, minic, norm, posic,
      &bkg, &nsites, &strat]
      (std::size_t i) {
        answers[i].reserve(vmots.size() - i);
        for (std::size_t j = i; j < vmots.size(); ++j) {
          answers[i].push_back(compare_motif_pair(vmots[i], vmots[j], method,
              minoverlap, RC, icscores[i], icscores[j], minic, norm, posic,
              bkg[i], bkg[j], nsites[i], nsites[j], strat));
        }
      }, nthreads);

  return answers;

}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_comparison_matrix(const std::vector<double> &ans,
    const std::vector<int> &index1, const std::vector<int> &index2,
    const std::string &method, const Rcpp::StringVector &motnames) {

  /* convert the output from compare_motifs_all_cpp() to a matrix */

  Rcpp::NumericMatrix out(motnames.size(), motnames.size());

  for (std::size_t i = 0; i < ans.size(); ++i) {
    out(index1[i], index2[i]) = ans[i];
    out(index2[i], index1[i]) = ans[i];
  }

  Rcpp::rownames(out) = motnames;
  Rcpp::colnames(out) = motnames;

  return out;

}

// [[Rcpp::export(rng = false)]]
Rcpp::List view_motifs_prep(const Rcpp::List &mots, const std::string &method,
    const bool RC, double minoverlap, const double minic, const double posic,
    std::vector<std::vector<double>> &bkg, const bool relative, const bool norm,
    const Rcpp::StringVector &rnames, const std::vector<double> &nsites,
    const std::string &strat) {

  if (minoverlap < 0) minoverlap = 1;
  if (minic < 0)
    Rcpp::stop("min.mean.ic must be positive");
  if (posic < 0)
    Rcpp::stop("min.position.ic must be positive");
  if (mots.size() == 0)
    Rcpp::stop("empty motif list");

  list_nmat_t vmots(mots.size());
  list_num_t icscores(vmots.size());
  vec_int_t motlens(mots.size());

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    Rcpp::NumericMatrix tmp = mots(i);
    vmots[i] = R_to_cpp_motif_num(tmp);
    if (vmots[i].size() == 0)
      Rcpp::stop("encountered an empty motif [compare_motifs_all_cpp()]");
    motlens[i] = vmots[i].size();
    fix_mot_bkg_zeros(vmots[i], bkg[i], method);
  }

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    icscores[i] = calc_ic_motif(vmots[i], bkg[i], relative);
  }

  vec_bool_t which_rc(vmots.size() - 1, false);
  vec_int_t offsets(vmots.size() - 1, 0);

  bool rc_ = false;
  for (std::size_t i = 1; i < vmots.size(); ++i) {

    rc_ = false;

    find_offsets(vmots[0], vmots[i],  rc_, method, minoverlap,
        icscores[0], icscores[i], norm, posic, minic, RC, nsites[0], nsites[i],
        bkg[0], bkg[i], offsets[i - 1], strat);

    which_rc[i - 1] = rc_;

  }

  for (std::size_t i = 0; i < which_rc.size(); ++i) {
    if (which_rc[i]) vmots[i + 1] = get_motif_rc(vmots[i + 1]);
  }

  vec_int_t toadd(vmots.size() - 1);
  list_nmat_t ttmots(vmots.size() - 1);

  for (std::size_t i = 1; i < vmots.size(); ++i) {

    list_num_t tmot1 = vmots[0], tmot2 = vmots[i];
    vec_num_t tic1 = icscores[0], tic2 = icscores[i];

    equalize_mot_cols(tmot1, tmot2, tic1, tic2, minoverlap);

    if (tmot1.size() > tmot2.size()) {
      tmot2 = add_motif_columns(tmot2, tmot1.size(),
          offsets[i - 1] % tmot1.size() - offsets[i - 1] / tmot1.size());
    } else if (tmot2.size() > tmot1.size()) {
      tmot1 = add_motif_columns(tmot1, tmot2.size(),
          offsets[i - 1] % tmot2.size() - offsets[i - 1] / tmot2.size());
    }

    trim_both_motifs(tmot1, tmot2);

    int left = count_left_empty(tmot1);
    toadd[i - 1] = left;

    ttmots[i - 1] = tmot2;

  }

  list_num_t mmot1 = vmots[0];

  int maxadd = *std::max_element(toadd.begin(), toadd.end());

  if (maxadd > 0) {
    mmot1 = add_motif_columns(mmot1, mmot1.size() + maxadd, maxadd);
  }

  neg_one_to_zero(mmot1);

  Rcpp::NumericMatrix mot1 = cpp_to_R_motif(mmot1);
  Rcpp::rownames(mot1) = rnames;

  Rcpp::List motlist(vmots.size());
  motlist[0] = mot1;
  for (std::size_t i = 0; i < vmots.size() - 1; ++i) {

    list_num_t tttmot = ttmots[i];
    if (maxadd - toadd[i] > 0) {
      tttmot = add_motif_columns(tttmot, tttmot.size() + maxadd - toadd[i],
          maxadd - toadd[i]);
    }

    neg_one_to_zero(tttmot);

    Rcpp::NumericMatrix tttmot2 = cpp_to_R_motif(tttmot);
    Rcpp::rownames(tttmot2) = rnames;

    motlist[i + 1] = tttmot2;

  }

  return Rcpp::List::create(
        Rcpp::_["motifs"] = motlist,
        Rcpp::_["motIsRC"] = which_rc
      );

}

// [[Rcpp::export(rng = false)]]
Rcpp::List merge_motifs_cpp(const Rcpp::List &mots,
    const std::string &method, const bool RC, double minoverlap, const double minic,
    const double posic, std::vector<std::vector<double>> &bkg,
    const bool relative, const bool norm, const std::vector<double> &nsites,
    const std::string &strat) {

  /* merge a list of motifs, as well as their backgrounds */

  if (minoverlap < 0) minoverlap = 1;
  if (minic < 0)
    Rcpp::stop("min.mean.ic must be positive");
  if (posic < 0)
    Rcpp::stop("min.position.ic must be positive");
  if (mots.size() == 0)
    Rcpp::stop("empty motif list");

  list_nmat_t vmots(mots.size());
  list_num_t icscores(vmots.size());

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    Rcpp::NumericMatrix tmp = mots(i);
    vmots[i] = R_to_cpp_motif_num(tmp);
    if (vmots[i].size() == 0)
      Rcpp::stop("encountered an empty motif [compare_motifs_all_cpp()]");
    fix_mot_bkg_zeros(vmots[i], bkg[i], method);
  }

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    icscores[i].reserve(vmots[i].size());
    for (std::size_t j = 0; j < vmots[i].size(); ++j) {
      icscores[i].push_back(internal_posIC(vmots[i][j], bkg[i], 1, relative));
    }
  }

  int weight = 1;
  // double nnsites = nsites[0] > nsites[1] ? nsites[0] : nsites[1];
  double nnsites = nsites[0] + nsites[1];
  list_num_t merged = merge_motif_pair(vmots[0], vmots[1], method, minoverlap,
      RC, icscores[0], icscores[1], weight, norm, posic, minic, nsites[0],
      nsites[1], bkg[0], bkg[1], strat);
  vec_num_t mergedbkg = merge_bkg_pair(bkg[0], bkg[1], weight);
  vec_num_t mergedic = calc_ic_motif(merged, mergedbkg, relative);

  if (vmots.size() > 2) {
    for (std::size_t i = 2; i < vmots.size(); ++i) {
      ++weight;
      merged = merge_motif_pair(merged, vmots[i], method, minoverlap, RC,
          mergedic, icscores[i], weight, norm, posic, minic, nnsites,
          nsites[i], mergedbkg, bkg[i], strat);
      mergedbkg = merge_bkg_pair(mergedbkg, bkg[i], weight);
      mergedic = calc_ic_motif(merged, mergedbkg, relative);
      nnsites += nsites[i];
      // nnsites = nnsites > nsites[i] ? nnsites : nsites[i];
    }
  }

  Rcpp::NumericMatrix outmotif = cpp_to_R_motif(merged);
  Rcpp::NumericVector outbkg = Rcpp::wrap(mergedbkg);

  return Rcpp::List::create(outmotif, outbkg);

}

// [[Rcpp::export(rng = false)]]
double compare_columns_cpp(const std::vector<double> &p1,
    const std::vector<double> &p2, const std::vector<double> &b1,
    const std::vector<double> &b2, const double n1 = 100,
    const double n2 = 100, const std::string &m = "PCC") {

  if (p1.size() < 2) Rcpp::stop("columns should have at least 2 entries");
  if (p1.size() != p2.size()) Rcpp::stop("both columns must be equal in size");

  list_num_t pp1(1, p1), pp2(1, p2);

  double ans;

  switch (::METRICS_enum[m]) {
    case EUCL   : ans = compare_eucl(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case KL     : ans = compare_kl(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case PCC    : ans = compare_pcc(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case SW     : ans = compare_sw(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case ALLR   : if (b1.size() != p1.size() || b2.size() != p1.size())
                    Rcpp::stop("incorrect background vector length");
                  if (n1 <= 1 || n2 <= 1)
                    Rcpp::stop("nsites1/nsites2 should be greater than 1");
                  ans = compare_allr(pp1, pp2, b1, b2, n1, n2, "sum", vec_num_t(),
                      vec_num_t());
                  break;
    case BHAT   : ans = compare_bhat(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case HELL   : ans = compare_hell(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case IS     : ans = compare_is(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case SEUCL  : ans = compare_seucl(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case MAN    : ans = compare_man(pp1, pp2, "sum", vec_num_t(), vec_num_t());
                  break;
    case ALLR_LL: if (b1.size() != p1.size() || b2.size() != p1.size())
                    Rcpp::stop("incorrect background vector length");
                  if (n1 <= 1 || n2 <= 1)
                    Rcpp::stop("nsites1/nsites2 should be greater than 1");
                  ans = compare_allr_ll(pp1, pp2, b1, b2, n1, n2, "sum", vec_num_t(),
                      vec_num_t());
                  break;
    default     : Rcpp::stop("unknown metric");
  }

  return ans;

}

// [[Rcpp::export(rng = false)]]
std::vector<double> pval_extractor(const std::vector<int> &ncols,
    const std::vector<double> &scores, const std::vector<int> &indices1,
    const std::vector<int> &indices2, const std::string &method,
    const std::vector<int> &subject, const std::vector<int> &target,
    const std::vector<double> &paramA, const std::vector<double> &paramB,
    const std::vector<std::string> &distribution, const int nthreads) {

  int ltail = 1;
  switch (::METRICS_enum[method]) {
    case PCC:
    case SW:
    case ALLR:
    case BHAT:
    case ALLR_LL: ltail = 0;
  }

  vec_num_t pvals(scores.size(), 0.0);

  std::size_t n = target.size() - 1;

  RcppThread::parallelFor(0, scores.size(),
      [&scores, &n, &pvals, &ncols, &indices1, &indices2, &subject, &target,
      &paramA, &paramB, &distribution, ltail]
      (std::size_t i) {

        int m1, m2, n1, n2;
        int row = -1;
        bool ok = false;

        /* Some notes:
         * - If the ncol (subject/target) combination is missing in db.scores,
         *   then +1 is added to each ncol in hopes to find the next possible
         *   combination. If it fails, then the P-value calculation is skipped
         *   and the final logPvalue is 0.
         */

        if (abs(scores[i]) != std::numeric_limits<double>::max()) {

          ok = false;

          m1 = ncols[indices1[i]];
          m2 = ncols[indices2[i]];

          n1 = std::min(m1, m2);
          n2 = std::max(m1, m2);

          if (n1 < subject[0])
            n1 = subject[0];
          else if (n1 > subject[n])
            n1 = subject[n];

          if (n2 < target[0])
            n2 = target[0];
          else if (n2 > target[n])
            n2 = target[n];

          // RcppThread::Rcout << n1 << ' ' << n2 << '\n';
          while (!ok) {
            for (std::size_t j = 0; j < subject.size(); ++j) {
          // RcppThread::Rcout << subject[j] << ' ' << target[j] << '\n';
              if (n1 == subject[j] && n2 == target[j]) {
                row = int(j);
                ok = true;
                break;
              }
            }
            ++n1;
            ++n2;
            if (!ok && (n1 > int(subject[n]) || n2 > int(target[n]))) {
              /* indicate failure to find right row */
              row = -1;
            }
          }
          // RcppThread::Rcout << row << '\n';

          if (row != -1) {
            /* only calculate pval if row was found */
            pvals[i] = pval_calculator(scores[i], paramA[row], paramB[row], ltail,
                distribution[row]);
          }

        } else {

          /* what to do when score is NA (represented as max double) */

          pvals[i] = 0.0;

        }

      }, nthreads);

  return pvals;

}
