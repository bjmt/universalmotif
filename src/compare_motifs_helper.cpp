#include <Rcpp.h>
#include <RcppThread.h>
#include <cmath>
#include "types.h"
#include "utils-internal.h"

std::unordered_map<std::string, int> metrics_enum = {
  {"EUCL",  1},
  {"MEUCL", 2},
  {"PCC",   3},
  {"MPCC",  4},
  {"KL",    5},
  {"MKL",   8},
  {"MSW",   6},
  {"SW",    7}
};

double compare_eucl(const list_num_t &mot1, const list_num_t &mot2,
    const bool norm = false) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  list_num_t diffmat(ncol, vec_num_t(nrow, 0.0));
  vec_bool_t good(ncol, false);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) good[i] = true;
  }

  for (std::size_t i = 0; i < ncol; ++i) {

    if (good[i]) {

      for (std::size_t j = 0; j < nrow; ++j) {
        diffmat[i][j] = mot1[i][j] - mot2[i][j];
        diffmat[i][j] = pow(diffmat[i][j], 2.0);
      }

    }

  }

  vec_num_t colsums(ncol, 0.0);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        colsums[i] += diffmat[i][j];
      }
      colsums[i] = sqrt(colsums[i]);
    }
  }

  double k = 0.0;
  int counter = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      ++counter;
      k += colsums[i];
    }
  }

  return norm ? k / double(counter) / sqrt(2.0) : k / sqrt(2.0);

}

double compare_pcc(const list_num_t &mot1, const list_num_t &mot2,
    const bool norm = false) {

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
      topbot[i] = bot[i] == 0 ? 0.0 : top[i] / bot[i];
      // topbot[i] = top[i] / bot[i];  // --> get Inf values if both mot1[i] and mot2[i]
    }                                  //     are uniform (e.g. {0.25, 0.25, 0.25, 0.25})
  }

  double total = std::accumulate(topbot.begin(), topbot.end(), 0.0);

  return norm ? total / double(n) : total;

}

double compare_kl(const list_num_t &mot1, const list_num_t &mot2,
    const bool norm = false) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  vec_bool_t good(ncol, false);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) good[i] = true;
  }

  list_num_t nmat(ncol, vec_num_t(nrow, 0.0));
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        nmat[i][j] = mot1[i][j] * log(mot1[i][j] / mot2[i][j]);
        nmat[i][j] += mot2[i][j] * log(mot2[i][j] / mot1[i][j]);
      }
    }
  }

  vec_num_t colsums(ncol, 0.0);
  int counter = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        colsums[i] += nmat[i][j];
      }
      ++counter;
    }
  }

  double total = 0.0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      total += colsums[i];
    }
  }

  return norm ? 0.5 / double(counter) * total : 0.5 * total;

}

double compare_sw(const list_num_t &mot1, const list_num_t &mot2,
    const bool norm = false) {

  std::size_t ncol = mot1.size(), nrow = mot1[0].size();
  vec_bool_t good(ncol, false);
  for (std::size_t i = 0; i < ncol; ++i) {
    if (mot1[i][0] >= 0 && mot2[i][0] >= 0) good[i] = true;
  }

  list_num_t nmat(ncol, vec_num_t(nrow, 0.0));
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        nmat[i][j] = pow(mot1[i][j] - mot2[i][j], 2.0);
      }
    }
  }

  vec_num_t colsums(ncol, 2.0);
  double total = 0.0;
  int counter = 0;
  for (std::size_t i = 0; i < ncol; ++i) {
    if (good[i]) {
      for (std::size_t j = 0; j < nrow; ++j) {
        colsums[i] -= nmat[i][j];
      }
      total += colsums[i];
      ++counter;
    }
  }

  return norm ? total / double(counter) : total;

}

void klfix(list_num_t &mot) {
  for (std::size_t i = 0; i < mot.size(); ++i) {
    for (std::size_t j = 0; j < mot[0].size(); ++j) {
      mot[i][j] += 0.01;
    }
  }
  return;
}

void equalize_mot_cols(list_num_t &mot1, list_num_t &mot2,
    vec_num_t &ic1, vec_num_t &ic2, const int overlap) {

  std::size_t nrow = mot1[0].size();
  std::size_t ncol1 = mot1.size();
  std::size_t ncol2 = mot2.size();
  std::size_t overlap1 = overlap, overlap2 = overlap;

  if (overlap < 1) {
    overlap1 *= ncol1;
    overlap2 *= ncol2;
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

  std::size_t ncol = mot.size(), nrow = mot[0].size();

  list_num_t rcmot(ncol, vec_num_t(nrow));
  for (std::size_t i = 0; i < nrow; ++i) {
    for (std::size_t j = 0; j < nrow; ++j) {
      rcmot[ncol - 1 - i][nrow - 1 - j] = mot[i][j];
    }
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

double compare_motif_pair(list_num_t mot1, list_num_t mot2,
    const std::string method, const int moverlap, const bool RC,
    vec_num_t ic1, vec_num_t ic2, const double minic, const bool norm,
    const double posic) {

  double ans_rc;
  if (RC) {
    list_num_t rcmot2 = get_motif_rc(mot2);
    vec_num_t rcic2(ic2.size());
    for (std::size_t i = 0; i < rcic2.size(); ++i) {
      rcic2[rcic2.size() - 1 - i] = ic2[i];
    }
    ans_rc = compare_motif_pair(mot1, rcmot2, method, moverlap, false,
        ic1, rcic2, minic, norm, posic);
  }

  switch (::metrics_enum[method]) {
    case 5:
    case 8: klfix(mot1);
            klfix(mot2);
  }

  std::size_t tlen = mot1.size() >= mot2.size() ? mot1.size() : mot2.size();

  equalize_mot_cols(mot1, mot2, ic1, ic2, moverlap);

  std::size_t minw = mot1.size() <= mot2.size() ? mot1.size() : mot2.size();
  list_num_t tmot1(minw), tmot2(minw);
  vec_num_t tic1(minw), tic2(minw);
  std::size_t fori = 1 + mot1.size() - minw, forj = 1 + mot2.size() - minw;
  vec_num_t ans;
  ans = RC ? vec_num_t(fori + forj) : vec_num_t(fori + forj - 1);
  int alignlen;
  double mic1, mic2;
  bool lowic = false;

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

      switch (::metrics_enum[method]) {

        case 1: ans[i + j] = lowic ? sqrt(2.0)
                                     : compare_eucl(tmot1, tmot2, false)
                                       * double(tlen) / double(alignlen);
                break;
        case 2: ans[i + j] = lowic ? sqrt(2.0) * double(alignlen)
                                     : compare_eucl(tmot1, tmot2, true)
                                       * double(tlen) / double(alignlen);
                break;
        case 3: ans[i + j] = lowic ? 0.0
                                     : compare_pcc(tmot1, tmot2, false)
                                       * double(alignlen) / double(tlen);
                break;
        case 4: ans[i + j] = lowic ? 0.0
                                     : compare_pcc(tmot1, tmot2, true)
                                       * double(alignlen) / double(tlen);
                break;
        case 5: ans[i + j] = lowic ? 5.0 * double(alignlen)
                                     : compare_kl(tmot1, tmot2, false)
                                       * double(tlen) / double(alignlen);
                break;
        case 6: ans[i + j] = lowic ? 0.0
                                     : compare_sw(tmot1, tmot2, true)
                                       * double(alignlen) / double(tlen);
                break;
        case 7: ans[i + j] = lowic ? 0.0
                                     : compare_sw(tmot1, tmot2, false)
                                       * double(alignlen) / double(tlen);
                break;
        case 8: ans[i + j] = lowic ? 5.0
                                     : compare_kl(tmot1, tmot2, true)
                                       * double(tlen) / double(alignlen);
                break;

      }

      lowic = false;

    }
  }

  if (RC) ans[ans.size() - 1] = ans_rc;

  switch (::metrics_enum[method]) {

    case 1: return *std::min_element(ans.begin(), ans.end());
    case 2: return *std::min_element(ans.begin(), ans.end());
    case 3: return *std::max_element(ans.begin(), ans.end());
    case 4: return *std::max_element(ans.begin(), ans.end());
    case 5: return *std::min_element(ans.begin(), ans.end());
    case 6: return *std::max_element(ans.begin(), ans.end());
    case 7: return *std::max_element(ans.begin(), ans.end());
    case 8: return *std::min_element(ans.begin(), ans.end());

  }

  return -1.0;

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

/* C++ ENTRY ---------------------------------------------------------------- */

// [[Rcpp::export(rng = false)]]
std::vector<double> compare_motifs_cpp(const Rcpp::List &mots,
    const std::vector<int> &index1, const std::vector<int> &index2,
    const std::string &method, int minoverlap, const bool RC,
    std::vector<std::vector<double>> &bkg, const int type, const bool relative,
    const double minic, const bool norm, const int nthreads, const double posic) {

  if (minoverlap < 1) minoverlap = 1;

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
  }

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    icscores[i].reserve(vmots[i].size());
    for (std::size_t j = 0; j < vmots[i].size(); ++j) {
      icscores[i].push_back(internal_posIC(vmots[i][j], bkg[i], type, relative));
    }
  }

  vec_num_t answers(index1.size());
  RcppThread::parallelFor(0, answers.size(),
      [&answers, &vmots, &index1, &index2, &icscores, &method, minoverlap, RC, minic, norm, posic]
      (std::size_t i) {
        answers[i] = compare_motif_pair(vmots[index1[i]], vmots[index2[i]], method,
            minoverlap, RC, icscores[index1[i]], icscores[index2[i]], minic, norm, posic);
      }, nthreads);

  return answers;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::vector<double>> compare_motifs_all_cpp(const Rcpp::List &mots,
    const std::string &method, int minoverlap, const bool RC,
    std::vector<std::vector<double>> &bkg, const int type, const bool relative,
    const double minic, const bool norm, const int nthreads, const double posic) {

  if (minoverlap < 1) minoverlap = 1;

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
  }

  for (std::size_t i = 0; i < vmots.size(); ++i) {
    icscores[i].reserve(vmots[i].size());
    for (std::size_t j = 0; j < vmots[i].size(); ++j) {
      icscores[i].push_back(internal_posIC(vmots[i][j], bkg[i], type, relative));
    }
  }

  list_num_t answers(vmots.size());
  RcppThread::parallelFor(0, answers.size() - 1,
      [&answers, &vmots, &icscores, &method, minoverlap, RC, minic, norm, posic]
      (std::size_t i) {
        answers[i].reserve(vmots.size() - i + 1);
        for (std::size_t j = i + 1; j < vmots.size(); ++j) {
          answers[i].push_back(compare_motif_pair(vmots[i], vmots[j], method,
              minoverlap, RC, icscores[i], icscores[j], minic, norm, posic));
        }
      }, nthreads);

  return answers;

}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_comparison_matrix(const std::vector<double> &ans,
    const std::vector<int> &index1, const std::vector<int> &index2,
    const std::string &method, const Rcpp::StringVector &motnames) {

  Rcpp::NumericMatrix out(motnames.size(), motnames.size());

  switch (::metrics_enum[method]) {
    case  4: out.fill_diag(1.0);
             break;
    case  6: out.fill_diag(2.0);
             break;
    default: out.fill_diag(0.0);
             break;
  }

  for (std::size_t i = 0; i < ans.size(); ++i) {
    out(index1[i], index2[i]) = ans[i];
    out(index2[i], index1[i]) = ans[i];
  }

  Rcpp::rownames(out) = motnames;
  Rcpp::colnames(out) = motnames;

  return out;

}
