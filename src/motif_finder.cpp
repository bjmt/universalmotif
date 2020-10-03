#include <Rcpp.h>
#include <RcppThread.h>
#include "types.h"
#include "utils-internal.h"
#include "shuffle_sequences.h"

vec_int_t seq_string_to_int(const str_t &seq1, const str_t &alph,
    const std::size_t &alphlen) {

  vec_int_t out(seq1.size());
  for (std::size_t i = 0; i < seq1.size(); ++i) {
    for (std::size_t a = 0; a < alphlen; ++a) {
      if (seq1[i] == alph[a]) {
        out[i] = a;
        break;
      }
    }
  }

  return out;

}

double calc_seq_prob(const str_t &seq1, const vec_num_t &bkg,
    const str_t &alph, const std::size_t &alphlen) {

  vec_int_t seq1i = seq_string_to_int(seq1, alph, alphlen);
  double out = 1;
  for (std::size_t i = 0; i < seq1.size(); ++i) {
    out *= bkg[seq1i[i]];
  }
  return out;
}

//------------------------------------------------------------------------------

// // [[Rcpp::export(rng = false)]]
// Rcpp::List motif_finder_no_bkg_cpp(const std::vector<std::string> &sequences,
//     const std::string &alph, const int &nmotifs, const double &pvalue,
//     const int &minsize, const int &maxsize, const int &nthreads) {
//
//
//
// }

// [[Rcpp::export(rng = false)]]
std::vector<double> calc_seq_probs_cpp(const std::vector<std::string> &seqs,
    const std::vector<double> &bkg, const std::string &alph,
    const int &nthreads) {

  std::size_t alphlen = alph.size();

  vec_num_t probs(seqs.size());
  RcppThread::parallelFor(0, seqs.size(),
      [&probs, &seqs, &bkg, &alph, &alphlen] (std::size_t i) {
        probs[i] = calc_seq_prob(seqs[i], bkg, alph, alphlen);
      }, nthreads);

  return probs;

}
