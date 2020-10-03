#include <Rcpp.h>
#include <RcppThread.h>
#include <cmath>
#include "types.h"
#include "shuffle_sequences.h"

vec_int_t klet_counter_with_alph(const str_t &single_seq, const str_t &alph,
    const int &k) {

  std::size_t alphlen = alph.size();
  std::size_t nlets = pow(alphlen, k);

  vec_int_t seq_ints(single_seq.size());
  for (std::size_t i = 0; i < single_seq.size(); ++i) {
    for (std::size_t a = 0; a < alphlen; ++a) {
      if (single_seq[i] == alph[a]) {
        seq_ints[i] = a;
        break;
      }
    }
  }

  vec_int_t counts = klet_counter(seq_ints, k, nlets, alphlen);

  return counts;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::vector<int>> count_klets_alph_cpp(const std::vector<std::string> &sequences,
    const std::string &alph, const int &k, const int &nthreads) {

  list_int_t counts(sequences.size());
  RcppThread::parallelFor(0, sequences.size(),
      [&counts, &sequences, &k, &alph] (std::size_t i) {
        counts[i] = klet_counter_with_alph(sequences[i], alph, k);
      }, nthreads);

  return counts;

}
