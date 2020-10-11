#include <Rcpp.h>
#include <RcppThread.h>
#include <cmath>
#include "types.h"

vec_int_t klet_counter_NA(const vec_int_t &single_seq, const int &k,
    const std::size_t &nlets, const std::size_t &alphlen) {

  vec_int_t klet_counts(nlets, 0);
  int l, counter;

  for (std::size_t i = 0; i < single_seq.size() - k + 1; ++i) {
    l = 0; counter = 0;
    for (int j = k - 1; j >= 0; --j) {
      l += pow(alphlen, j) * single_seq[i + counter];
      ++counter;
    }

    if (l >= 0) ++klet_counts[l];
    // Un-comment this line if you decide to use remainder count in probability
    // calculations.
    // else ++klet_counts[nlets];

  }

  return klet_counts;

}

vec_int_t klet_counter_with_alph(const str_t &single_seq, const str_t &alph,
    const int &k) {

  std::size_t alphlen = alph.size();
  std::size_t nlets = pow(alphlen, k);

  vec_int_t seq_ints(single_seq.size(), -1000000);
  for (std::size_t i = 0; i < single_seq.size(); ++i) {
    for (std::size_t a = 0; a < alphlen; ++a) {
      if (single_seq[i] == alph[a]) {
        seq_ints[i] = a;
        break;
      }
    }
  }

  vec_int_t counts = klet_counter_NA(seq_ints, k, nlets, alphlen);

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
