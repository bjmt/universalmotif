#include <Rcpp.h>
#include <RcppThread.h>
#include <cmath>
#include "types.h"

vec_int_t klet_counter_NA(const vec_int_t &single_seq, const int &k,
    const std::size_t &nlets, const std::size_t &alphlen) {

  vec_int_t klet_counts(nlets, 0);
  if (single_seq.size() < std::size_t(k)) return klet_counts;

  // Precompute alphlen^j for j = 0..k-1 once instead of calling pow() at every
  // sequence position.
  std::vector<long long> pow_table(k);
  pow_table[0] = 1;
  for (int j = 1; j < k; ++j) pow_table[j] = pow_table[j - 1] * (long long)(alphlen);

  // k-mers containing a letter outside the alphabet (marked as negative) are
  // skipped.
  for (std::size_t i = 0; i < single_seq.size() - k + 1; ++i) {
    long long l = 0;
    bool known = true;
    for (int j = 0; j < k; ++j) {
      int let = single_seq[i + j];
      if (let < 0) {
        known = false;
        break;
      }
      l += pow_table[k - 1 - j] * let;
    }
    if (known) ++klet_counts[l];
  }

  return klet_counts;

}

vec_int_t klet_counter_with_alph(const str_t &single_seq, const str_t &alph,
    const int &k) {

  std::size_t alphlen = alph.size();
  std::size_t nlets = pow(alphlen, k);

  // For a string with 25 million characters, about 150 ms can be saved by using
  // char indexing on an int array. However, I'm not sure how to deal with
  // non-ASCII chars in that case. Maybe don't bother. Either way for k < 5 the
  // output is returned in less than a second.
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
