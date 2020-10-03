#ifndef _SHUFFLE_SEQUENCES_
#define _SHUFFLE_SEQUENCES_

#include "types.h"

vec_str_t get_klet_strings(const vec_str_t &alph, const int &k);

vec_int_t klet_counter(const vec_int_t &single_seq, const int &k,
    const std::size_t &nlets, const std::size_t &alphlen);

#endif
