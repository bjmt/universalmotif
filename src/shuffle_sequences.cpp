#include <Rcpp.h>
#include <RcppThread.h>
#include <cmath>
#include <random>
#include <set>
#include "types.h"
#include "utils-internal.h"

vec_int_t klet_counter(const vec_int_t &single_seq, const int &k,
    const std::size_t &nlets, const std::size_t &alphlen) {

  vec_int_t klet_counts(nlets, 0);
  int l, counter;

  for (std::size_t i = 0; i < single_seq.size() - k + 1; ++i) {
    l = 0; counter = 0;
    for (int j = k - 1; j >= 0; --j) {
      l += pow(alphlen, j) * single_seq[i + counter];
      ++counter;
    }
    ++klet_counts[l];
  }

  return klet_counts;

}

int get_lastlet(const vec_int_t &single_seq, const int &k,
    const std::size_t &alphlen) {

  int lastlet = 0;
  for (int i = k - 2; i >= 0; --i) {
    lastlet += pow(alphlen, i) * single_seq[single_seq.size() - 1 - i];
  }

  return lastlet;

}

vec_int_t get_firstlet(const vec_int_t &single_seq, const int &k) {

  vec_int_t firstlet;
  firstlet.reserve(k - 1);

  for (int i = 0; i < k - 1; ++i) {
    firstlet.push_back(single_seq[i]);
  }

  return firstlet;

}

list_int_t get_edgecounts(const vec_int_t &klet_counts, const std::size_t &mlets,
    const std::size_t &alphlen) {

  list_int_t edgecounts(mlets, vec_int_t(alphlen));
  std::size_t counter = 0;

  for (std::size_t i = 0; i < mlets; ++i) {
    for (std::size_t j = 0; j < alphlen; ++j) {
      edgecounts[i][j] = klet_counts[counter];
      ++counter;
    }
  }

  return edgecounts;

}

vec_bool_t get_emptyvertices(const std::size_t &mlets, const std::size_t &alphlen,
    const list_int_t &edgelist) {

  vec_bool_t emptyvertices;
  emptyvertices.reserve(mlets);

  for (std::size_t i = 0; i < mlets; ++i) {
    emptyvertices.push_back(true);
    for (std::size_t j = 0; j < alphlen; ++j) {
      if (edgelist[i][j] > 0) {
        emptyvertices[i] = false;
        break;
      }
    }
  }

  return emptyvertices;

}

vec_int_t get_eulerpath(const list_int_t &edgelist, const int &lastlet,
    const std::size_t &mlets, const std::size_t &alphlen, const int &k,
    const vec_bool_t &emptyvertices, std::mt19937 gen) {

  vec_int_t eulerpath(mlets, 0);
  vec_bool_t vertices(mlets, false);
  vec_int_t next_index(mlets);
  int counter = 0, vertex_ok = 0;
  std::size_t mmlets = pow(alphlen, k - 2);
  int u;

  vertices[lastlet] = true;

  for (std::size_t i = 0; i < mlets; ++i) {
    next_index[i] = counter * alphlen;
    if (counter == int(mmlets) - 1)
      counter = 0;
    else
      ++counter;
  }

  for (std::size_t i = 0; i < mlets; ++i) {
    if (emptyvertices[i])
      vertices[i] = true;
    else
      ++vertex_ok;
  }

  for (std::size_t i = 0; i < mlets; ++i) {

    u = i;

    while (!vertices[u]) {
      std::discrete_distribution<int> next_u(edgelist[u].begin(), edgelist[u].end());
      eulerpath[u] = next_u(gen);
      u = k == 2 ? eulerpath[u] : eulerpath[u] + next_index[u];
    }

    u = i;

    while (!vertices[u]) {
      vertices[u] = true;
      u = k == 2 ? eulerpath[u] : eulerpath[u] + next_index[u];
    }

  }

  return eulerpath;

}

list_int_t get_edgelist(const list_int_t &edgecounts, const vec_int_t &eulerpath,
    const std::size_t &mlets, const std::size_t &alphlen, const int &lastlet,
    std::mt19937 gen, const vec_bool_t &emptyvertices) {

  list_int_t edgelist(mlets);
  int b;

  for (std::size_t i = 0; i < mlets; ++i) {
    if (emptyvertices[i]) continue;
    edgelist[i].reserve(std::accumulate(edgecounts[i].begin(), edgecounts[i].end(), 0));
    for (int j = 0; j < int(alphlen); ++j) {
      b = edgecounts[i][j];
      for (int h = 0; h < b; ++h) {
        edgelist[i].push_back(j);
      }
    }
    std::shuffle(edgelist[i].begin(), edgelist[i].end(), gen);
    if (int(i) != lastlet) edgelist[i].push_back(eulerpath[i]);
  }

  return edgelist;

}

vec_int_t eulerian_walk(const list_int_t &edgelist, const std::size_t &seqlen,
    const vec_int_t &firstlet, const std::size_t &mlets, const std::size_t &alphlen) {

  vec_int_t shuffled_seq_ints;
  shuffled_seq_ints.reserve(seqlen);
  vec_int_t edge_index(mlets, 0);
  int which_vertex;

  for (std::size_t i = 0; i < firstlet.size(); ++i) {
    shuffled_seq_ints.push_back(firstlet[i]);
  }

  for (std::size_t i = firstlet.size() - 1; i < seqlen - 1; ++i) {
    which_vertex = 0;
    for (int j = int(firstlet.size()) - 1; j >= 0; --j) {
      which_vertex += pow(alphlen, j) * (shuffled_seq_ints[i - j] % alphlen);
    }
    shuffled_seq_ints.push_back(edgelist[which_vertex][edge_index[which_vertex]]);
    ++edge_index[which_vertex];
  }

  return shuffled_seq_ints;

}

std::string make_new_seq(const vec_int_t &shuffled_seq_ints,
    const std::string &alph) {

  std::string out;
  out.reserve(shuffled_seq_ints.size());

  for (std::size_t i = 0; i < shuffled_seq_ints.size(); ++i) {
    out += alph[shuffled_seq_ints[i]];
  }

  return out;

}

std::string shuffle_euler_one(const std::string &single_seq, const int &k,
    std::mt19937 gen) {

  std::set<int> alph_s;
  for (std::size_t i = 0; i < single_seq.size(); ++i) {
    alph_s.insert(single_seq[i]);
  }
  std::string alph;
  alph.assign(alph_s.begin(), alph_s.end());
  std::size_t alphlen = alph.size();
  std::size_t nlets = pow(alph.size(), k);
  std::size_t mlets = pow(alph.size(), k - 1);

  vec_int_t seq_ints(single_seq.size());
  for (std::size_t i = 0; i < single_seq.size(); ++i) {
    for (std::size_t a = 0; a < alph.size(); ++a) {
      if (single_seq[i] == alph[a]) {
        seq_ints[i] = a;
        break;
      }
    }
  }

  vec_int_t klet_counts = klet_counter(seq_ints, k, nlets, alphlen);
  vec_int_t firslet = get_firstlet(seq_ints, k);
  int lastlet = get_lastlet(seq_ints, k, alphlen);
  list_int_t edgecounts = get_edgecounts(klet_counts, mlets, alphlen);
  vec_bool_t emptyvertices = get_emptyvertices(mlets, alphlen, edgecounts);
  vec_int_t eulerpath = get_eulerpath(edgecounts, lastlet, mlets, alphlen, k,
      emptyvertices, gen);

  for (std::size_t i = 0; i < eulerpath.size(); ++i) {
    if (int(i) != lastlet) --edgecounts[i][eulerpath[i]];
  }

  list_int_t edgelist = get_edgelist(edgecounts, eulerpath, mlets, alphlen,
      lastlet, gen, emptyvertices);

  vec_int_t shuffled_seq_ints = eulerian_walk(edgelist, single_seq.size(),
      firslet, mlets, alphlen);

  std::string out = make_new_seq(shuffled_seq_ints, alph);
  return out;

}

list_int_t make_klet_lists(const std::size_t &nlets, const int &k,
    const std::size_t &alphlen) {

  list_int_t klet_list(nlets);
  for (std::size_t i = 0; i < nlets; ++i) {
    klet_list[i].reserve(k);
  }
  int counter, let, step;

  for (int i = k - 1; i >= 0; --i) {
    counter = 0;
    let = 0;
    step = pow(alphlen, i);
    while (counter < int(nlets)) {
      for (int j = 0; j < step; ++j) {
        klet_list[counter].push_back(let);
        ++counter;
      }
      if (let == int(alphlen) - 1)
        let = 0;
      else
        ++let;
    }
  }

  return klet_list;

}

/*
vec_int_t get_next_let_list(const list_int_t &nlet_lists, const int &k,
    const std::size_t &alphlen, const std::size_t &nlets) {

  vec_int_t out(nlets, 0);

  for (std::size_t i = 0; i < nlets; ++i) {
    for (int j = 0; j < k; ++j) {
      out[i] += nlet_lists[i][j] * pow(alphlen, k - j - 1);
    }
  }

  return out;

}
*/

vec_int_t markov_generator(const std::size_t &seqsize, const vec_int_t &nlet_counts,
    const list_int_t &transitions, std::mt19937 gen,
    const std::size_t &nlets, const int &k, const std::size_t &alphlen) {

  vec_int_t out;
  out.reserve(seqsize);

  list_int_t nlet_lists = make_klet_lists(nlets, k, alphlen);

  std::discrete_distribution<int> first_gen(nlet_counts.begin(), nlet_counts.end());
  int firstletters = first_gen(gen);
  for (int i = 0; i < k; ++i) {
    out.push_back(nlet_lists[firstletters][i]);
  }

  int previous_mlet;
  for (std::size_t i = k - 1; i < seqsize; ++i) {
    previous_mlet = 0;
    for (int j = k - 1; j >= 1; --j) {
      previous_mlet += out[i - j] * pow(alphlen, j - 1);
    }
    std::discrete_distribution<int> next_gen(transitions[previous_mlet].begin(),
        transitions[previous_mlet].end());
    out.push_back(next_gen(gen));
  }

  return out;

}

std::string shuffle_markov_one(const std::string &single_seq, const int &k,
    std::mt19937 gen) {

  std::set<int> alph_s;
  for (std::size_t i = 0; i < single_seq.size(); ++i) {
    alph_s.insert(single_seq[i]);
  }
  std::string alph;
  alph.assign(alph_s.begin(), alph_s.end());
  std::size_t alphlen = alph.size();
  std::size_t nlets = pow(alphlen, k);
  std::size_t mlets = pow(alphlen, k - 1);

  vec_int_t seq_ints(single_seq.size());
  for (std::size_t i = 0; i < single_seq.size(); ++i) {
    for (std::size_t a = 0; a < alphlen; ++a) {
      if (single_seq[i] == alph[a]) {
        seq_ints[i] = a;
        break;
      }
    }
  }

  vec_int_t nlet_counts = klet_counter(seq_ints, k, nlets, alphlen);
  list_int_t transitions = get_edgecounts(nlet_counts, mlets, alphlen);

  vec_int_t out_ints = markov_generator(seq_ints.size(), nlet_counts, transitions,
      gen, nlets, k, alphlen);

  std::string out = make_new_seq(out_ints, alph);

  return out;

}

std::string shuffle_linear_one (const std::string &single_seq, const int &k,
    std::mt19937 gen) {

  std::size_t seqlen = single_seq.size();
  std::size_t seqlen_k = seqlen / k;
  std::size_t remaind = seqlen % k;
  std::size_t remaindlen = seqlen - remaind;

  std::string out;
  out.reserve(seqlen);

  vec_int_t indices;
  indices.reserve(seqlen_k);
  for (std::size_t i = 0; i < seqlen_k; ++i) {
    indices.push_back(i * k);
  }

  shuffle(indices.begin(), indices.end(), gen);

  for (std::size_t i = 0; i < seqlen_k; ++i) {
    for (int j = 0; j < k; ++j) {
      out += single_seq[indices[i] + j];
    }
  }

  if (remaind > 0) {
    for (std::size_t i = remaindlen; i < seqlen; ++i) {
      out += single_seq[i];
    }
  }

  return out;

}

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

vec_int_t klet_counter_from_string(const str_t &single_seq, const int &k) {

  std::set<int> alph_s;
  for (std::size_t i = 0; i < single_seq.size(); ++i) {
    alph_s.insert(single_seq[i]);
  }
  str_t alph;
  alph.assign(alph_s.begin(), alph_s.end());
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

vec_str_t get_klet_strings(const vec_str_t &alph, const int &k) {

  int alphlen = alph.size();
  int nlets = pow(alphlen, k);
  vec_str_t out(nlets, "");
  int let, counter, step;

  for (int i = k - 1; i >= 0; --i) {
    counter = 0; let = 0;
    step = pow(alphlen, i);
    while (counter < nlets) {
      for (int j = 0; j < step; ++j) {
        out[counter] += alph[let];
        ++counter;
      }
      if (let == alphlen - 1)
        let = 0;
      else
        ++let;
    }
  }

  return out;

}

/* C++ ENTRY ---------------------------------------------------------------- */

// [[Rcpp::export(rng = false)]]
std::vector<std::string> shuffle_markov_cpp(const std::vector<std::string> &sequences,
    const int &k, const int &nthreads, const int &seed) {

  unsigned int useed = seed;

  vec_str_t out(sequences.size());
  RcppThread::parallelFor(0, sequences.size(),
      [&out, &sequences, &useed, &k] (std::size_t i) {

        std::mt19937 gen(useed * (int(i) + 1));
        out[i] = shuffle_markov_one(sequences[i], k, gen);

      }, nthreads);

  return out;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> shuffle_euler_cpp(const std::vector<std::string> &sequences,
    const int &k, const int &nthreads, const int &seed) {

  unsigned int useed = seed;

  vec_str_t out(sequences.size());
  RcppThread::parallelFor(0, sequences.size(),
      [&out, &sequences, &k, &useed] (std::size_t i) {

        std::mt19937 gen(useed * (int(i) + 1));
        out[i] = shuffle_euler_one(sequences[i], k, gen);

      }, nthreads);

  return out;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> shuffle_linear_cpp(const std::vector<std::string> &sequences,
    const int &k, const int &nthreads, const int &seed) {

  unsigned int useed = seed;

  vec_str_t out(sequences.size());
  RcppThread::parallelFor(0, sequences.size(),
      [&out, &sequences, &k, &useed] (std::size_t i) {

        std::mt19937 gen(useed * (int(i) + 1));
        out[i] = shuffle_linear_one(sequences[i], k, gen);

      }, nthreads);

  return out;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> shuffle_k1_cpp(const std::vector<std::string> &sequences,
    const int &nthreads, const int &seed) {

  unsigned int useed = seed;

  vec_str_t out(sequences.size());
  RcppThread::parallelFor(0, sequences.size(),
      [&out, &sequences, &useed] (std::size_t i) {
        std::mt19937 gen(useed * (int(i) + 1));
        out[i] = sequences[i];
        shuffle(out[i].begin(), out[i].end(), gen);
      }, nthreads);

  return out;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::vector<int>> count_klets_cpp(const std::vector<std::string> &sequences,
    const int &k, const int &nthreads) {

  list_int_t counts(sequences.size());
  RcppThread::parallelFor(0, sequences.size(),
      [&counts, &sequences, &k] (std::size_t i) {
        counts[i] = klet_counter_from_string(sequences[i], k);
      }, nthreads);

  return counts;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> split_seq_by_win(std::string &seq1,
    const std::vector<int> &start, const std::vector<int> &stop) {

  std::vector<std::string> out(start.size());

  for (R_xlen_t i = 0; i < start.size(); ++i) {
    out[i] = seq1.substr(start[i] - 1, stop[i] - start[i] + 1);
  }

  return out;

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

// [[Rcpp::export(rng = false)]]
std::vector<std::string> get_klets_cpp(std::vector<std::string> &alph,
    const int &k) {

  vec_str_t klets = get_klet_strings(alph, k);

  return klets;

}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> create_sequences_cpp(const int seqlen,
    const int seqnum, const std::vector<std::string> &alph, const int k,
    const std::vector<double> &freqs, const int nthreads, const int seed,
    const Rcpp::NumericMatrix &transitions) {

  unsigned int useed = seed;

  std::size_t alphlen = alph.size();
  std::size_t nlets = pow(alphlen, k);
  list_int_t nlet_lists = make_klet_lists(nlets, k, alphlen);

  list_num_t trans = R_to_cpp_motif_num(transitions);

  vec_str_t out(seqnum, "");

  if (k == 1) {

    RcppThread::parallelFor(0, out.size(),
        [&seqlen, &alph, &useed, &out, &freqs] (std::size_t i) {

          out[i].reserve(seqlen);
          std::mt19937 gen(useed * (int(i) + 1));
          std::discrete_distribution<int> nextlet(freqs.begin(), freqs.end());

          for (int j = 0; j < seqlen; ++j) {
            out[i] += alph[nextlet(gen)];
          }

        }, nthreads);

  } else if (k > 1) {

    RcppThread::parallelFor(0, out.size(),
        [&seqlen, &alph, &useed, &out, &freqs, &trans, &nlet_lists, &k, &alphlen]
        (std::size_t i) {

          vec_int_t out_i;
          out_i.reserve(seqlen);
          std::mt19937 gen(useed * (int(i) + 1));

          std::discrete_distribution<int> let1(freqs.begin(), freqs.end());
          int firstletters = let1(gen);
          for (int j = 0; j < k; ++j) {
            out_i.push_back(nlet_lists[firstletters][j]);
          }

          int prevmlet;
          for (int j = k - 1; j < seqlen; ++j) {
            prevmlet = 0;
            for (int b = k - 1; b >= 1; --b) {
              prevmlet += out_i[j - b] * pow(alphlen, b - 1);
            }
            std::discrete_distribution<int> nextlet(trans[prevmlet].begin(),
                trans[prevmlet].end());
            out_i.push_back(nextlet(gen));
          }

          out[i].reserve(seqlen);
          for (int j = 0; j < seqlen; ++j) {
            out[i] += alph[out_i[j]];
          }

        }, nthreads);

  }

  return out;

}
