#include <Rcpp.h>
#include <RcppThread.h>
#include <cmath>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_map>
#include <numeric>

// Timings on 20 char string, maxWordSize=7:
// WF:       3.16 us
// WF-fast:  2.17 us
// T:       12.14 us
// T-fast:  12.14 us
// DUST:     3.53 us

std::unordered_map<std::string, int> COMPLEXITY_METRICS = {
  {"WoottonFederhen", 1},
  {"WoottonFederhenFast", 2},
  {"Trifonov", 3},
  {"TrifonovFast", 4},
  {"DUST", 5}
};

enum COMPLEXITY_METRICS_ENUM {
  WoottonFederhen = 1,
  WoottonFederhenFast = 2,
  Trifonov = 3,
  TrifonovFast = 4,
  DUST = 5
};

std::vector<double> get_complexity_state_vector(const std::string &x, const std::string &alph) {
  std::map<char, double> S;
  for (std::size_t i = 0; i < alph.size(); ++i) {
    S[alph[i]] = 0.0;
  }
  for (std::size_t i = 0; i < x.size(); ++i) {
    S[x[i]] += 1.0;
  }
  std::vector<double> S_d(alph.size());
  for (std::size_t i = 0; i < S_d.size(); ++i) {
    S_d[i] = S[alph[i]];
  }
  return S_d;
}

// [[Rcpp::export(rng = false)]]
std::string get_alphabet_cpp(const std::string &x) {
  std::set<char> alph(x.begin(), x.end());
  return std::string(alph.begin(), alph.end());
}

double prod_cpp(const std::vector<double> &x) {
  // x = 12 is the max before overflowing ints; use doubles instead.
  double y = 1;
  for (int i = 0; i < int(x.size()); ++i) {
    y *= x[i];
  }
  return y;
}

double prod_cpp(const int x) {
  double y = 1;
  for (double i = 1; i <= double(x); ++i) {
    y *= i;
  }
  return y;
}

std::vector<double> count_unique_strings(const std::vector<std::string> &y) {
  std::set<std::string> y_unique_s(y.begin(), y.end());
  std::vector<std::string> y_unique(y_unique_s.begin(), y_unique_s.end());
  std::vector<double> counts(y_unique.size(), 0.0);
  for (std::size_t i = 0; i < y.size(); ++i) {
    for (std::size_t j = 0; j < y_unique.size(); ++j) {
      if (y[i] == y_unique[j]) {
        counts[j] += 1.0;
        break;
      }
    }
  }
  return counts;
}

// [[Rcpp::export]]
std::vector<std::vector<std::size_t>> calc_wins_cpp2(const std::size_t seqlen, const std::size_t window, const std::size_t overlap, const bool return_incomplete_window = false) {

  if (window > seqlen) return std::vector<std::vector<std::size_t>>();
  if (overlap >= window) return std::vector<std::vector<std::size_t>>();

  std::vector<std::size_t> starts;
  std::vector<std::size_t> stops;

  starts.push_back(1);
  stops.push_back(window);
  
  std::size_t i = 0;
  while (true) {
    if (starts[i] == seqlen) break;
    std::size_t next_start = starts[i] + window - overlap;
    std::size_t next_stop = next_start + window - 1;
    if (next_stop > seqlen) {
      next_stop = seqlen;
      if (next_start > seqlen) {
        next_start = seqlen;
      }
      if (next_stop == stops[i]) break;
      starts.push_back(next_start);
      stops.push_back(next_stop);
      break;
    }
    starts.push_back(next_start);
    stops.push_back(next_stop);
    ++i;
  }

  std::vector<std::vector<std::size_t>> wins(2);
  wins[0] = starts;
  wins[1] = stops;

  return wins;

}

std::vector<std::string> split_every_n_cpp(const std::string &x, const std::size_t n = 1) {
  std::vector<std::string> x_split(x.size() - n + 1);
  for (std::size_t i = 0; i < x_split.size(); ++i) {
    x_split[i] = x.substr(i, n);
  }
  return x_split;
}

// [[Rcpp::export]]
double dust_cpp(const std::string &x) {
  double l = x.size() - 2;
  std::vector<std::string> aSplit = split_every_n_cpp(x, 3);
  std::vector<double> aSplitCount = count_unique_strings(aSplit);
  std::vector<double> Sa(aSplitCount.size());
  for (std::size_t i = 0; i < Sa.size(); ++i) {
    Sa[i] = aSplitCount[i] * (aSplitCount[i] - 1.0) / 2.0;
  }
  return std::accumulate(Sa.begin(), Sa.end(), 0.0) / (l - 1.0);
}

// [[Rcpp::export]]
double trifonov_fast_cpp(const std::string &x, int maxWordSize, std::string alph = "") {
  // And not actually any faster...
  if (!alph.size()) alph = get_alphabet_cpp(x); 
  maxWordSize = std::min(maxWordSize, int(x.size()));
  std::size_t N = x.size(), K = alph.size();
  std::vector<double> V(maxWordSize);
  std::vector<double> V_max(maxWordSize);
  for (std::size_t i = 0; i < V.size(); ++i) {
    std::vector<std::string> x_split = split_every_n_cpp(x, i + 1);
    V[i] = double(std::set<std::string>(x_split.begin(), x_split.end()).size());
    V_max[i] = std::min(double(std::pow(K, i + 1)), double(N - i));
  }
  return std::accumulate(V.begin(), V.end(), 0.0) / std::accumulate(V_max.begin(), V_max.end(), 0.0);
}

// [[Rcpp::export]]
double trifonov_cpp(const std::string &x, int maxWordSize, std::string alph = "") {
  if (!alph.size()) alph = get_alphabet_cpp(x); 
  maxWordSize = std::min(maxWordSize, int(x.size()));
  std::size_t N = x.size(), K = alph.size();
  std::vector<double> CT(maxWordSize);
  for (std::size_t i = 0; i < CT.size(); ++i) {
    std::vector<std::string> x_split = split_every_n_cpp(x, i + 1);
    CT[i] = double(std::set<std::string>(x_split.begin(), x_split.end()).size());
    CT[i] /= std::min(double(std::pow(K, i + 1)), double(N - i));
  }
  return prod_cpp(CT);
}

// [[Rcpp::export]]
double wootton_federhen_fast_cpp(const std::string &x, std::string alph = "") {
  // Not that much faster tbh, but can be used with larger strings
  if (!alph.size()) alph = get_alphabet_cpp(x); 
  double answer = 0, L = x.size(), N = alph.size();
  const std::vector<double> S = get_complexity_state_vector(x, alph);
  for (std::size_t i = 0; i < S.size(); ++i) {
    if (!S[i]) continue;
    answer -= (S[i] / L) * (log(S[i] / L) / log(N));
  }
  return answer;
}

// [[Rcpp::export]]
double wootton_federhen_cpp(const std::string &x, std::string alph = "") {
  // Stops working for DNA strings length 171+
  if (!alph.size()) alph = get_alphabet_cpp(x); 
  double L = x.size(), N = alph.size();
  const std::vector<double> S = get_complexity_state_vector(x, alph);
  const double Pi_top = prod_cpp(x.size());
  std::vector<double> Pi_bot(S.size(), 1);
  for (std::size_t i = 0; i < S.size(); ++i) {
    if (!S[i]) continue;
    Pi_bot[i] = prod_cpp(S[i]);
  }
  double Pi = Pi_top / prod_cpp(Pi_bot);
  return (log(Pi) / log(N)) / L;
}

// [[Rcpp::export]]
std::vector<std::string> slide_windows_cpp(const std::string &x, const std::size_t window, const std::size_t overlap, const bool return_incomplete_window = false, const int nthreads = 1) {
  std::vector<std::vector<std::size_t>> win_locs = calc_wins_cpp2(x.size(), window, overlap, return_incomplete_window);
  std::vector<std::string> win_strs(win_locs[0].size());
  RcppThread::parallelFor(0, win_strs.size(),
      [&win_strs, &x, &win_locs] (std::size_t i) {
        win_strs[i] = x.substr(win_locs[0][i], win_locs[1][i] - win_locs[0][i] + 1);
      }, nthreads);
  return win_strs;
}

// [[Rcpp::export]]
std::vector<double> sliding_complexity_cpp(const std::string &x, const std::size_t window, const std::size_t overlap, const std::string metric, std::string alph = "", int maxWordSize = 7, const int nthreads = 1) {
  if (!alph.size()) alph = get_alphabet_cpp(x);
  std::vector<std::vector<std::size_t>> wins = calc_wins_cpp2(x.size(), window, overlap);
  if (!wins.size()) return std::vector<double>();
  std::vector<double> complexities(wins[0].size());
  switch (::COMPLEXITY_METRICS[metric]) {
    case WoottonFederhen: {
      RcppThread::parallelFor(0, complexities.size(),
          [&complexities, &x, &wins, &alph] (std::size_t i) {
            complexities[i] = wootton_federhen_cpp(x.substr(wins[0][i], wins[1][i] - wins[0][i] + 1), alph);
          }, nthreads);
      break;
    }
    case WoottonFederhenFast: {
      RcppThread::parallelFor(0, complexities.size(),
          [&complexities, &x, &wins, &alph] (std::size_t i) {
            complexities[i] = wootton_federhen_fast_cpp(x.substr(wins[0][i], wins[1][i] - wins[0][i] + 1), alph);
          }, nthreads);
      break;
    }
    case Trifonov: {
      if (maxWordSize > int(window)) maxWordSize = int(window);
      RcppThread::parallelFor(0, complexities.size(),
          [&complexities, &x, &wins, &alph, &maxWordSize] (std::size_t i) {
            complexities[i] = trifonov_cpp(x.substr(wins[0][i], wins[1][i] - wins[0][i] + 1), maxWordSize, alph);
          }, nthreads);
      break;
    }
    case TrifonovFast: {
      if (maxWordSize > int(window)) maxWordSize = int(window);
      RcppThread::parallelFor(0, complexities.size(),
          [&complexities, &x, &wins, &alph, &maxWordSize] (std::size_t i) {
            complexities[i] = trifonov_fast_cpp(x.substr(wins[0][i], wins[1][i] - wins[0][i] + 1), maxWordSize, alph);
          }, nthreads);
      break;
    }
    case DUST: {
      // Warning: DUST only works for DNA
      RcppThread::parallelFor(0, complexities.size(),
          [&complexities, &x, &wins] (std::size_t i) {
            complexities[i] = dust_cpp(x.substr(wins[0][i], wins[1][i] - wins[0][i] + 1));
          }, nthreads);
      break;
    }
    default: return std::vector<double>();
  }
  return complexities;
}
