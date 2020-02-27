#include <Rcpp.h>
#include "types.h"

int peakfinder_single_cpp(int i, const Rcpp::NumericVector &x, int m) {

  int z = i - m;
  if (z < 0) z = 0;

  int w = i + m;
  if (w > x.size() - 1) w = x.size() - 1;

  Rcpp::IntegerVector x1 = Rcpp::seq(z, i);
  Rcpp::IntegerVector x2 = Rcpp::seq(i, w);

  Rcpp::NumericVector x_1 = x[x1];
  Rcpp::NumericVector x_2 = x[x2];

  if (Rcpp::is_true(Rcpp::all(x_1 <= x[i])) && Rcpp::is_true(Rcpp::all(x_2 <= x[i])))
    return i + 1;
  else
    return NA_INTEGER;

}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector peakfinder_cpp(const Rcpp::NumericVector &x, int m = 3) {

  Rcpp::IntegerVector shape = Rcpp::diff(Rcpp::sign(Rcpp::diff(x)));

  Rcpp::IntegerVector shape_count = Rcpp::seq(0, shape.size() - 1);
  Rcpp::IntegerVector shape_which = shape_count[shape < 0];

  Rcpp::IntegerVector pks(shape_which.size());

  for (R_xlen_t i = 0; i < shape_which.size(); ++i) {
    pks[i] = peakfinder_single_cpp(shape_which[i] + 1, x, m);
  }

  return pks[!is_na(pks)];

}

// [[Rcpp::export(rng = false)]]
std::vector<double> linbin_cpp(const std::vector<int> &x, const std::vector<int> &gpoints) {

  double M = gpoints.size(), b = gpoints.size(), n = x.size(), a = 1;

  std::vector<double> gcnts(M, 0.0);
  double delta = (b - a) / (M - 1);

  double lxi, rem, li;

  for (int i = 0; i < n; ++i) {
    lxi = ((x[i] - a) / delta) + 1;
    li = trunc(lxi);
    rem = lxi - li;
    if (li > 1 && li < M) {
      gcnts[li - 1] += (1 - rem);
      gcnts[li] += rem;
    }
  }

  return gcnts;

}
