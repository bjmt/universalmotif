#include <Rcpp.h>
using namespace Rcpp;

int peakfinder_single_cpp(int i, NumericVector x, int m) {

  int z = i - m;
  if (z < 0) z = 0;

  int w = i + m;
  if (w > x.length() - 1) w = x.length() - 1;

  IntegerVector x1 = seq(z, i);
  IntegerVector x2 = seq(i, w);

  NumericVector x_1 = x[x1];
  NumericVector x_2 = x[x2];

  if (is_true(all(x_1 <= x[i])) && is_true(all(x_2 <= x[i])))
    return i + 1;
  else
    return NA_INTEGER;

}

// [[Rcpp::export(rng = false)]]
IntegerVector peakfinder_cpp(NumericVector x, int m = 3) {

  IntegerVector shape = diff(sign(diff(x)));

  IntegerVector shape_count = seq(0, shape.length() - 1);
  IntegerVector shape_which = shape_count[shape < 0];

  IntegerVector pks(shape_which.length());

  for (int i = 0; i < shape_which.length(); ++i) {
    pks[i] = peakfinder_single_cpp(shape_which[i] + 1, x, m);
  }

  return pks[!is_na(pks)];

}

// [[Rcpp::export(rng = false)]]
IntegerVector linbin_cpp(IntegerVector x, IntegerVector gpoints) {

  int M = gpoints.length();
  IntegerVector gcnts(M);

  double delta = (gpoints[M - 1] - gpoints[0]) / (M - 1);

  for (int i = 0; i < x.length(); ++i) {

    double lxi = (x[i] - gpoints[0]) / delta;
    int li = lxi;
    double rem = lxi - li;

    if (li >= 0 && li < M) {
      gcnts[li] += 1 - rem;
      gcnts[li + 1] += rem;
    }

  }

  return gcnts;

}
