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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
IntegerVector linbin_cpp(IntegerVector x, IntegerVector gpoints) {

  int M = gpoints.length();
  IntegerVector gcnts(M);

  double delta = (gpoints[M - 1] - gpoints[0]) / (M - 1);

  for (int i = 0; i < x.length(); ++i) {

    double lxi = ((x[i] - gpoints[0]) / delta);
    int li = lxi;
    double rem = lxi - li;

    if (li >= 0 && li < M) {
      gcnts[li] += 1 - rem;
      gcnts[li + 1] += rem;
    }

  }

  return gcnts;

}

// // [[Rcpp::export]]
// List kern_cpp(IntegerVector x, int bandwidth, int gridsize, int range_x) {
//
  // Rcpp::Environment base("package:stats");
  // Rcpp::Function fft_r = base["fft"];
//
  // int n = x.length();
  // int M = gridsize;
  // int tau = 4;
  // int h = bandwidth;
  // int a = 1;
  // int b = range_x;
//
  // IntegerVector gpoints = seq(a, b);
  // IntegerVector gcounts = linbin_cpp(x, gpoints);
//
  // double delta = (b - a) / (h * (M - 1));
  // int L = min(IntegerVector::create(floor(tau / delta), M));
//
  // IntegerVector lvec = seq(0, L);
  // NumericVector kappa = dnorm(lvec * delta) / (n * h);
//
  // NumericVector P_2 = ceiling(NumericVector::create(log(M + L + 1) / log(2)));
  // double P_2_ = P_2[0];
  // double P = pow(2, P_2_);
//
  // NumericVector kappa_last = rev(kappa[seq(1, kappa.length() - 1)]);
  // int kappa_mid_len = P - 2 * L - 1;
  // NumericVector kappa_new(kappa.length() + kappa_mid_len + kappa_last.length());
  // for (int i = 0; i < kappa.length(); ++i) {
    // kappa_new[i] = kappa[i];
  // }
  // for (int i = 0; i < kappa_last.length(); ++i) {
    // kappa_new[i + kappa.length() + kappa_mid_len] = kappa_last[i];
  // }
//
  // double tot = sum(kappa_new) * (b - a) / (M - 1) * n;
//
  // kappa = fft_r(kappa_new / tot, false);
  // gcounts = fft_r(gcounts, false);
//
  // NumericVector kappa_gcounts(kappa_new.length());
  // for (int i = 0; i < kappa_new.length(); ++i) {
    // kappa_gcounts[i] = kappa_new[i] * gcounts[i];
  // }
//
// // Re() not working here
  // NumericVector gcounts_y = Re(fft_r(kappa_gcounts, true)) / P;
//
  // List out = List::create(_["x"] = gpoints, _["y"] = gcounts_y[seq(1, M)]);
  // return out;
//
// }
