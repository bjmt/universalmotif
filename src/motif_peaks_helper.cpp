#include <Rcpp.h>

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
Rcpp::IntegerVector linbin_cpp(const Rcpp::IntegerVector &x,
    const Rcpp::IntegerVector &gpoints) {

  R_xlen_t M = gpoints.size();
  Rcpp::IntegerVector gcnts(M);

  double delta = (gpoints[M - 1] - gpoints[0]) / (M - 1);

  for (R_xlen_t i = 0; i < x.size(); ++i) {

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
