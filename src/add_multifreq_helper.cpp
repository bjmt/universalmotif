#include <Rcpp.h>
using namespace Rcpp;

StringVector paste_cpp(StringVector string_in) {

  int n = string_in.length();
  std::string out = as<std::string>(string_in(0));

  for (int i = 1; i < n; ++i) {
    out += as<std::string>(string_in[i]);
  }

  return out;

}

// [[Rcpp::export]]
StringVector single_to_k(StringVector seq, int k) {

  int n = seq.length();
  StringMatrix seq_mat(k, n - k + 1);
  StringVector out(n - k + 1);
  StringVector out_(k);
  StringVector out_i;

  for (int i = 0; i < seq_mat.ncol(); ++i) {
    seq_mat(0, i) = seq[i];
  }

  for (int i = 1; i < k; ++i) {
    for (int j = 0; j < seq_mat.ncol(); ++j) {
      seq_mat(i, j) = seq(j + i);
    }
  }

  for (int i = 0; i < out.length(); ++i) {
    for (int j = 0; j < k; ++j) {
      out_[j] = seq_mat(j, i);
    }
    out_i = paste_cpp(out_);
    out[i] = out_i[0];
  }

  return out;

}
