#include <Rcpp.h>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <limits>
#include "utils-internal.h"
#include "types.h"

std::unordered_map<std::string, int> TYPES_e = {
  {"PCM", 1},
  {"PPM", 2},
  {"PWM", 3},
  {"ICM", 4}
};

std::unordered_map<std::string, int> RDNA_e = {
  {"A",  1},
  {"C",  2},
  {"G",  3},
  {"T",  4},
  {"U",  5},
  {"R",  6},
  {"Y",  7},
  {"M",  8},
  {"K",  9},
  {"S", 10},
  {"W", 11},
  {"H", 12},
  {"B", 13},
  {"V", 14},
  {"D", 15},
  {"N", 16}
};

std::unordered_map<std::string, int> AA_e = {
  {"A",  1},
  {"C",  2},
  {"D",  3},
  {"E",  4},
  {"F",  5},
  {"G",  6},
  {"H",  7},
  {"I",  8},
  {"K",  9},
  {"L", 10},
  {"M", 11},
  {"N", 12},
  {"P", 13},
  {"Q", 14},
  {"R", 15},
  {"S", 16},
  {"T", 17},
  {"V", 18},
  {"W", 19},
  {"Y", 20}
};

const str_t FAIL1 = " * Incorrect type for '";
const str_t FAIL2 = " * Incorrect vector length for '";
const str_t FAIL3 = "': expected ";
const str_t FAIL4 = "; got ";

vec_str_t clean_up_check(const vec_str_t &fails) {

  std::size_t faillen = fails.size();
  vec_str_t out;
  out.reserve(faillen);

  for (std::size_t i = 0; i < faillen; ++i) {
    if (fails[i] != "") out.push_back(fails[i]);
  }

  return out;

}

Rcpp::NumericVector generate_pos(const std::vector<double> &bkg) {

  Rcpp::NumericVector rgam(bkg.size());
  for (std::size_t i = 0; i < bkg.size(); ++i) {
    rgam[i] = R::rgamma(bkg[i], 1.0);
  }

  double rgam_s = std::accumulate(rgam.begin(), rgam.end(), 0.0);

  for (std::size_t i = 0; i < bkg.size(); ++i) {
    rgam[i] /= rgam_s;
  }

  return rgam;

}

/* C++ ENTRY ---------------------------------------------------------------- */

// [[Rcpp::export(rng = false)]]
Rcpp::List split_gapped(const Rcpp::NumericMatrix &mot,
    const std::vector<int> &gaploc) {

  R_xlen_t nrow = mot.nrow(), ncol = mot.ncol();

  Rcpp::StringVector rnames = Rcpp::rownames(mot), cnames = Rcpp::colnames(mot);

  bool has_rnames = true, has_cnames = true;
  if (rnames[0] == R_NilValue) has_rnames = false;
  if (cnames[0] == R_NilValue) has_cnames = false;

  vec_int_t gaploc1;
  vec_int_t gaploc2;

  gaploc1.reserve(gaploc.size() + 1);
  gaploc2.reserve(gaploc.size() + 1);

  gaploc1.push_back(0);

  for (std::size_t i = 0; i < gaploc.size(); ++i) {
    gaploc1.push_back(gaploc[i]);
    gaploc2.push_back(gaploc[i] - 1);
  }

  gaploc2.push_back(ncol - 1);

  vec_int_t sizes(gaploc.size() + 1);
  for (std::size_t i = 0; i < gaploc.size() + 1; ++i) {
    sizes[i] = gaploc2[i] - gaploc1[i] + 1;
  }

  Rcpp::List out(gaploc.size() + 1);

  int counter;

  for (std::size_t i = 0; i < gaploc.size() + 1; ++i) {

    Rcpp::NumericMatrix tmp(nrow, sizes[i]);
    counter = 0;
    for (int j = gaploc1[i]; j <= gaploc2[i]; ++j) {
      tmp(Rcpp::_, counter) = mot(Rcpp::_, j);
      ++counter;
    }

    if (has_rnames) {
      Rcpp::rownames(tmp) = rnames;
    }

    if (has_cnames) {
      Rcpp::StringVector tmp_cnames = cnames[Rcpp::seq(gaploc1[i], gaploc2[i])];
      Rcpp::colnames(tmp) = tmp_cnames;
    }

    out[i] = tmp;

  }

  return out;

}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix generate_motif(const int ncol, const std::vector<double> &bkg) {

  Rcpp::NumericMatrix out(bkg.size(), ncol);
  for (int i = 0; i < ncol; ++i) {
    out(Rcpp::_, i) = generate_pos(bkg);
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
Rcpp::List min_max_doubles() {

  return Rcpp::List::create(
        Rcpp::_["min"] = -std::numeric_limits<double>::max(),
        Rcpp::_["max"] = std::numeric_limits<double>::max()
      );

}

// [[Rcpp::export(rng = false)]]
std::vector<std::vector<int>> comb2_cpp(const int n) {

  int outlen = pow(n, 2) / 2 + n / 2 + 1;
  list_int_t out(2);
  out[0].reserve(outlen);
  out[1].reserve(outlen);

  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      out[0].push_back(i);
      out[1].push_back(j);
    }
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector table_cpp(const Rcpp::StringVector &x) {
  return Rcpp::table(x);
}

// [[Rcpp::export(rng = false)]]
Rcpp::StringVector sort_unique_cpp(const Rcpp::StringVector &x) {
  return Rcpp::sort_unique(x);
}

// [[Rcpp::export(rng = false)]]
Rcpp::StringVector collapse_rows_mat(const Rcpp::CharacterMatrix &seqs_k) {

  Rcpp::StringVector out(seqs_k.nrow());

  for (R_xlen_t i = 0; i < seqs_k.nrow(); ++i) {
    out[i] = Rcpp::collapse(seqs_k(i, Rcpp::_));
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
Rcpp::StringVector collapse_cols_mat(const Rcpp::CharacterMatrix &seqs_k) {

  Rcpp::StringVector out(seqs_k.ncol());

  for (R_xlen_t i = 0; i < seqs_k.ncol(); ++i) {
    out[i] = Rcpp::collapse(seqs_k(Rcpp::_, i));
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
Rcpp::StringVector collapse_rows_df(const Rcpp::DataFrame &seqs_k) {

  Rcpp::CharacterMatrix seqs_k_mat(seqs_k.nrow(), seqs_k.size());

  for (R_xlen_t i = 0; i < seqs_k.size(); ++i)
    seqs_k_mat(Rcpp::_, i) = Rcpp::StringVector(seqs_k[i]);

  return collapse_rows_mat(seqs_k_mat);

}

// [[Rcpp::export(rng = false)]]
Rcpp::String collapse_cpp(const Rcpp::StringVector &x) {
  return Rcpp::collapse(x);
}

// [[Rcpp::export(rng = false)]]
void print_pb(const R_xlen_t &out) {
  if (out >= 10 && out < 100) {
    Rprintf("\b\b\b\b %i%%", out);
    return;
  }
  if (out > 0 && out < 10) {
    Rprintf("\b\b\b\b  %i%%", out);
    return;
  }
  switch (out) {
    case   0: Rprintf("   0%%");            return;
    case 100: Rprintf("\b\b\b\b%i%%", out); return;
    case  -1: Rprintf("\b\b\b\b100%%\n");   return;
  }
  Rcpp::stop("Input must be an integer in between -1 and 100");
}

// [[Rcpp::export(rng = false)]]
void update_pb(const R_xlen_t &i, const R_xlen_t &max) {

  R_xlen_t out, prev = i - 1;
  if (i == max)
    out = -1;
  else
    out = 100 * i / max;

  if (prev > 0 && out != -1) {
    prev = 100 * prev / max;
    if (prev != out) print_pb(out);
  } else if (out == -1) {
    print_pb(out);
  }

}

// [[Rcpp::export(rng = false)]]
Rcpp::String all_checks_collapse(const Rcpp::StringVector &checks) {

  R_xlen_t n = checks.size();

  Rcpp::StringVector out_pre(n * 2);
  R_xlen_t i_ = 0;
  for (R_xlen_t i = 0; i < n * 2; ++i) {
    if (i % 2 == 0) {
      out_pre[i] = "\n";
    } else {
      out_pre[i] = checks[i_];
      i_ += 1;
    }
  }

  return Rcpp::collapse(out_pre);

}

// [[Rcpp::export(rng = false)]]
std::vector<double> pcm_to_ppmC(std::vector<double> pos,
    const double pseudocount = 0.0) {

  double spos = std::accumulate(pos.begin(), pos.end(), 0.0);
  double num = pos.size();

  if (pseudocount > 0) {
    for (std::size_t i = 0; i < pos.size(); ++i) {
      pos[i] = (pos[i] + pseudocount / num) / (spos + pseudocount);
    }
  } else {
    for (std::size_t i = 0; i < pos.size(); ++i) {
      pos[i] = pos[i] / spos;
    }
  }

  return pos;

}

// [[Rcpp::export(rng = false)]]
std::vector<double> ppm_to_pcmC(std::vector<double> pos, double nsites = 0) {

  if (std::isnan(nsites) || nsites <= 1) nsites = 100;

  std::size_t n = pos.size();

  for (std::size_t i = 0; i < n; ++i) {
    pos[i] = round(pos[i] * nsites);
  }

  double spos = std::accumulate(pos.begin(), pos.end(), 0.0);

  if (spos != nsites) {
    double fix = nsites - spos;
    std::size_t tochange = std::distance(pos.begin(),
        std::max_element(pos.begin(), pos.end()));
    pos[tochange] += fix;
  }

  return pos;

}

// [[Rcpp::export(rng = false)]]
std::vector<double> ppm_to_pwmC(std::vector<double> pos,
    std::vector<double> bkg, const double pseudocount = 0.0,
    double nsites = 100) {

  std::size_t n_pos = pos.size(), n_bkg = bkg.size();

  if (std::isnan(nsites) || nsites <= 1) nsites = 100;
  if (n_bkg != n_pos) {
    bkg = vec_num_t(n_pos, 1 / double(n_pos));
  }

  if (pseudocount > 0) {
    pos = ppm_to_pcmC(pos, nsites);
    pos = pcm_to_ppmC(pos, pseudocount);
  }

  for (std::size_t i = 0; i < n_pos; ++i) {
    pos[i] = log2(pos[i] / bkg[i]);
  }

  return pos;

}

// [[Rcpp::export(rng = false)]]
std::vector<double> pwm_to_ppmC(std::vector<double> pos,
    std::vector<double> bkg) {

  std::size_t n_pos = pos.size(), n_bkg = bkg.size();

  if (n_bkg != n_pos) {
    bkg = vec_num_t(n_pos, 1 / double(n_pos));
  }

  for (std::size_t i = 0; i < n_pos; ++i) {
    pos[i] = pow(2.0, pos[i]);
  }

  double spos = std::accumulate(pos.begin(), pos.end(), 0.0);

  if (spos > 0.99 && spos < 1.01) return pos;

  for (std::size_t i = 0; i < n_pos; ++i) {
    pos[i] *= bkg[i];
  }

  spos = std::accumulate(pos.begin(), pos.end(), 0.0);

  if (spos > 0.99 && spos < 1.01) return pos;

  for (std::size_t i = 0; i < n_pos; ++i) {
    pos[i] /= spos;
  }

  return pos;

}

// [[Rcpp::export(rng = false)]]
std::vector<double> ppm_to_icmC(std::vector<double> pos,
    std::vector<double> bkg, const bool relative_entropy = false) {

  std::size_t n_pos = pos.size(), n_bkg = bkg.size();

  if (n_bkg != n_pos) {
    bkg = vec_num_t(n_pos, 1 / double(n_pos));
  } else {
    double sbkg = std::accumulate(bkg.begin(), bkg.end(), 0.0);
    if (sbkg > 1.01 || sbkg < 0.99) {
      for (std::size_t i = 0; i < n_pos; ++i) {
        bkg[i] /= sbkg;
      }
    }
  }

  if (relative_entropy) {
    for (std::size_t i = 0; i < n_pos; ++i) {
      double tmp = pos[i] / bkg[i];
      pos[i] *= tmp >= 0 ? log2(tmp) : 0.0;
      if (pos[i] < 0) pos[i] = 0.0;
    }
  } else {
    vec_num_t heights(n_pos);
    for (std::size_t i = 0; i < n_pos; ++i) {
      heights[i] = pos[i] > 0 ? -pos[i] * log2(pos[i]) : 0.0;
    }
    double height_after = std::accumulate(heights.begin(), heights.end(), 0.0);
    double total_ic = log2(double(n_pos)) - height_after;
    for (std::size_t i = 0; i < n_pos; ++i) {
      pos[i] *= total_ic;
    }
  }

  return pos;

}

// [[Rcpp::export(rng = false)]]
double position_icscoreC(std::vector<double> pos,
    std::vector<double> bkg, const std::string &type = "PPM",
    const double pseudocount = 1.0, double nsites = 100.0,
    const bool relative_entropy = false) {

  if (nsites < 1) nsites = 100;

  std::size_t n_pos = pos.size(), n_bkg = bkg.size();

  if (n_bkg != n_pos) {
    bkg = vec_num_t(n_pos, 1 / double(n_pos));
  }

  switch (::TYPES_e[type]) {
    case 1: pos = pcm_to_ppmC(pos, pseudocount);
            break;
    case 2: pos = ppm_to_pcmC(pos, nsites);
            pos = pcm_to_ppmC(pos, pseudocount);
            break;
    case 3: pos = pwm_to_ppmC(pos, bkg);
            break;
    case 4: return std::accumulate(pos.begin(), pos.end(), 0.0);
    default: Rcpp::stop("incorrect type");
  }

  pos = ppm_to_icmC(pos, bkg, relative_entropy);

  return std::accumulate(pos.begin(), pos.end(), 0.0);

}

// [[Rcpp::export(rng = false)]]
std::vector<double> icm_to_ppmC(std::vector<double> pos) {

  double total_ic = std::accumulate(pos.begin(), pos.end(), 0.0);
  for (std::size_t i = 0; i < pos.size(); ++i) {
    pos[i] /= total_ic;
  }

  return pos;

}

// [[Rcpp::export(rng = false)]]
std::string get_consensusC(std::vector<double> pos,
    const std::string &alphabet = "DNA", const std::string &type = "PPM",
    const double pseudocount = 1.0) {

  switch (::TYPES_e[type]) {
    case 1: pos = pcm_to_ppmC(pos, pseudocount); break;
    case 3: pos = pwm_to_ppmC(pos, vec_num_t()); break;
    case 4: pos = icm_to_ppmC(pos); break;
  }

  // single letter consensus

  if (pos[0] > 0.5 && pos[1] <= 0.25 && pos[2] <= 0.25 && pos[3] <= 0.25) return "A";
  if (pos[1] > 0.5 && pos[0] <= 0.25 && pos[2] <= 0.25 && pos[3] <= 0.25) return "C";
  if (pos[2] > 0.5 && pos[0] <= 0.25 && pos[1] <= 0.25 && pos[3] <= 0.25) return "G";
  if (pos[3] > 0.5 && pos[0] <= 0.25 && pos[1] <= 0.25 && pos[2] <= 0.25) {
    if (alphabet == "DNA") return "T"; else return "U";
  }

  // two letter consensus

  if ((pos[0] + pos[1]) > 0.75) return "M";
  if ((pos[0] + pos[2]) > 0.75) return "R";
  if ((pos[0] + pos[3]) > 0.75) return "W";
  if ((pos[1] + pos[2]) > 0.75) return "S";
  if ((pos[1] + pos[3]) > 0.75) return "Y";
  if ((pos[2] + pos[3]) > 0.75) return "K";

  // three letter consensus

  if (pos[0] > 0.25 && pos[1] > 0.25 && pos[3] > 0.25) return "H";
  if (pos[1] > 0.25 && pos[2] > 0.25 && pos[3] > 0.25) return "B";
  if (pos[0] > 0.25 && pos[1] > 0.25 && pos[2] > 0.25) return "V";
  if (pos[0] > 0.25 && pos[2] > 0.25 && pos[3] > 0.25) return "D";

  // no consensus

  return "N";

}

// [[Rcpp::export(rng = false)]]
std::vector<double> consensus_to_ppmC(const std::string &letter) {

  switch (::RDNA_e[letter]) {

    case  1: return { 0.997, 0.001, 0.001, 0.001 };  // A
    case  2: return { 0.001, 0.997, 0.001, 0.001 };  // C
    case  3: return { 0.001, 0.001, 0.997, 0.001 };  // G
    case  4: return { 0.001, 0.001, 0.001, 0.997 };  // T
    case  5: return { 0.001, 0.001, 0.001, 0.997 };  // U
    case  6: return { 0.499, 0.001, 0.499, 0.001 };  // R
    case  7: return { 0.001, 0.499, 0.001, 0.449 };  // Y
    case  8: return { 0.499, 0.499, 0.001, 0.001 };  // M
    case  9: return { 0.001, 0.001, 0.499, 0.499 };  // K
    case 10: return { 0.001, 0.499, 0.499, 0.001 };  // S
    case 11: return { 0.499, 0.001, 0.001, 0.499 };  // W
    case 12: return { 0.333, 0.333, 0.001, 0.333 };  // H
    case 13: return { 0.001, 0.333, 0.333, 0.333 };  // B
    case 14: return { 0.333, 0.333, 0.333, 0.001 };  // V
    case 15: return { 0.333, 0.001, 0.333, 0.333 };  // D
    case 16: return {  0.25,  0.25,  0.25,  0.25 };  // N
    default: return {  0.25,  0.25,  0.25,  0.25 };

  }

}

// [[Rcpp::export(rng = false)]]
std::vector<double> consensus_to_ppmAAC(const std::string &letter) {

  int let = ::AA_e[letter];
  vec_num_t out(20, 0.001);

  if (let == 0) {

    if (letter == "B") {
      out[2] = 0.491;
      out[11] = 0.491;
    } else if (letter == "Z") {
      out[3] = 0.491;
      out[13] = 0.491;
    } else if (letter == "J") {
      out[7] = 0.491;
      out[9] = 0.491;
    } else {
      return vec_num_t(20, 0.05);
    }

  } else {

    out[let - 1] = 0.981;

  }

  return out;

}

// [[Rcpp::export(rng = false)]]
std::string get_consensusAAC(std::vector<double> pos,
    const std::string &type = "PPM", const double pseudocount = 0.0) {

  switch (::TYPES_e[type]) {
    case 1: pos = pcm_to_ppmC(pos, pseudocount); break;
    case 3: pos = pwm_to_ppmC(pos, vec_num_t()); break;
    case 4: pos = icm_to_ppmC(pos); break;
  }

  if      (pos[2] >= 0.4 && pos[11] >= 0.4) return "B";
  else if (pos[3] >= 0.4 && pos[13] >= 0.4) return "Z";
  else if (pos[7] >= 0.4 && pos[9]  >= 0.4) return "J";

  if (*std::max_element(pos.begin(), pos.end()) < 0.1) return "X";

  vec_num_t pos2 = pos;
  std::sort(pos2.begin(), pos2.end());
  if (pos2[19] == pos2[18]) return "X";

  return ::AMINOACIDS2[std::distance(pos.begin(), std::max_element(pos.begin(), pos.end()))];

}

// [[Rcpp::export(rng = false)]]
std::vector<std::string> check_fun_params(const Rcpp::List &param_args,
    std::vector<int> param_len, std::vector<bool> param_null,
    int expected_type) {

  // param_args: A named list of arguments.
  //
  // param_len:  A vector of integers, each entry corresponding to an arg. The
  //             values represent length of the argument. For arguments which
  //             can be of any length, 0 can be used. The default is length 1
  //             for all args, which is set as param_len = numeric(). Note that
  //             length checks are not performed on S4 objects.
  //
  // param_null: A vector of booleans, each entry entry corresponding to an
  //             arg. For args which allow NULL or missing values, TRUE is used.
  //             The default is not allowing NULL or missing for all args,
  //             which is set as param_null = logical().
  //
  // expected_type: Expected type integer. One of 16 (character),
  //             14 (numeric), 10 (logical), 25 (S4).

  R_xlen_t arglen = param_args.size();

  if (param_len.size() > 1 && int(param_len.size()) != int(arglen)) {
    Rcpp::stop("incorrect param_len");
  }
  if (param_null.size() > 1 && int(param_null.size()) != int(arglen)) {
    Rcpp::stop("incorrect param_null");
  }

  vec_str_t fails(arglen * 2, "");
  vec_str_t parnames = param_args.names();

  if (param_len.size() == 0)
    param_len = vec_int_t(arglen, 1);
  else if (param_len.size() == 1 && arglen > 1)
    param_len = vec_int_t(arglen, param_len[0]);

  if (param_null.size() == 0)
    param_null = vec_bool_t(arglen, false);
  else if (param_null.size() == 1 && arglen > 1)
    param_null = vec_bool_t(arglen, param_null[0]);

  str_t arg_name, exp_type, obs_type, exp_len_c, obs_len_c;
  bool null_check, argfail;
  for (R_xlen_t i = 0; i < arglen; ++i) {

    Rcpp::RObject arg = param_args[i];
    int argtype = arg.sexp_type();
    if (argtype == 13) argtype = 14;
    null_check = param_null[i];
    argfail = false;
    arg_name = parnames[i];

    switch (expected_type) {
      case 10: exp_type = "`logical`"; break;
      case 14: exp_type = "`numeric`"; break;
      case 16: exp_type = "`character`"; break;
      case 25: exp_type = "`S4`"; break;
      default: Rcpp::stop("incorrect expected_type");
    }

    switch (argtype) {
      case  0: obs_type = "`NULL`"; break;
      case  1: obs_type = "`missing`"; break;
      case 10: obs_type = "`logical`"; break;
      case 14: obs_type = "`numeric`"; break;
      case 16: obs_type = "`character`"; break;
      case 25: obs_type = "`S4`"; break;
      default: obs_type = "`unknown`"; break;
    }

    if (null_check) {
      if (argtype != 0 && argtype != 1 && argtype != expected_type) {
        fails[i] += ::FAIL1 + arg_name + ::FAIL3 + exp_type + ::FAIL4 + obs_type;
        argfail = true;
      }
    } else if (argtype != expected_type) {
      fails[i] += ::FAIL1 + arg_name + ::FAIL3 + exp_type + ::FAIL4 + obs_type;
      argfail = true;
    }

    if (!argfail && argtype != 0 && !null_check && argtype != 25) {
      R_xlen_t sarglen;
      switch (argtype) {
        case 10: {
                   Rcpp::LogicalVector arg_ = Rcpp::as<Rcpp::LogicalVector>(arg);
                   sarglen = arg_.size();
                   break;
                 }
        case 14: {
                   Rcpp::NumericVector arg_ = Rcpp::as<Rcpp::NumericVector>(arg);
                   sarglen = arg_.size();
                   break;
                 }
        case 16: {
                   Rcpp::StringVector arg_ = Rcpp::as<Rcpp::StringVector>(arg);
                   sarglen = arg_.size();
                   break;
                 }
        default: Rcpp::stop("unrecognised parameter object type");
      }
      int explen = param_len[i];
      exp_len_c = std::to_string(explen);
      obs_len_c = std::to_string(sarglen);
      if (sarglen != explen && explen != 0) {
        fails[((i + 1) * 2) - 1] += ::FAIL2 + arg_name + ::FAIL3 + exp_len_c
          + ::FAIL4 + obs_len_c;
      }
    }

  }

  return clean_up_check(fails);

}
