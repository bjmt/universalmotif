#include <Rcpp.h>
#include <sstream>  // needed for `to_string` work around
using namespace Rcpp;

StringVector dna = StringVector::create(

    "A", "C", "G", "T"

);

StringVector rna = StringVector::create(

    "A", "C", "G", "U"

);

StringVector types = StringVector::create(

    "PCM", "PPM", "PWM", "ICM"

);

StringVector rdna = StringVector::create(

    "A", "C", "G", "T",
    "U", "R", "Y", "M",
    "K", "S", "W", "H",
    "B", "V", "D", "N"

);

StringVector aa = StringVector::create(

    "A", "C", "D", "E",
    "F", "G", "H", "I",
    "K", "L", "M", "N",
    "P", "Q", "R", "S",
    "T", "V", "W", "Y"

);

namespace std {
  template<typename T>
  std::string to_string(const T &n) {

  // This is required to get around `to_string` g++ compiler error

  // https://stackoverflow.com/questions/19122574/to-string-isnt-a-member-of-std/19122592
  // https://stackoverflow.com/questions/12975341/to-string-is-not-a-member-of-std-says-g-mingw

    std::ostringstream s;
    s << n;
    return s.str();
  }
}

// [[Rcpp::export(rng = false)]]
IntegerVector table_cpp(StringVector x) {
  return table(x);
}

// [[Rcpp::export(rng = false)]]
StringVector collapse_rows_mat(CharacterMatrix seqs_k) {

  // ~10 times as fast as apply(seqs_k_matrix, 1, collapse_cpp)
  // ~15 faster than collapse_rows_df(seqs_k_df)

  StringVector out(seqs_k.nrow());

  for (int i = 0; i < seqs_k.nrow(); ++i) {
    out[i] = collapse(seqs_k(i, _));
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
StringVector collapse_cols_mat(CharacterMatrix seqs_k) {

  StringVector out(seqs_k.ncol());

  for (int i = 0; i < seqs_k.ncol(); ++i) {
    out[i] = collapse(seqs_k(_, i));
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
StringVector collapse_rows_df(DataFrame seqs_k) {

  // >2 times as fast as apply(seqs_k_df, 1, collapse_cpp)

  CharacterMatrix seqs_k_mat(seqs_k.nrow(), seqs_k.size());

  for (int i = 0; i < seqs_k.size(); ++i)
    seqs_k_mat(_, i) = StringVector(seqs_k[i]);

  return collapse_rows_mat(seqs_k_mat);

}

//' @rdname utilities
//' @export
// [[Rcpp::export(rng = false)]]
StringVector get_klets(StringVector lets, int k) {

  // ~5 times faster and ~1/2 the memory allocations versus:
  //
  // sort(collapse_rows_df(expand.grid(rep(list(lets), k),
  //                                   stringsAsFactors = FALSE)))
  //
  // ~11 times faster and ~1/3 the memory allocations versus:
  //
  // sort(apply(expand.grid(rep(list(lets),k), stringsAsFactors= FALSE),
  //            1, paste0, collapse = ""))
  //
  // Slightly faster than collpase_rows_mat(RcppAlgos::permuteGeneral) with low
  // k, and slightly slower with high k (similar timings for lets=DNA_BASES and
  // k=4)

  if (k <= 0) stop("`k` must be greater than 0");

  lets = sort_unique(lets);

  int n1 = lets.length();
  int n2 = pow(n1, k);
  StringMatrix out(n2, k + 1);

  for (int i = 0; i < k ; ++i) {
    out(_, i) = rep(rep_each(lets, pow(n1, k - i - 1)), pow(n1, i + 1));
  }

  return collapse_rows_mat(out);

}

// [[Rcpp::export(rng = false)]]
String collapse_cpp(StringVector x) {

  // collapse_cpp(x) is about 3 times faster than base::paste(x, collapse = "")
  // collapse_cpp(c(x, y)) about 2 times faster than base::paste0(x, y)

  return collapse(x);

}

// [[Rcpp::export(rng = false)]]
void print_pb(int out) {
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
  stop("Input must be an integer in between -1 and 100");
}

// [[Rcpp::export(rng = false)]]
void update_pb(int i, int max) {

  int out, prev = i - 1;
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

StringVector strsplit_cpp(std::string x) {  // slightly slower than
  int n = x.size();                         // strsplit(x, "")[[1]]
  StringVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = x.substr(i, 1);
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
String all_checks_collapse(StringVector checks) {

  int n = checks.length();

  StringVector out_pre(n * 2);
  int i_ = 0;
  for (int i = 0; i < n * 2; ++i) {
    if (i % 2 == 0) {
      out_pre[i] = "\n";
    } else {
      out_pre[i] = checks[i_];
      i_ += 1;
    }
  }

  return collapse(out_pre);

}

// [[Rcpp::export(rng = false)]]
NumericVector pcm_to_ppmC(NumericVector position, double pseudocount=0) {

  double possum = sum(position);
  int num = position.size();
  NumericVector out(num);

  if (pseudocount != 0) {

    out = (position + (pseudocount / num)) / (possum + pseudocount);

  } else {

    out = position / possum;

  }

  return out;

}

// [[Rcpp::export(rng = false)]]
NumericVector ppm_to_pcmC(NumericVector position, double nsites=0) {

  if (nsites == 0) nsites = 100;
  int n = position.size();
  for (int i = 0; i < n; ++i) {
    position[i] = round(position[i] * nsites);
  }

  double possum = sum(position);
  if (possum != nsites) {
    double fix = nsites - possum;
    int tochange = which_max(position);
    position[tochange] += fix;
  }

  return position;

}

// [[Rcpp::export(rng = false)]]
NumericVector ppm_to_pwmC(NumericVector position, NumericVector bkg=0,
    double pseudocount=0, NumericVector nsites=NumericVector::create()) {

  int n_pos = position.size();
  int n_bkg = bkg.size();

  if (nsites.length() == 0) nsites = 100;
  else if (nsites[0] == 0) nsites = 100;

  double n_pos2 = n_pos;
  NumericVector bkg2(n_pos, 1.0 / n_pos2);
  if (n_pos != n_bkg) bkg = bkg2;

  if (pseudocount != 0) {
    position = ppm_to_pcmC(position, nsites[0]);
    position = pcm_to_ppmC(position, pseudocount);
  }

  for (int i = 0; i < n_pos; ++i) {
    position[i] = log2(position[i] / bkg[i]);
  }

  return position;

}

// [[Rcpp::export(rng = false)]]
NumericVector pwm_to_ppmC(NumericVector position, NumericVector bkg=0) {

  int n_pos = position.size();
  int n_bkg = bkg.size();
  double n_pos2 = n_pos;

  NumericVector bkg2(n_pos, 1.0 / n_pos2);
  if (n_pos != n_bkg) bkg = bkg2;

  for (int i = 0; i < n_pos; ++i) {
    position[i] = pow(2.0, position[i]);
  }

  double possum = sum(position);

  if (possum > 0.99 && possum < 1.01) return position;

  position = position * bkg;
  possum = sum(position);

  if (possum > 0.99 && possum < 1.01) return position;

  position = position / possum;

  return position;

}

// [[Rcpp::export(rng = false)]]
NumericVector ppm_to_icmC(NumericVector position, NumericVector bkg=0,
    bool relative_entropy=false) {

  int n_pos = position.size();
  int n_bkg = bkg.size();

  double n_pos2 = n_pos;
  if (n_pos != n_bkg) bkg = rep(1.0 / n_pos2, n_pos);
  else {
    double s_bkg = sum(bkg);
    if (s_bkg > 1.01 || s_bkg < 0.99)
      bkg = bkg / s_bkg;
  }

  if (relative_entropy) {
    for (int i = 0; i < n_pos; ++i) {
      position[i] = position[i] * log2(position[i] / bkg[i]);
      if (NumericVector::is_na(position[i])) position[i] = 0;
      if (position[i] < 0) position[i] = 0;
    }
    return position;
  } else {
    NumericVector heights(n_pos);
    for (int i = 0; i < n_pos; ++i) {
      heights[i] = -position[i] * log2(position[i]);
      if (NumericVector::is_na(heights[i])) heights[i] = 0;
    }
    double height_after = sum(heights);
    double total_ic = log2(n_pos2) - height_after;
    position = position * total_ic;
    return position;
  }

}

// [[Rcpp::export(rng = false)]]
double position_icscoreC(NumericVector position, NumericVector bkg=0,
    String type="PPM", double pseudocount=1, double nsites=100,
    bool relative_entropy=false) {

  if (nsites == 1) nsites = 100;

  int n_pos = position.size();
  int n_bkg = bkg.size();

  double n_pos2 = n_pos;
  NumericVector bkg2(n_pos, 1.0 / n_pos2);
  if (n_pos != n_bkg) bkg = bkg2;

  int ntype = match(StringVector::create(type), ::types)[0];
  switch (ntype) {
    case  1: position = pcm_to_ppmC(position, pseudocount);
             break;
    case  2: position = ppm_to_pcmC(position, nsites);
             position = pcm_to_ppmC(position, pseudocount);
             break;
    case  3: position = pwm_to_ppmC(position, bkg);
             break;
    case  4: return sum(position);
    default: stop("Incorrect type");
  }

  if (relative_entropy) {
    for (int i = 0; i < n_pos; ++i) {
      position[i] *= log2(position[i] / bkg[i]);
      if (NumericVector::is_na(position[i])) position[i] = 0.0;
      if (position[i] < 0) position[i] = 0.0;
    }
    return sum(position);
  } else {
    NumericVector heights(n_pos);
    for (int i = 0; i < n_pos; ++i) {
      heights[i] = -position[i] * log2(position[i]);
      if (NumericVector::is_na(heights[i])) heights[i] = 0.0;
    }
    double height_after = sum(heights);
    double total_ic = log2(n_pos2) - height_after;
    return total_ic;
  }

}

// [[Rcpp::export(rng = false)]]
NumericVector icm_to_ppmC(NumericVector position) {

  double total_ic = sum(position);
  position = position / total_ic;

  return position;

}

// [[Rcpp::export(rng = false)]]
String get_consensusC(NumericVector pos, String alphabet="DNA",
    String type="PPM", double pseudocount=1) {

  int ntype = match(StringVector::create(type), ::types)[0];
  switch (ntype) {
    case 1: pos = pcm_to_ppmC(pos, pseudocount); break;
    case 2: break;
    case 3: pos = pwm_to_ppmC(pos); break;
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
NumericVector consensus_to_ppmC(String letter) {

  int let = match(StringVector::create(letter), ::rdna)[0];
  switch (let) {

    case  1: return NumericVector::create(0.997, 0.001, 0.001, 0.001);  // A
    case  2: return NumericVector::create(0.001, 0.997, 0.001, 0.001);  // C
    case  3: return NumericVector::create(0.001, 0.001, 0.997, 0.001);  // G
    case  4: return NumericVector::create(0.001, 0.001, 0.001, 0.997);  // T
    case  5: return NumericVector::create(0.001, 0.001, 0.001, 0.997);  // U
    case  6: return NumericVector::create(0.499, 0.001, 0.499, 0.001);  // R
    case  7: return NumericVector::create(0.001, 0.499, 0.001, 0.449);  // Y
    case  8: return NumericVector::create(0.499, 0.499, 0.001, 0.001);  // M
    case  9: return NumericVector::create(0.001, 0.001, 0.499, 0.499);  // K
    case 10: return NumericVector::create(0.001, 0.499, 0.499, 0.001);  // S
    case 11: return NumericVector::create(0.499, 0.001, 0.001, 0.499);  // W
    case 12: return NumericVector::create(0.333, 0.333, 0.001, 0.333);  // H
    case 13: return NumericVector::create(0.001, 0.333, 0.333, 0.333);  // B
    case 14: return NumericVector::create(0.333, 0.333, 0.333, 0.001);  // V
    case 15: return NumericVector::create(0.333, 0.001, 0.333, 0.333);  // D
    case 16: return NumericVector::create( 0.25,  0.25,  0.25,  0.25);  // N
    default: return NumericVector::create( 0.25,  0.25,  0.25,  0.25);

  }

}

// [[Rcpp::export(rng = false)]]
NumericVector consensus_to_ppmAAC(String letter) {

  int let = match(StringVector::create(letter), ::aa)[0];
  NumericVector out(20, 0.001);

  if (IntegerVector::is_na(let)) {

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
      return rep(0.05, 20);
    }

  } else {

    out[let - 1] = 0.981;

  }

  return out;

}

// [[Rcpp::export(rng = false)]]
String get_consensusAAC(NumericVector position, String type="PPM",
    double pseudocount=0.0) {

  int ntype = match(StringVector::create(type), ::types)[0];
  switch (ntype) {
    case 1: position = pcm_to_ppmC(position, pseudocount);
            break;
    case 2: break;
    case 3: position = pwm_to_ppmC(position);
            break;
    case 4: position = icm_to_ppmC(position);
            break;
  }

  if      (position[2] >= 0.4 && position[11] >= 0.4) return "B";
  else if (position[3] >= 0.4 && position[13] >= 0.4) return "Z";
  else if (position[7] >= 0.4 && position[9]  >= 0.4) return "J";

  if (is_true(all(position < 0.1))) return "X";

  NumericVector position2 = Rcpp::clone(position);
  std::sort(position2.begin(), position2.end());
  if (position2[19] == position2[18]) return "X";

  return ::aa[which_max(position)];

}

// [[Rcpp::export(rng = false)]]
StringVector clean_up_check(StringVector fails) {

  int fails_len = fails.length();
  LogicalVector fails_keep(fails_len, true);

  for (int i = 0; i < fails_len; ++i) {
    if (fails[i] == "") fails_keep[i] = false;
  }

  return fails[fails_keep];

}

// [[Rcpp::export(rng = false)]]
StringVector check_fun_params(List param_args, IntegerVector param_len,
    LogicalVector param_null, int expected_type) {

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

  int arg_len = param_args.length();
  StringVector fails(arg_len * 2);
  StringVector param_names = param_args.names();

  IntegerVector param_len2;
  LogicalVector param_null2;

  if (param_len.length() == 0) {
    param_len2 = rep(1, arg_len);
  } else if (param_len.length() == 1 && arg_len > 1) {
    param_len2 = rep(param_len[0], arg_len);
  } else {
    param_len2 = param_len;
  }

  if (param_null.length() == 0) {
    param_null2 = rep(false, arg_len);
  } else if (param_null.length() == 1 && arg_len > 1) {
    param_null2 = rep(param_null[0], arg_len);
  } else {
    param_null2 = param_null;
  }

  StringVector fail_i;
  String arg_name, exp_type, obs_type, exp_len_c, obs_len_c;
  String fail_string1_1 = " * Incorrect type for '";
  String fail_string1_2 = " * Incorrect vector length for: '";
  String fail_string2 = "': expected ";
  String fail_string3 = "; got ";

  for (int i = 0; i < arg_len; ++i) {

    RObject arg = param_args[i];
    int arg_type = arg.sexp_type();
    if (arg_type == 13) arg_type = 14;  // May change this in future
    bool null_check = param_null2[i];
    bool arg_fail = false;
    arg_name = param_names[i];

    switch (expected_type) {
      case 10: exp_type = "`logical`";   break;
      case 14: exp_type = "`numeric`";   break;
      case 16: exp_type = "`character`"; break;
      case 25: exp_type = "`S4`";        break;
      default: stop("utils.cpp: unknown 'expected_type'");
    }

    switch (arg_type) {
      case  0: obs_type = "`NULL`";      break;
      case  1: obs_type = "`missing`";   break;
      case 10: obs_type = "`logical`";   break;
      case 14: obs_type = "`numeric`";   break;
      case 16: obs_type = "`character`"; break;
      case 25: obs_type = "`S4`";        break;
      default: obs_type = "`unknown`";   break;
    }

    if (null_check) {
      if (arg_type != 0 && arg_type != 1 && arg_type != expected_type) {
        fail_i = StringVector::create(fail_string1_1, arg_name, fail_string2,
            exp_type, fail_string3, obs_type);
        fails[i] = collapse(fail_i);
        arg_fail = true;
      }
    } else if (arg_type != expected_type) {
        fail_i = StringVector::create(fail_string1_1, arg_name, fail_string2,
            exp_type, fail_string3, obs_type);
        fails[i] = collapse(fail_i);
        arg_fail = true;
    }

    // S4 objects can't be length checked
    if (!arg_fail && arg_type != 0 && !null_check && arg_type != 25) {

      int arg_len;

      // {} are required for each case here or otherwise gives error at comp
      switch (arg_type) {
        case 10: {LogicalVector arg_ = as<LogicalVector>(arg);
                  arg_len = arg_.length();
                  break;}
        case 14: {NumericVector arg_ = as<NumericVector>(arg);
                  arg_len = arg_.length();
                  break;}
        case 16: {StringVector arg_ = as<StringVector>(arg);
                  arg_len = arg_.length();
                  break;}
        default: stop("utils.cpp: Unrecognised param type [INTERNAL ERROR]");
      }

      int exp_len = param_len2[i];
      // g++ compiler error with to_string! (but not with clang++)
      exp_len_c = std::to_string(exp_len);
      obs_len_c = std::to_string(arg_len);

      if (arg_len != exp_len && exp_len != 0) {
        fail_i = StringVector::create(fail_string1_2, arg_name, fail_string2,
            exp_len_c, fail_string3, obs_len_c);
        fails[((i + 1) * 2) - 1] = collapse(fail_i);
      }

    }

  }

  fails = clean_up_check(fails);
  return fails;

}
