#include <Rcpp.h>
#include <sstream>  // needed for `to_string` work around
using namespace Rcpp;

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

// [[Rcpp::export]]
void print_pb(int out) {
       if (out == 0)               Rprintf("   0%%");
  else if (out > 0 && out < 10)    Rprintf("\b\b\b\b  %i%%", out);
  else if (out >= 10 && out < 100) Rprintf("\b\b\b\b %i%%", out);
  else if (out == 100)             Rprintf("\b\b\b\b%i%%", out);
  else if (out == -1)              Rprintf("\b\b\b\b100%%\n");
}

// [[Rcpp::export]]
void update_pb(int i, int max) {

  int prev = i - 1;
  int out = 100 * i / max;
  if (i == max) out = -1;

  if (prev > 0 && out != -1) {
    prev = 100 * prev / max;
    if (prev != out) print_pb(out);
  } else if (out == -1) {
    print_pb(out);
  }

}

// [[Rcpp::export]]
StringVector strsplit_cpp(std::string x) {  // slightly slower than
  int n = x.size();                         // strsplit(x, "")[[1]]
  StringVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = x.substr(i, 1);
  }
  return out;
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
String collapse_cpp(StringVector x) {
  // collapse_cpp(x) is about 3 times faster than base::paste(x, collapse = "")
  // collapse_cpp(c(x, y)) about 2 times faster than base::paste0(x, y)
  return collapse(x);
}

// [[Rcpp::export]]
StringVector sample_string_cpp(StringVector x, int size, bool replace = false,
    sugar::probs_t prob = R_NilValue) {
  // about 1.5 times faster than base::sample
  return sample(x, size, replace, prob);
}

// NumericVector sample_numeric_cpp(NumericVector x, int size, bool replace = false,
    // sugar::probs_t prob = R_NilValue) {
  // // about 1.25 times faster than base::sample
  // return sample(x, size, replace, prob);
// }

// [[Rcpp::export]]
NumericVector pcm_to_ppmC(NumericVector position, double pseudocount=0) {

  double possum = sum(position);
  int num = position.size();
  NumericVector out(num);

  if (pseudocount != 0) {
    for (int i = 0; i < num; ++i) {
      out[i] = (position[i] + (pseudocount / num)) / (possum + pseudocount);
    }
  } else {
    for (int i = 0; i < num; ++i) {
      out[i] = position[i] / possum;
    }
  }

  return out;

}

// [[Rcpp::export]]
NumericVector ppm_to_pcmC(NumericVector position, double nsites) {

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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
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

  for (int i = 0; i < n_pos; ++i) {
    position[i] *= bkg[i];
  }

  possum = sum(position);
  if (possum > 0.99 && possum < 1.01) return position;

  possum = sum(position);
  for (int i = 0; i < n_pos; ++i) {
    position[i] /= possum;
  }

  return position;

}

// [[Rcpp::export]]
NumericVector ppm_to_icmC(NumericVector position, NumericVector bkg=0,
    bool relative_entropy=false) {

  int n_pos = position.size();
  int n_bkg = bkg.size();

  double n_pos2 = n_pos;
  NumericVector bkg2(n_pos, 1.0 / n_pos2);
  if (n_pos != n_bkg) bkg = bkg2;

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
    NumericVector ic(1, total_ic);
    for (int i = 0; i < n_pos; ++i) {
      position[i] = position[i] * total_ic;
    }
    return position;
  }

}

// [[Rcpp::export]]
double position_icscoreC(NumericVector position, NumericVector bkg=0,
    String type="PPM", double pseudocount=0.8, double nsites=100,
    bool relative_entropy=false) {

  int n_pos = position.size();
  int n_bkg = bkg.size();

  double n_pos2 = n_pos;
  NumericVector bkg2(n_pos, 1.0 / n_pos2);
  if (n_pos != n_bkg) bkg = bkg2;

  if (type == "PCM") position = pcm_to_ppmC(position, pseudocount);
  else if (type == "PWM") position = pwm_to_ppmC(position, bkg);
  else if (type == "ICM") return sum(position);
  else if (type == "PPM") {
    position = ppm_to_pcmC(position, nsites);
    position = pcm_to_ppmC(position, pseudocount);
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

// [[Rcpp::export]]
NumericVector icm_to_ppmC(NumericVector position) {
  double total_ic = sum(position);
  int n = position.size();
  for (int i = 0; i < n; ++i) {
    position[i] /= total_ic;
  }
  return position;
}

// [[Rcpp::export]]
String get_consensusC(NumericVector position, String alphabet="DNA",
    String type="PPM", double pseudocount=0.8) {

  if (type == "PCM") {
    position = pcm_to_ppmC(position, pseudocount);
  } else if (type == "PWM") {
    position = pwm_to_ppmC(position);
  } else if (type == "ICM") {
    position = icm_to_ppmC(position);
  }


  // single letter consensus
       if (position[0] > 0.5 &&
           position[0] > position[1] * 2.0 &&
           position[0] > position[2] * 2.0 &&
           position[0] > position[3] * 2.0) return "A";
  else if (position[1] > 0.5 &&
           position[1] > position[0] * 2.0 &&
           position[1] > position[2] * 2.0 &&
           position[1] > position[3] * 2.0) return "C";
  else if (position[2] > 0.5 &&
           position[2] > position[0] * 2.0 &&
           position[2] > position[1] * 2.0 &&
           position[2] > position[3] * 2.0) return "G";
  else if (position[3] > 0.5 &&
           position[3] > position[0] * 2.0 &&
           position[3] > position[1] * 2.0 &&
           position[3] > position[2] * 2.0) {
    if (alphabet == "DNA") return "T";
    else return "U";
  }

  // two letter consensus
  else if (position[0] > 0.5) {
         if (position[1] > 0.25) return "M";
    else if (position[2] > 0.25) return "R";
    else if (position[3] > 0.25) return "W";
  }
  else if (position[1] > 0.5) {
         if (position[0] > 0.25) return "M";
    else if (position[2] > 0.25) return "S";
    else if (position[3] > 0.25) return "Y";
  }
  else if (position[2] > 0.5) {
         if (position[0] > 0.25) return "R";
    else if (position[1] > 0.25) return "S";
    else if (position[3] > 0.25) return "Y";
  }
  else if (position[3] > 0.5) {
         if (position[0] > 0.25) return "W";
    else if (position[1] > 0.25) return "Y";
    else if (position[2] > 0.25) return "K";
  }
  else if ((position[0] + position[1]) > 0.75) return "M";
  else if ((position[0] + position[2]) > 0.75) return "R";
  else if ((position[0] + position[3]) > 0.75) return "W";
  else if ((position[1] + position[2]) > 0.75) return "S";
  else if ((position[1] + position[3]) > 0.75) return "Y";
  else if ((position[2] + position[3]) > 0.75) return "K";

  // three letter consensus
  else if (position[0] > 0.25 &&
           position[1] > 0.25 &&
           position[3] > 0.25) return "H";
  else if (position[1] > 0.25 &&
           position[2] > 0.25 &&
           position[3] > 0.25) return "B";
  else if (position[0] > 0.25 &&
           position[1] > 0.25 &&
           position[2] > 0.25) return "V";
  else if (position[0] > 0.25 &&
           position[2] > 0.25 &&
           position[3] > 0.25) return "D";

  // no consensus
  return "N";

}

// [[Rcpp::export]]
NumericVector consensus_to_ppmC(String letter) {

       if (letter == "A") return NumericVector::create(0.997, 0.001, 0.001, 0.001);
  else if (letter == "C") return NumericVector::create(0.001, 0.997, 0.001, 0.001);
  else if (letter == "G") return NumericVector::create(0.001, 0.001, 0.997, 0.001);
  else if (letter == "T") return NumericVector::create(0.001, 0.001, 0.001, 0.997);
  else if (letter == "U") return NumericVector::create(0.001, 0.001, 0.001, 0.997);
  else if (letter == "R") return NumericVector::create(0.499, 0.001, 0.499, 0.001);
  else if (letter == "Y") return NumericVector::create(0.001, 0.499, 0.001, 0.449);
  else if (letter == "M") return NumericVector::create(0.499, 0.499, 0.001, 0.001);
  else if (letter == "K") return NumericVector::create(0.001, 0.001, 0.499, 0.499);
  else if (letter == "S") return NumericVector::create(0.001, 0.499, 0.499, 0.001);
  else if (letter == "W") return NumericVector::create(0.499, 0.001, 0.001, 0.499);
  else if (letter == "H") return NumericVector::create(0.333, 0.333, 0.001, 0.333);
  else if (letter == "B") return NumericVector::create(0.001, 0.333, 0.333, 0.333);
  else if (letter == "V") return NumericVector::create(0.333, 0.333, 0.333, 0.001);
  else if (letter == "D") return NumericVector::create(0.333, 0.001, 0.333, 0.333);
  else if (letter == "N") return NumericVector::create( 0.25,  0.25,  0.25,  0.25);
  else if (letter == "+") return NumericVector::create( 0.25,  0.25,  0.25,  0.25);
  else if (letter == "-") return NumericVector::create( 0.25,  0.25,  0.25,  0.25);
  else if (letter == ".") return NumericVector::create( 0.25,  0.25,  0.25,  0.25);

  return NumericVector::create(0.0, 0.0, 0.0, 0.0);

}

// [[Rcpp::export]]
NumericVector consensus_to_ppmAAC(String letter) {

  NumericVector let_n(20, 0.05);
       if (letter == "X") return let_n;
  else if (letter == ".") return let_n;
  else if (letter == "-") return let_n;
  else if (letter == "+") return let_n;

  NumericVector let_2(20, 0.001);
  if (letter == "B") {
    let_2[2] = 0.491;
    let_2[11] = 0.491;
    return let_2;
  } else if (letter == "Z") {
    let_2[3] = 0.491;
    let_2[13] = 0.491;
    return let_2;
  } else if (letter == "J") {
    let_2[7] = 0.491;
    let_2[9] = 0.491;
    return let_2;
  }

  let_2.names() = CharacterVector::create("A", "C", "D", "E", "F",
                                          "G", "H", "I", "K", "L",
                                          "M", "N", "P", "Q", "R",
                                          "S", "T", "V", "W", "Y");
  int n = let_2.findName(letter);
  let_2[n] = 0.981;
  return let_2;

}

// [[Rcpp::export]]
String get_consensusAAC(NumericVector position, String type="PPM",
    double pseudocount=0.0) {

  if (type == "PCM") {
    position = pcm_to_ppmC(position, pseudocount);
  } else if (type == "PWM") {
    position = pwm_to_ppmC(position);
  } else if (type == "ICM") {
    position = icm_to_ppmC(position);
  }

       if (position[2] >= 0.4 && position[11] >= 0.4) return "B";
  else if (position[3] >= 0.4 && position[13] >= 0.4) return "Z";
  else if (position[7] >= 0.4 && position[9]  >= 0.4) return "J";

  bool test = all(position < 0.1).is_true();
  if (test) return "X";

  NumericVector position2 = Rcpp::clone(position);
  std::sort(position2.begin(), position2.end());
  if (position2[19] == position2[18]) return "X";

  CharacterVector lets = CharacterVector::create("A", "C", "D", "E", "F",
                                                 "G", "H", "I", "K", "L",
                                                 "M", "N", "P", "Q", "R",
                                                 "S", "T", "V", "W", "Y");

  return lets[which_max(position)];

}

// [[Rcpp::export]]
StringVector clean_up_check(StringVector fails) {

  int fails_len = fails.length();
  LogicalVector fails_keep(fails_len, true);

  for (int i = 0; i < fails_len; ++i) {
    if (fails[i] == "") fails_keep[i] = false;
  }

  return fails[fails_keep];

}

// [[Rcpp::export]]
StringVector check_fun_params(List param_args, IntegerVector param_len,
    LogicalVector param_null, String expected_type_string) {

  // This function WILL fail if param_args is not a NAMED list

  int expected_type;
  if (expected_type_string == "character") expected_type = 16;
  else if (expected_type_string == "numeric") expected_type = 14;
  else if (expected_type_string == "logical") expected_type = 10;
  else if (expected_type_string == "S4") expected_type = 25;
  else stop("Unknown 'expected_type_string'");
  int arg_len = param_args.length();
  StringVector fails(arg_len * 2);
  StringVector param_names = param_args.names();

  IntegerVector param_len2;
  LogicalVector param_null2;

  if (param_len.length() == 0) {
    param_len2 = rep(1, arg_len);
  } else param_len2 = param_len;

  if (param_null.length() == 0) {
    param_null2 = rep(false, arg_len);
  } else param_null2 = param_null;

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
      // compiler error with to_string!
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
