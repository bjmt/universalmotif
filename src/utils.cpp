#include <Rcpp.h>
using namespace Rcpp;

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
