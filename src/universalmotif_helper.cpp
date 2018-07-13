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
    double pseudocount=0, double nsites=0) {

  int n_pos = position.size();
  int n_bkg = bkg.size();

  if (nsites == 0) nsites = 100;

  double n_pos2 = n_pos;
  NumericVector bkg2(n_pos, 1.0 / n_pos2);
  if (n_pos != n_bkg) bkg = bkg2;

  if (pseudocount != 0) {
    position = ppm_to_pcmC(position, nsites);
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
S4 universalmotif_cpp(
    NumericMatrix motif,
    String name = "new motif",
    StringVector altname = NA_STRING,
    StringVector family = NA_STRING,
    StringVector organism = NA_STRING,
    StringVector alphabet = "DNA",
    StringVector type = NA_STRING,
    NumericVector icscore = NumericVector::create(),
    NumericVector nsites = NumericVector::create(),
    double pseudocount = 0.8,
    NumericVector bkg = NumericVector::create(),
    NumericVector bkgsites = NumericVector::create(),
    StringVector consensus = NA_STRING,
    String strand = "+-",
    NumericVector pval = NumericVector::create(),
    NumericVector qval = NumericVector::create(),
    NumericVector eval = NumericVector::create(),
    StringVector extrainfo = NA_STRING) {

  S4 x("universalmotif");

  NumericMatrix m_motif = Rcpp::clone(motif);

  x.slot("name") = name;

  if (StringVector::is_na(altname[0]) || altname.length() == 0)
    altname = StringVector::create();
  x.slot("altname") = altname;

  if (StringVector::is_na(family[0]) || family.length() == 0)
    family = StringVector::create();
  x.slot("family") = family;

  if (StringVector::is_na(organism[0]) || organism.length() == 0)
    organism = StringVector::create();
  x.slot("organism") = organism;
  
  if (alphabet[0] == "DNA")
    rownames(m_motif) = CharacterVector::create("A", "C", "G", "T");
  else if (alphabet[0] == "RNA")
    rownames(m_motif) = CharacterVector::create("A", "C", "G", "U");
  else if (alphabet[0] == "AA")
    rownames(m_motif) = CharacterVector::create("A", "C", "D", "E", "F",
                                                "G", "H", "I", "K", "L",
                                                "M", "N", "P", "Q", "R",
                                                "S", "T", "V", "W", "Y");
  else if (alphabet[0] != "custom") {
    StringVector alph_split;
    for (int i = 0; i < alphabet[0].size(); ++i) {
      alph_split.push_back(alphabet[0][i]);
    }
    StringVector mat_rownames = rownames(m_motif);
    rownames(m_motif) = sort_unique(alph_split);
  }

  if (StringVector::is_na(alphabet[0]) || alphabet.length() == 0) {
    switch(m_motif.nrow()) {
      case 4: alphabet = "DNA";
              break;
      case 20: alphabet = "AA";
               break;
      default: alphabet = "custom";
               break;
    }
  }
  x.slot("alphabet") = alphabet;

  NumericVector motif_colsums = colSums(m_motif);
  LogicalVector mat_1_check = m_motif >= 1;
  LogicalVector mat_0_check = m_motif == 0;
  LogicalVector mat_1_0_check = mat_1_check | mat_0_check;
  LogicalVector mat_099_101_check = motif_colsums > 0.99 & motif_colsums < 1.01;
  LogicalVector mat_pos_check = m_motif >= 0;
  if (StringVector::is_na(type[0]) || type.length() == 0) {
    if (is_true(all(mat_1_0_check)))
      type = "PCM";
    else if (is_true(all(mat_099_101_check)) && is_true(all(mat_pos_check)))
      type = "PPM";
    else if (is_true(all(mat_pos_check)))
      type = "ICM";
    else
      type = "PWM";
  } else if (type[0] == "PPM" && is_false(all(mat_099_101_check)) &&
      is_true(all(mat_pos_check))) {
    for (int i = 0; i < m_motif.ncol(); ++i) {
      for (int j = 0; j < m_motif.nrow(); ++j) {
        m_motif(j, i) = m_motif(j, i) / motif_colsums[i];
      }
    }
  }
  x.slot("type") = type;
  
  if (NumericVector::is_na(nsites[0]) || nsites.length() == 0) {
    if (type[0] == "PCM") nsites[0] = sum(m_motif(_, 1));
    else nsites = NumericVector::create();
  } else if (type[0] == "PCM" && is_true(any(motif_colsums != nsites[0]))) {
    double possum, fix;
    int tochange;
    for (int i = 0; i < m_motif.ncol(); ++i) {
      for (int j = 0; j < motif.nrow(); ++j) {
        m_motif(j, i) = m_motif(j, i) / motif_colsums[i];  
        m_motif(j, i) = round(m_motif(j, i) * nsites[0]);
      }
      possum = sum(m_motif(_, i));
      if (possum != nsites[0]) {
        fix = nsites[0] - possum;
        tochange = which_max(m_motif(_, i));
        m_motif(tochange, i) += fix;
      }
    }
  }
  x.slot("nsites") = nsites;
  
  x.slot("pseudocount") = pseudocount;
  
  if (NumericVector::is_na(bkg[0]) || bkg.length() == 0) {
    bkg = rep(1.0 / m_motif.nrow(), m_motif.nrow());
  }
  x.slot("bkg") = bkg;
  
  if (NumericVector::is_na(icscore[0]) || icscore.length() == 0) {
    double tmp_nsites = 0;
    if (nsites.length() != 0) tmp_nsites = nsites[0];
    NumericVector icscore_tmp(m_motif.ncol());
    for (int i = 0; i < m_motif.ncol(); ++i) {
      icscore_tmp[i] = position_icscoreC(m_motif(_, i), bkg, type[0], pseudocount,
          tmp_nsites);
    }
    icscore[0] = sum(icscore_tmp);
  }
  x.slot("icscore") = icscore[0];

  if (NumericVector::is_na(bkgsites[0])) bkgsites = NumericVector::create();
  x.slot("bkgsites") = bkgsites;
  
  StringVector consensus_tmp(m_motif.ncol());
  if (alphabet[0] == "DNA") {
    for (int i = 0; i < m_motif.ncol(); ++i) {
      consensus_tmp[i] = get_consensusC(m_motif(_, i), "DNA", type[0], pseudocount);
    }
    x.slot("consensus") = collapse(consensus_tmp);
    colnames(m_motif) = consensus_tmp;
  } else if (alphabet[0] == "RNA") {
    for (int i = 0; i < motif.ncol(); ++i) {
      consensus_tmp[i] = get_consensusC(m_motif(_, i), "RNA", type[0], pseudocount);
    }
    x.slot("consensus") = collapse(consensus_tmp);
    colnames(m_motif) = consensus_tmp;
  } else if (alphabet[0] == "AA") {
    for (int i =0; i < m_motif.ncol(); ++i) {
      consensus_tmp[i] = get_consensusAAC(m_motif(_, i), type[0], pseudocount);
    }
    x.slot("consensus") = collapse(consensus_tmp);
    colnames(m_motif) = consensus_tmp;
  } else {
    StringVector mot_rownames = rownames(m_motif);
    colnames(m_motif) = StringVector::create();
    rownames(m_motif) = mot_rownames;
    x.slot("consensus") = StringVector::create();
  }
  
  x.slot("motif") = m_motif;
  
  x.slot("strand") = strand;
  
  if (NumericVector::is_na(pval[0])) pval = NumericVector::create();
  x.slot("pval") = pval;
  
  if (NumericVector::is_na(qval[0])) qval = NumericVector::create();
  x.slot("qval") = qval;
  
  if (NumericVector::is_na(eval[0])) eval = NumericVector::create();
  x.slot("eval") = eval;
  
  if (StringVector::is_na(extrainfo[0])) extrainfo = StringVector::create();
  x.slot("extrainfo") = extrainfo;

  return x;

}

// [[Rcpp::export]]
StringVector validObject_universalmotif(S4 motif) {

  StringVector msg;

  StringVector m_name = motif.slot("name");
  StringVector m_altname = motif.slot("altname");
  StringVector m_family = motif.slot("family");
  StringVector m_organism = motif.slot("organism");
  NumericMatrix m_motif = motif.slot("motif");
  StringVector m_alphabet = motif.slot("alphabet");
  StringVector m_type = motif.slot("type");
  NumericVector m_icscore = motif.slot("icscore");
  NumericVector m_nsites = motif.slot("nsites");
  NumericVector m_pseudocount = motif.slot("pseudocount");
  NumericVector m_bkg = motif.slot("bkg");
  NumericVector m_bkgsites = motif.slot("bkgsites");
  StringVector m_consensus = motif.slot("consensus");
  StringVector m_strand = motif.slot("strand");
  NumericVector m_pval = motif.slot("pval");
  NumericVector m_qval = motif.slot("qval");
  NumericVector m_eval = motif.slot("eval");

  // slot length checks
  
  if (m_name.length() != 1) msg.push_back("name must be length 1");
  if (m_altname.length() > 1) msg.push_back("altname cannot be longer than 1");
  if (m_family.length() > 1) msg.push_back("family cannot be longer than 1");
  if (m_organism.length() > 1) msg.push_back("organism cannot be longer than 1");
  if (m_alphabet.length() != 1) msg.push_back("alphabet must be a single string");
  if (m_type.length() != 1) msg.push_back("type must be length 1");
  if (m_icscore.length() != 1) msg.push_back("icscore must be length 1");
  if (m_nsites.length() > 1) msg.push_back("nsites cannot be longer than 1");
  if (m_pseudocount.length() != 1) msg.push_back("pseudocount must be length 1");
  if (m_bkgsites.length() > 1) msg.push_back("bkgsites cannot be longer than 1");
  if (m_consensus.length() > 1) msg.push_back("consensus cannot be longer than 1");
  if (m_strand.length() != 1) msg.push_back("strand must be length 1");
  if (m_pval.length() > 1) msg.push_back("pval cannot be longer than 1");
  if (m_qval.length() > 1) msg.push_back("qval cannot be longer than 1");
  if (m_eval.length() > 1) msg.push_back("eval cannot be longer than 1");

  // character slot checks
  if (m_type[0] != "PCM" && m_type[0] != "PPM" && m_type[0] != "PWM" &&
      m_type[0] != "ICM") msg.push_back("type must be one of PCM, PPM, PWM, ICM");
  if (m_strand[0] != "+" && m_strand[0] != "-" && m_strand[0] != "+-" &&
      m_strand[0] != "-+") msg.push_back("strand must be one of +, -, +-");

  // bkg slot check
  if (m_bkg.length() != m_motif.nrow())
    msg.push_back("bkg length must be equal to alphabet size");

  // check motif and type
  NumericVector motif_colsums = colSums(m_motif);

  if (m_type[0] == "PCM" && m_nsites.length() > 0) {
    NumericVector colsums_unique = unique(m_nsites);
    if (colsums_unique.length() > 1) msg.push_back("motif colSums must equal nsites");
  } else if (m_type[0] == "PPM") {
    LogicalVector colsums_1_check = motif_colsums > 0.99 & motif_colsums < 1.01;
    if (is_false(all(colsums_1_check)))
      msg.push_back("for type PPM colSums must equal 1");
    LogicalVector motif_pos_check = m_motif >= 0;
    if (is_false(all(motif_pos_check)))
      msg.push_back("for type PPM only positive values are allowed");
  }

  // check motif matches alphabet
  StringVector m_rownames = rownames(m_motif);
  if (m_alphabet[0] == "DNA") {

    if (m_motif.nrow() != 4) msg.push_back("DNA/RNA motifs must have 4 rows");
    StringVector ref_rownames = StringVector::create("A", "C", "G", "T");
    LogicalVector rownames_check = m_rownames == ref_rownames;
    if (is_false(all(rownames_check))) msg.push_back("rownames must be A, C, G, T");

  } else if (m_alphabet[0] == "RNA") {

    if (m_motif.nrow() != 4) msg.push_back("DNA/RNA motifs must have 4 rows");
    StringVector ref_rownames = StringVector::create("A", "C", "G", "U");
    LogicalVector rownames_check = m_rownames == ref_rownames;
    if (is_false(all(rownames_check))) msg.push_back("rownames must be A, C, G, U");

  } else if (m_alphabet[0] == "AA") {

    if (m_motif.nrow() != 20) msg.push_back("AA motifs must have 20 rows");
    StringVector ref_rownames = StringVector::create("A", "C", "D", "E", "F",
                                                     "G", "H", "I", "K", "L",
                                                     "M", "N", "P", "Q", "R",
                                                     "S", "T", "V", "W", "Y");
    LogicalVector rownames_check = m_rownames == ref_rownames;
    if (is_false(all(rownames_check)))
      msg.push_back("rownames must be A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y");

  } else if (m_alphabet[0] != "custom") {

    if (m_motif.nrow() != m_alphabet[0].size())
      msg.push_back("alphabet length does not match number of rows in motif");
    StringVector alph_split;
    for (int i = 0; i < m_alphabet[0].size(); ++i) {
      alph_split.push_back(m_alphabet[0][i]);
    }
    LogicalVector rownames_check = sort_unique(alph_split) == sort_unique(m_rownames);
    if (is_false(all(rownames_check))) msg.push_back("rownames must match alphabet");

  }

  // consensus check
  
  if (m_consensus.length() > 0) {
    if (m_consensus[0].size() != m_motif.ncol())
      msg.push_back("consensus string must have the same number of letters as motif positions");
    StringVector consensus_split;
    StringVector motif_colnames = colnames(m_motif);
    for (int i = 0; i < m_consensus[0].size(); ++i) {
      consensus_split.push_back(m_consensus[0][i]);
      if (consensus_split[i] != motif_colnames[i])
        msg.push_back("consensus string must match colnames");
    }
  }

  return msg;

}
