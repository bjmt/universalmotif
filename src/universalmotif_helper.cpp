#include <Rcpp.h>
#include "utils.h"
using namespace Rcpp;

std::unordered_map<String, int> alphs_e = {
  {"DNA", 1}, {"RNA", 2},
  {"AA", 3}, {"custom", 4}
};

std::unordered_map<String, int> types2_e = {
  {"PCM", 1}, {"PPM", 2},
  {"PWM", 3}, {"ICM", 4}
};


StringVector universalmotif_alphabet(StringVector alphabet,
    NumericMatrix &m_motif) {

  // NOTE: assumes motif rows are properly alphabetically sorted

  switch (::alphs_e[alphabet[0]]) {
    case 1: rownames(m_motif) = ::dna;
            break;
    case 2: rownames(m_motif) = ::rna;
            break;
    case 3: rownames(m_motif) = ::aa;
            break;
    case 4: {
              StringVector mat_rownames = rownames(m_motif);
              if (mat_rownames.size() == 0)
                stop("Error creating universalmotif object; missing alphabet");
              mat_rownames = sort_unique(mat_rownames);
              alphabet[0] = collapse(mat_rownames);
              break;
             }
    default: {
               StringVector alph_split;
               for (R_xlen_t i = 0; i < alphabet[0].size(); ++i) {
                 alph_split.push_back(alphabet[0][i]);
               }
               if (alph_split.size() != m_motif.nrow())
                 stop("Alphabet length does not match matrix rows");
               rownames(m_motif) = sort_unique(alph_split);
               alphabet[0] = collapse(sort_unique(alph_split));
               break;
            }
  }

  return alphabet;

}

StringVector universalmotif_type(NumericMatrix &m_motif, StringVector type,
    const NumericVector &motif_colsums) {

  LogicalVector mat_1_check = m_motif >= 1;
  LogicalVector mat_0_check = m_motif == 0;
  LogicalVector mat_1_0_check = mat_1_check | mat_0_check;
  LogicalVector mat_099_101_check = (motif_colsums > 0.99) & (motif_colsums < 1.01);
  LogicalVector mat_pos_check = m_motif >= 0;

  if (StringVector::is_na(type[0]) || type.size() == 0) {
    if (is_true(all(mat_1_0_check)))
      type = StringVector::create("PCM");
    else if (is_true(all(mat_099_101_check)) && is_true(all(mat_pos_check)))
      type = StringVector::create("PPM");
    else if (is_true(all(mat_pos_check)))
      type = StringVector::create("ICM");
    else
      type = StringVector::create("PWM");
  } else if (type[0] == "PPM" && is_false(all(mat_099_101_check)) &&
      is_true(all(mat_pos_check))) {
    for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {
      m_motif(_, i) = m_motif(_, i) / motif_colsums[i];
    }
  }

  return type;

}

NumericVector universalmotif_nsites(NumericVector nsites, StringVector type,
    NumericMatrix &m_motif, const NumericVector &motif_colsums) {

  if (NumericVector::is_na(nsites[0]) || nsites.size() == 0) {

    if (type[0] == "PCM") nsites[0] = sum(m_motif(_, 1));
    else nsites = NumericVector::create();

  } else if (type[0] == "PCM" && is_true(any(motif_colsums != nsites[0]))) {

    double possum, fix;
    R_xlen_t tochange;

    for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {

      for (R_xlen_t j = 0; j < m_motif.nrow(); ++j) {
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

  return nsites;

}

NumericVector universalmotif_icscore(NumericVector icscore,
    NumericVector nsites, NumericMatrix m_motif, NumericVector bkg,
    StringVector type, double pseudocount) {

  R_xlen_t alph_len = m_motif.nrow();
  IntegerVector bkg_i = seq_len(alph_len) - 1;

  if (NumericVector::is_na(icscore[0]) || icscore.size() == 0) {
    double tmp_nsites = 0;
    if (nsites.size() != 0) tmp_nsites = nsites[0];
    NumericVector icscore_tmp(m_motif.ncol());
    for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {
      icscore_tmp[i] = position_icscoreC(m_motif(_, i), bkg[bkg_i], type[0],
          pseudocount, tmp_nsites);
    }
    icscore[0] = sum(icscore_tmp);
  }

  return icscore;

}

NumericVector universalmotif_bkg(NumericVector bkg, NumericMatrix m_motif) {

  // NOTE: Assumes the vector is already properly alphabetically sorted.

  R_xlen_t alph_len = m_motif.nrow();
  R_xlen_t bkg_len = bkg.size();

  if (NumericVector::is_na(bkg[0]) || bkg_len == 0) {

    bkg = rep(1.0 / m_motif.nrow(), m_motif.nrow());
    bkg.attr("names") = rownames(m_motif);
    return bkg;

  }

  SEXP bnames = bkg.attr("names");
  if (bkg_len == alph_len && Rf_isNull(bnames)) {

    bkg.attr("names") = rownames(m_motif);
    return bkg;

  }

  if (bkg_len < alph_len) stop("'bkg' vector is too short");

  return bkg;

}

StringVector universalmotif_consensus(NumericMatrix &m_motif, StringVector alphabet,
    StringVector type, double pseudocount, StringVector consensus) {

  StringVector consensus_tmp(m_motif.ncol());

  switch (::alphs_e[alphabet[0]]) {
    case 1: for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {
              consensus_tmp[i] = get_consensusC(m_motif(_, i), "DNA", type[0], pseudocount);
            }
            colnames(m_motif) = consensus_tmp;
            consensus = collapse(consensus_tmp);
            break;
    case 2: for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {
              consensus_tmp[i] = get_consensusC(m_motif(_, i), "RNA", type[0], pseudocount);
            }
            colnames(m_motif) = consensus_tmp;
            consensus = collapse(consensus_tmp);
            break;
    case 3: for (R_xlen_t i =0; i < m_motif.ncol(); ++i) {
              consensus_tmp[i] = get_consensusAAC(m_motif(_, i), type[0], pseudocount);
            }
            colnames(m_motif) = consensus_tmp;
            consensus = collapse(consensus_tmp);
            break;
    default: {
               StringVector mot_rownames = rownames(m_motif);
               colnames(m_motif) = StringVector::create();
               rownames(m_motif) = mot_rownames;
               consensus = StringVector::create();
               break;
             }
  }

  return consensus;

}

// [[Rcpp::export(rng = false)]]
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
    double pseudocount = 1.0,
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
  NumericVector motif_colsums = colSums(m_motif);
  // sometimes nsites can slip through (?) as nan (not R_NaN)
  if (std::isnan(nsites[0])) nsites = NumericVector::create();

  // name
  x.slot("name") = name;

  // altname
  if (!StringVector::is_na(altname[0]) && altname.length() > 0)
    x.slot("altname") = altname;

  // family
  if (!StringVector::is_na(family[0]) && family.length() > 0)
    x.slot("family") = family;

  // organism
  if (!StringVector::is_na(organism[0]) && organism.length() > 0)
    x.slot("organism") = organism;

  // alphabet (&m_motif)
  alphabet = universalmotif_alphabet(alphabet, m_motif);
  x.slot("alphabet") = alphabet;

  // type (&m_motif)
  type = universalmotif_type(m_motif, type, motif_colsums);
  x.slot("type") = type;

  // nsites (&m_motif)
  nsites = universalmotif_nsites(nsites, type, m_motif, motif_colsums);
  x.slot("nsites") = nsites;

  // pseudocount
  x.slot("pseudocount") = pseudocount;

  // bkg
  bkg = universalmotif_bkg(bkg, m_motif);
  x.slot("bkg") = bkg;

  // icscore
  icscore = universalmotif_icscore(icscore, nsites, m_motif, bkg, type, pseudocount);
  x.slot("icscore") = icscore[0];

  // bkgsites
  if (!NumericVector::is_na(bkgsites[0]))
    x.slot("bkgsites") = bkgsites;

  // consensus (&m_motif)
  consensus = universalmotif_consensus(m_motif, alphabet, type, pseudocount, consensus);
  x.slot("consensus") = consensus;

  // motif
  x.slot("motif") = m_motif;

  // strand
  x.slot("strand") = strand;

  // pval
  if (!NumericVector::is_na(pval[0]))
    x.slot("pval") = pval;

  // qval
  if (!NumericVector::is_na(qval[0]))
    x.slot("qval") = qval;

  // eval
  if (!NumericVector::is_na(eval[0]))
    x.slot("eval") = eval;

  // extrainfo
  if (!StringVector::is_na(extrainfo[0]))
    x.slot("extrainfo") = extrainfo;

  return x;

}

StringVector check_length(StringVector m_name, StringVector m_altname,
    StringVector m_family, StringVector m_organism, StringVector m_alphabet,
    StringVector m_type, NumericVector m_icscore, NumericVector m_nsites,
    NumericVector m_pseudocount, NumericVector m_bkgsites,
    StringVector m_consensus, StringVector m_strand, NumericVector m_pval,
    NumericVector m_qval, NumericVector m_eval, StringVector msg) {

  if (m_name.size()        != 1) msg.push_back("* name must be length 1");
  if (m_altname.size()      > 1) msg.push_back("* altname cannot be longer than 1");
  if (m_family.size()       > 1) msg.push_back("* family cannot be longer than 1");
  if (m_organism.size()     > 1) msg.push_back("* organism cannot be longer than 1");
  if (m_alphabet.size()    != 1) msg.push_back("* alphabet must be a single string");
  if (m_type.size()        != 1) msg.push_back("* type must be length 1");
  if (m_icscore.size()     != 1) msg.push_back("* icscore must be length 1");
  if (m_nsites.size()       > 1) msg.push_back("* nsites cannot be longer than 1");
  if (m_pseudocount.size() != 1) msg.push_back("* pseudocount must be length 1");
  if (m_bkgsites.size()     > 1) msg.push_back("* bkgsites cannot be longer than 1");
  if (m_consensus.size()    > 1) msg.push_back("* consensus cannot be longer than 1");
  if (m_strand.size()      != 1) msg.push_back("* strand must be length 1");
  if (m_pval.size()         > 1) msg.push_back("* pval cannot be longer than 1");
  if (m_qval.size()         > 1) msg.push_back("* qval cannot be longer than 1");
  if (m_eval.size()         > 1) msg.push_back("* eval cannot be longer than 1");

  return msg;

}

// [[Rcpp::export]]
bool check_bkg_names(StringVector alph, std::string blet) {

  LogicalVector failed(blet.size(), true);

  for (R_xlen_t i = 0; i < blet.length(); ++i) {

    for (R_xlen_t j = 0; j < alph.size(); ++j) {

      std::string alph_j = as<std::string>(alph[j]);

      if (alph_j[0] == blet[i]) {
        failed[i] = false;
        break;
      }

    }

  }

  return is_true(any(failed));

}

StringVector check_bkg(NumericVector bkg, StringVector alph, StringVector msg) {

  R_xlen_t blen = bkg.size();
  R_xlen_t alen = alph.size();

  SEXP bnames = bkg.attr("names");

  if (blen > alen) {
    if (Rf_isNull(bnames))
      msg.push_back("* bkg must be a named vector");
  }

  if (blen < alen)
    msg.push_back("* bkg vector length is too short");

  if (!Rf_isNull(bnames)) {

    StringVector bnames = bkg.names();
    bool zero_check = false;
    LogicalVector name_comp;
    StringVector alph_i;
    for (R_xlen_t i = 0; i < alen; ++i) {
      alph_i = rep(StringVector::create(alph[i]), blen);
      name_comp = bnames == alph_i;
      if (is_false(any(name_comp))) zero_check = true;
    }

    if (zero_check)

      msg.push_back("* bkg must contain 0-order possibilities for all letters");

    else {

      LogicalVector low_check = bkg < 0;
      LogicalVector high_check = bkg > 1;

      if (is_true(any(low_check)))
        msg.push_back("* bkg does not allow values less than 0");

      if (is_true(any(high_check)))
        msg.push_back("* bkg does not allow values higher than 1");

      NumericVector bkg_zero = bkg[alph];
      double bkg_sum = sum(bkg_zero);
      if (bkg_sum < 0.99 || bkg_sum > 1.01)
        msg.push_back("* 0-order bkg probabilities must add up to 1");

    }

    for (R_xlen_t i = 0; i < blen; ++i) {

      if (check_bkg_names(alph, as<std::string>(bnames[i]))) {
        msg.push_back("* unknown letters found in bkg names");
        break;
      }

    }

  }

  return msg;

}

StringVector check_char_slots(StringVector m_type, StringVector m_strand,
    StringVector msg) {

  SEXP s_type = m_type[0];
  if (Rf_isNull(s_type)) {
    msg.push_back("* type cannot be NULL");
    return msg;
  }

  if (m_type[0] != "PCM" &&
      m_type[0] != "PPM" &&
      m_type[0] != "PWM" &&
      m_type[0] != "ICM")
    msg.push_back("* type must be one of PCM, PPM, PWM, ICM");

  if (m_strand[0] != "+"  &&
      m_strand[0] != "-"  &&
      m_strand[0] != "+-" &&
      m_strand[0] != "-+")
    msg.push_back("* strand must be one of +, -, +-");

  return msg;

}

StringVector check_motif_and_type(NumericMatrix m_motif, StringVector m_type,
    NumericVector m_nsites, StringVector msg) {

  SEXP s_type = m_type[0];
  if (Rf_isNull(s_type)) return msg;

  NumericVector motif_colsums = colSums(m_motif);

  switch (::types2_e[m_type[0]]) {
    case 1: {
              if (m_nsites.size() > 0) {
                NumericVector colsums_unique = unique(m_nsites);
                if (colsums_unique.size() > 1)
                  msg.push_back("* for type PCM motif colSums must equal nsites");
              }
              LogicalVector int_check = (m_motif < 1.0) & (m_motif > 0.0);
              if (is_true(any(int_check)))
                msg.push_back("* type PCM motifs cannot contain values between 0 and 1");
              break;
            }
    case 2: {
              LogicalVector colsums_1_check = (motif_colsums > 0.99)
                & (motif_colsums < 1.01);
              if (is_false(all(colsums_1_check)))
                msg.push_back("* for type PPM colSums must equal 1");
              LogicalVector motif_pos_check = m_motif >= 0;
              if (is_false(all(motif_pos_check)))
                msg.push_back("* for type PPM only positive values are allowed");
              break;
            }
    case 4: if (is_true(any(m_motif < 0)))
              msg.push_back("* type ICM motifs cannot contain negative values");
            break;
  }

  return msg;

}

StringVector check_alphabet(NumericMatrix m_motif, StringVector m_alphabet,
    StringVector msg) {

  StringVector m_rownames = rownames(m_motif);

  switch (::alphs_e[m_alphabet[0]]) {
    case 1: {
              if (m_motif.nrow() != 4)
                msg.push_back("* DNA/RNA motifs must have 4 rows");
              LogicalVector rownames_check = m_rownames == ::dna;
              if (is_false(all(rownames_check)))
                msg.push_back("* rownames must be A, C, G, T");
              break;
            }
    case 2: {
              if (m_motif.nrow() != 4)
                msg.push_back("* DNA/RNA motifs must have 4 rows");
              LogicalVector rownames_check = m_rownames == ::rna;
              if (is_false(all(rownames_check)))
                msg.push_back("* rownames must be A, C, G, U");
              break;
            }
    case 3: {
              if (m_motif.nrow() != 20)
                msg.push_back("* AA motifs must have 20 rows");
              LogicalVector rownames_check = m_rownames == ::aa;
              if (is_false(all(rownames_check)))
                msg.push_back("* rownames must be ACDEFGHIKLMNPQRSTVWY");
              break;
            }
    default: {
               if (m_alphabet[0] == "custom") break;
               if (m_motif.nrow() != m_alphabet[0].size())
                 msg.push_back("* alphabet length does not match number of rows in motif");
               StringVector alph_split;
               for (R_xlen_t i = 0; i < m_alphabet[0].size(); ++i) {
                 alph_split.push_back(m_alphabet[0][i]);
               }
               LogicalVector rownames_check = sort_unique(alph_split) == m_rownames;
               if (is_false(all(rownames_check)))
                 msg.push_back("* rownames must match alphabet and be in alphabetical order");
               break;
            }
  }

  return msg;

}

StringVector check_consensus(StringVector m_consensus, NumericMatrix m_motif,
    StringVector msg) {

  if (m_consensus.size() > 0) {

    if (m_consensus[0].size() != m_motif.ncol()) {

      msg.push_back("* consensus string must have the same number of letters as motif positions");

    } else {

      StringVector consensus_split;
      StringVector motif_colnames = colnames(m_motif);

      for (R_xlen_t i = 0; i < m_consensus[0].size(); ++i) {
        consensus_split.push_back(m_consensus[0][i]);
        if (consensus_split[i] != motif_colnames[i])
          msg.push_back("* consensus string must match colnames");
      }

    }

  }

  return msg;

}

// [[Rcpp::export(rng = false)]]
StringVector validObject_universalmotif(S4 motif, bool throw_error = true) {

  StringVector msg;

  StringVector  m_name        = motif.slot("name");
  StringVector  m_altname     = motif.slot("altname");
  StringVector  m_family      = motif.slot("family");
  StringVector  m_organism    = motif.slot("organism");
  NumericMatrix m_motif       = motif.slot("motif");
  StringVector  m_alphabet    = motif.slot("alphabet");
  StringVector  m_type        = motif.slot("type");
  NumericVector m_icscore     = motif.slot("icscore");
  NumericVector m_nsites      = motif.slot("nsites");
  NumericVector m_pseudocount = motif.slot("pseudocount");
  NumericVector m_bkg         = motif.slot("bkg");
  NumericVector m_bkgsites    = motif.slot("bkgsites");
  StringVector  m_consensus   = motif.slot("consensus");
  StringVector  m_strand      = motif.slot("strand");
  NumericVector m_pval        = motif.slot("pval");
  NumericVector m_qval        = motif.slot("qval");
  NumericVector m_eval        = motif.slot("eval");

  // slot length checks
  msg = check_length(m_name, m_altname, m_family, m_organism, m_alphabet, m_type,
      m_icscore, m_nsites, m_pseudocount, m_bkgsites, m_consensus, m_strand,
      m_pval, m_qval, m_eval, msg);

  // character slot checks
  msg = check_char_slots(m_type, m_strand, msg);

  // bkg slot check
  msg = check_bkg(m_bkg, rownames(m_motif), msg);

  // check motif and type
  msg = check_motif_and_type(m_motif, m_type, m_nsites, msg);

  // check motif matches alphabet
  msg = check_alphabet(m_motif, m_alphabet, msg);

  // consensus check
  msg = check_consensus(m_consensus, m_motif, msg);

  if (msg.length() > 0) {

    msg = StringVector::create(all_checks_collapse(msg));

    if (throw_error) stop(as<std::string>(msg[0]));

  }

  return msg;

}

// [[Rcpp::export(rng = false)]]
DataFrame summarise_motifs_cpp(List motifs) {

  R_xlen_t n = motifs.size();
  StringVector name(n);
  StringVector altname(n);
  StringVector family(n);
  StringVector organism(n);
  StringVector alphabet(n);
  NumericVector icscore(n);
  IntegerVector nsites(n);
  IntegerVector bkgsites(n);
  StringVector consensus(n);
  StringVector strand(n);
  NumericVector pval(n);
  NumericVector qval(n);
  NumericVector eval(n);

  for (R_xlen_t i = 0; i < n; ++i) {

    S4 mot = motifs[i];

    StringVector name_i = mot.slot("name");
    name[i] = name_i[0];

    StringVector altname_i = mot.slot("altname");
    altname[i] = altname.size() == 1 ? altname_i[0] : NA_STRING;

    StringVector family_i = mot.slot("family");
    family[i] = family_i.size() == 1 ? family_i[0] : NA_STRING;

    StringVector organism_i = mot.slot("organism");
    organism[i] = organism_i.size() == 1 ? organism_i[0] : NA_STRING;

    StringVector alphabet_i = mot.slot("alphabet");
    alphabet[i] = alphabet_i[0];

    NumericVector icscore_i = mot.slot("icscore");
    icscore[i] = icscore_i[0];

    IntegerVector nsites_i = mot.slot("nsites");
    nsites[i] = nsites_i.size() == 1 ? nsites_i[0] : NA_INTEGER;

    IntegerVector bkgsites_i = mot.slot("bkgsites");
    bkgsites[i] = bkgsites_i.size() == 1 ? bkgsites_i[0] : NA_INTEGER;

    StringVector consensus_i = mot.slot("consensus");
    consensus[i] = consensus_i.size() == 1 ? consensus_i[0] : NA_STRING;

    StringVector strand_i = mot.slot("strand");
    strand[i] = strand_i[0];

    NumericVector pval_i = mot.slot("pval");
    pval[i] = pval_i.size() == 1 ? pval_i[0] : NA_REAL;

    NumericVector qval_i = mot.slot("qval");
    qval[i] = qval_i.size() == 1 ? qval_i[0] : NA_REAL;

    NumericVector eval_i = mot.slot("eval");
    eval[i] = eval_i.size() == 1 ? eval_i[0] : NA_REAL;

  }

  return DataFrame::create(

      _["name"]      = name,
      _["altname"]   = altname,
      _["family"]    = family,
      _["organism"]  = organism,
      _["alphabet"]  = alphabet,
      _["icscore"]   = icscore,
      _["nsites"]    = nsites,
      _["bkgsites"]  = bkgsites,
      _["consensus"] = consensus,
      _["strand"]    = strand,
      _["pval"]      = pval,
      _["qval"]      = qval,
      _["eval"]      = eval,

      _["stringsAsFactors"] = false);

}

// [[Rcpp::export(rng = false)]]
List universalmotif_to_list(S4 motif) {

  return List::create(

      _["name"]        = motif.slot("name"),
      _["altname"]     = motif.slot("altname"),
      _["family"]      = motif.slot("family"),
      _["organism"]    = motif.slot("organism"),
      _["motif"]       = motif.slot("motif"),
      _["alphabet"]    = motif.slot("alphabet"),
      _["type"]        = motif.slot("type"),
      _["icscore"]     = motif.slot("icscore"),
      _["nsites"]      = motif.slot("nsites"),
      _["pseudocount"] = motif.slot("pseudocount"),
      _["bkg"]         = motif.slot("bkg"),
      _["bkgsites"]    = motif.slot("bkgsites"),
      _["consensus"]   = motif.slot("consensus"),
      _["strand"]      = motif.slot("strand"),
      _["pval"]        = motif.slot("pval"),
      _["qval"]        = motif.slot ("qval"),
      _["eval"]        = motif.slot("eval"),
      _["multifreq"]   = motif.slot("multifreq"),
      _["extrainfo"]   = motif.slot("extrainfo")

      );

}
