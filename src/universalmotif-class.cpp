#include <Rcpp.h>
#include <unordered_map>
#include "utils-exported.h"
#include "utils-internal.h"
#include "types.h"

std::unordered_map<Rcpp::String, int> ALPHS_e = {
  {"DNA",    1},
  {"RNA",    2},
  {"AA",     3},
  {"custom", 4}
};

enum SEQUENCE_ALPHS {

  DNAe    = 1,
  RNAe    = 2,
  AAe     = 3,
  CUSTOMe = 4

};

std::unordered_map<Rcpp::String, int> TYPES2_e = {
  {"PCM", 1},
  {"PPM", 2},
  {"PWM", 3},
  {"ICM", 4}
};

enum MATRIX_TYPES {

  PCMe = 1,
  PPMe = 2,
  PWMe = 3,
  ICMe = 4

};

Rcpp::StringVector universalmotif_alphabet(Rcpp::StringVector alphabet,
    Rcpp::NumericMatrix &m_motif) {

  // NOTE: assumes motif rows are properly alphabetically sorted

  switch (::ALPHS_e[alphabet[0]]) {
    case DNAe: Rcpp::rownames(m_motif) = ::DNA;
            break;
    case RNAe: Rcpp::rownames(m_motif) = ::RNA;
            break;
    case AAe: Rcpp::rownames(m_motif) = ::AMINOACIDS;
            break;
    case CUSTOMe: {
              Rcpp::StringVector mat_rownames = Rcpp::rownames(m_motif);
              if (mat_rownames.size() == 0)
                Rcpp::stop("Error creating universalmotif object; missing alphabet");
              mat_rownames = Rcpp::sort_unique(mat_rownames);
              alphabet[0] = Rcpp::collapse(mat_rownames);
              break;
             }
    default: {
               Rcpp::StringVector alph_split;
               for (R_xlen_t i = 0; i < alphabet[0].size(); ++i) {
                 alph_split.push_back(alphabet[0][i]);
               }
               if (alph_split.size() != m_motif.nrow())
                 Rcpp::stop("Alphabet length does not match matrix rows");
               Rcpp::rownames(m_motif) = Rcpp::sort_unique(alph_split);
               alphabet[0] = Rcpp::collapse(Rcpp::sort_unique(alph_split));
               break;
            }
  }

  return alphabet;

}

Rcpp::StringVector universalmotif_type(Rcpp::NumericMatrix &m_motif,
    Rcpp::StringVector type, const Rcpp::NumericVector &motif_colsums) {

  Rcpp::LogicalVector mat_1_check = m_motif >= 1;
  Rcpp::LogicalVector mat_0_check = m_motif == 0;
  Rcpp::LogicalVector mat_1_0_check = mat_1_check | mat_0_check;
  Rcpp::LogicalVector mat_099_101_check = (motif_colsums > 0.99) & (motif_colsums < 1.01);
  Rcpp::LogicalVector mat_pos_check = m_motif >= 0;

  if (Rcpp::StringVector::is_na(type[0]) || type.size() == 0) {
    if (Rcpp::is_true(Rcpp::all(mat_1_0_check)))
      type = Rcpp::StringVector::create("PCM");
    else if (Rcpp::is_true(Rcpp::all(mat_099_101_check))
             && Rcpp::is_true(Rcpp::all(mat_pos_check)))
      type = Rcpp::StringVector::create("PPM");
    else if (Rcpp::is_true(Rcpp::all(mat_pos_check)))
      type = Rcpp::StringVector::create("ICM");
    else
      type = Rcpp::StringVector::create("PWM");
  } else if (type[0] == "PPM" && Rcpp::is_false(Rcpp::all(mat_099_101_check)) &&
      Rcpp::is_true(Rcpp::all(mat_pos_check))) {
    for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {
      m_motif(Rcpp::_, i) = m_motif(Rcpp::_, i) / motif_colsums[i];
    }
  }

  return type;

}

Rcpp::NumericVector universalmotif_nsites(Rcpp::NumericVector nsites,
    const Rcpp::StringVector &type, Rcpp::NumericMatrix &m_motif,
    const Rcpp::NumericVector &motif_colsums) {

  if (Rcpp::NumericVector::is_na(nsites[0]) || nsites.size() == 0) {

    if (type[0] == "PCM") nsites[0] = Rcpp::sum(m_motif(Rcpp::_, 0));
    else nsites = Rcpp::NumericVector::create();

  } else if (type[0] == "PCM" && Rcpp::is_true(Rcpp::any(motif_colsums != nsites[0]))) {

    double possum, fix;
    R_xlen_t tochange;

    for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {

      for (R_xlen_t j = 0; j < m_motif.nrow(); ++j) {
        m_motif(j, i) = m_motif(j, i) / motif_colsums[i];
        m_motif(j, i) = round(m_motif(j, i) * nsites[0]);
      }

      possum = Rcpp::sum(m_motif(Rcpp::_, i));
      if (possum != nsites[0]) {
        fix = nsites[0] - possum;
        tochange = Rcpp::which_max(m_motif(Rcpp::_, i));
        m_motif(tochange, i) += fix;
      }

    }

  }

  return nsites;

}

Rcpp::NumericVector universalmotif_icscore(Rcpp::NumericVector icscore,
    const Rcpp::NumericVector &nsites, const Rcpp::NumericMatrix &m_motif,
    const Rcpp::NumericVector &bkg, const Rcpp::StringVector &type,
    double pseudocount) {

  R_xlen_t alph_len = m_motif.nrow();
  Rcpp::IntegerVector bkg_i = Rcpp::seq_len(alph_len) - 1;

  if (Rcpp::NumericVector::is_na(icscore[0]) || icscore.size() == 0) {
    double tmp_nsites = 0;
    if (nsites.size() != 0) tmp_nsites = nsites[0];
    Rcpp::NumericVector icscore_tmp(m_motif.ncol());
    for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {
      Rcpp::NumericVector tmp = m_motif(Rcpp::_, i);
      icscore_tmp[i] = position_icscoreC(Rcpp::as<std::vector<double>>(tmp),
          Rcpp::as<std::vector<double>>(bkg[bkg_i]), Rcpp::as<std::string>(type[0]),
          pseudocount, tmp_nsites);
    }
    icscore[0] = Rcpp::sum(icscore_tmp);
  }

  return icscore;

}

Rcpp::NumericVector reorder_named_num_vec_cpp(const Rcpp::NumericVector x,
    const Rcpp::IntegerVector index) {
  if (x.size() != index.size()) {
    Rcpp::stop("[reorder_named_num_vec_cpp] x.size() != index.size()");
  }
  SEXP x_names = x.attr("names");
  if (Rf_isNull(x_names)) {
    Rcpp::stop("[reorder_named_num_vec_cpp] x is not named");
  }
  Rcpp::CharacterVector x_names2 = Rcpp::as<Rcpp::CharacterVector>(x_names);
  Rcpp::CharacterVector x_names3(x_names2.size());
  Rcpp::NumericVector x2(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    int j = index[i] - 1;
    x_names3[i] = x_names2[j];
    x2[i] = x[j];
  }
  x2.attr("names") = x_names3;
  return x2;
}

Rcpp::NumericVector universalmotif_bkg(Rcpp::NumericVector bkg,
    const Rcpp::NumericMatrix &m_motif) {

  R_xlen_t alph_len = m_motif.nrow();
  R_xlen_t bkg_len = bkg.size();

  if (Rcpp::NumericVector::is_na(bkg[0]) || bkg_len == 0) {

    bkg = Rcpp::rep(1.0 / m_motif.nrow(), m_motif.nrow());
    bkg.attr("names") = Rcpp::rownames(m_motif);
    return bkg;

  }

  SEXP bnames = bkg.attr("names");
  if (bkg_len == alph_len && Rf_isNull(bnames)) {

    bkg.attr("names") = Rcpp::rownames(m_motif);
    return bkg;

  }

  Rcpp::IntegerVector bkg_order = order_char_cpp(bkg.attr("names"));
  bkg = reorder_named_num_vec_cpp(bkg, bkg_order);

  if (bkg_len < alph_len) Rcpp::stop("'bkg' vector is too short");

  return bkg;

}

Rcpp::StringVector universalmotif_consensus(Rcpp::NumericMatrix &m_motif,
    const Rcpp::StringVector &alphabet, const Rcpp::StringVector &type,
    double pseudocount, Rcpp::StringVector consensus) {

  Rcpp::StringVector consensus_tmp(m_motif.ncol());

  switch (::ALPHS_e[alphabet[0]]) {
    case DNAe: for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {
              Rcpp::NumericVector tmp = m_motif(Rcpp::_, i);
              consensus_tmp[i] = get_consensusC(Rcpp::as<std::vector<double>>(tmp),
                  "DNA", Rcpp::as<std::string>(type[0]), pseudocount);
            }
            Rcpp::colnames(m_motif) = consensus_tmp;
            consensus = Rcpp::collapse(consensus_tmp);
            break;
    case RNAe: for (R_xlen_t i = 0; i < m_motif.ncol(); ++i) {
              Rcpp::NumericVector tmp = m_motif(Rcpp::_, i);
              consensus_tmp[i] = get_consensusC(Rcpp::as<std::vector<double>>(tmp),
                 "RNA", Rcpp::as<std::string>(type[0]), pseudocount);
            }
            Rcpp::colnames(m_motif) = consensus_tmp;
            consensus = Rcpp::collapse(consensus_tmp);
            break;
    case AAe: for (R_xlen_t i =0; i < m_motif.ncol(); ++i) {
              Rcpp::NumericVector tmp = m_motif(Rcpp::_, i);
              consensus_tmp[i] = get_consensusAAC(Rcpp::as<std::vector<double>>(tmp),
                  Rcpp::as<std::string>(type[0]), pseudocount);
            }
            Rcpp::colnames(m_motif) = consensus_tmp;
            consensus = Rcpp::collapse(consensus_tmp);
            break;
    default: {
               Rcpp::StringVector mot_rownames = Rcpp::rownames(m_motif);
               Rcpp::colnames(m_motif) = Rcpp::StringVector::create();
               Rcpp::rownames(m_motif) = mot_rownames;
               consensus = Rcpp::StringVector::create();
               break;
             }
  }

  return consensus;

}

Rcpp::StringVector check_length(const Rcpp::StringVector &m_name,
    const Rcpp::StringVector &m_altname, const Rcpp::StringVector &m_family,
    const Rcpp::StringVector &m_organism, const Rcpp::StringVector &m_alphabet,
    const Rcpp::StringVector &m_type, const Rcpp::NumericVector &m_icscore,
    const Rcpp::NumericVector &m_nsites, const Rcpp::NumericVector &m_pseudocount,
    const Rcpp::NumericVector &m_bkgsites, const Rcpp::StringVector &m_consensus,
    const Rcpp::StringVector &m_strand, const Rcpp::NumericVector &m_pval,
    const Rcpp::NumericVector &m_qval, const Rcpp::NumericVector &m_eval,
    Rcpp::StringVector msg) {

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

bool check_bkg_names(const Rcpp::StringVector &alph, const std::string &blet) {

  Rcpp::LogicalVector failed(blet.size(), true);

  for (size_t i = 0; i < blet.length(); ++i) {

    for (R_xlen_t j = 0; j < alph.size(); ++j) {

      std::string alph_j = Rcpp::as<std::string>(alph[j]);

      if (alph_j[0] == blet[i]) {
        failed[i] = false;
        break;
      }

    }

  }

  return Rcpp::is_true(Rcpp::any(failed));

}

Rcpp::StringVector check_bkg(const Rcpp::NumericVector &bkg,
    const Rcpp::StringVector &alph, Rcpp::StringVector msg) {

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

    Rcpp::StringVector bnames = bkg.names();
    bool zero_check = false;
    Rcpp::LogicalVector name_comp;
    Rcpp::StringVector alph_i;
    for (R_xlen_t i = 0; i < alen; ++i) {
      alph_i = Rcpp::rep(Rcpp::StringVector::create(alph[i]), blen);
      name_comp = bnames == alph_i;
      if (Rcpp::is_false(Rcpp::any(name_comp))) zero_check = true;
    }

    if (zero_check)

      msg.push_back("* bkg must contain 0-order possibilities for all letters");

    else {

      Rcpp::LogicalVector low_check = bkg < 0;
      Rcpp::LogicalVector high_check = bkg > 1;

      if (Rcpp::is_true(Rcpp::any(low_check)))
        msg.push_back("* bkg does not allow values less than 0");

      if (Rcpp::is_true(Rcpp::any(high_check)))
        msg.push_back("* bkg does not allow values higher than 1");

      Rcpp::NumericVector bkg_zero = bkg[alph];
      double bkg_sum = Rcpp::sum(bkg_zero);
      if (bkg_sum < 0.99 || bkg_sum > 1.01)
        msg.push_back("* 0-order bkg probabilities must add up to 1");

    }

    for (R_xlen_t i = 0; i < blen; ++i) {

      if (check_bkg_names(alph, Rcpp::as<std::string>(bnames[i]))) {
        msg.push_back("* unknown letters found in bkg names");
        break;
      }

    }

  }

  return msg;

}

Rcpp::StringVector check_char_slots(const Rcpp::StringVector &m_type,
    const Rcpp::StringVector &m_strand, Rcpp::StringVector msg) {

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

Rcpp::StringVector check_motif_and_type(const Rcpp::NumericMatrix &m_motif,
    const Rcpp::StringVector &m_type, const Rcpp::NumericVector &m_nsites,
    Rcpp::StringVector msg) {

  SEXP s_type = m_type[0];
  if (Rf_isNull(s_type)) return msg;

  Rcpp::NumericVector motif_colsums = Rcpp::colSums(m_motif);

  switch (::TYPES2_e[m_type[0]]) {
    case PCMe: {
              if (m_nsites.size() > 0) {
                Rcpp::NumericVector colsums_unique = Rcpp::unique(m_nsites);
                if (colsums_unique.size() > 1)
                  msg.push_back("* for type PCM motif colSums must equal nsites");
              }
              Rcpp::LogicalVector int_check = (m_motif < 1.0) & (m_motif > 0.0);
              if (Rcpp::is_true(any(int_check)))
                msg.push_back("* type PCM motifs cannot contain values between 0 and 1");
              break;
            }
    case PPMe: {
              Rcpp::LogicalVector colsums_1_check = (motif_colsums > 0.99)
                & (motif_colsums < 1.01);
              if (Rcpp::is_false(Rcpp::all(colsums_1_check)))
                msg.push_back("* for type PPM colSums must equal 1");
              Rcpp::LogicalVector motif_pos_check = m_motif >= 0;
              if (Rcpp::is_false(Rcpp::all(motif_pos_check)))
                msg.push_back("* for type PPM only positive values are allowed");
              break;
            }
    case ICMe: if (Rcpp::is_true(Rcpp::any(m_motif < 0)))
              msg.push_back("* type ICM motifs cannot contain negative values");
            break;
  }

  return msg;

}

Rcpp::StringVector check_alphabet(const Rcpp::NumericMatrix &m_motif,
    const Rcpp::StringVector &m_alphabet, Rcpp::StringVector msg) {

  Rcpp::StringVector m_rownames = Rcpp::rownames(m_motif);

  switch (::ALPHS_e[m_alphabet[0]]) {
    case DNAe: {
              if (m_motif.nrow() != 4)
                msg.push_back("* DNA/RNA motifs must have 4 rows");
              Rcpp::LogicalVector rownames_check = m_rownames == ::DNA;
              if (Rcpp::is_false(Rcpp::all(rownames_check)))
                msg.push_back("* rownames must be A, C, G, T");
              break;
            }
    case RNAe: {
              if (m_motif.nrow() != 4)
                msg.push_back("* DNA/RNA motifs must have 4 rows");
              Rcpp::LogicalVector rownames_check = m_rownames == ::RNA;
              if (Rcpp::is_false(Rcpp::all(rownames_check)))
                msg.push_back("* rownames must be A, C, G, U");
              break;
            }
    case AAe: {
              if (m_motif.nrow() != 20)
                msg.push_back("* AA motifs must have 20 rows");
              Rcpp::LogicalVector rownames_check = m_rownames == ::AMINOACIDS;
              if (Rcpp::is_false(Rcpp::all(rownames_check)))
                msg.push_back("* rownames must be ACDEFGHIKLMNPQRSTVWY");
              break;
            }
    default: {
               if (m_alphabet[0] == "custom") break;
               if (m_motif.nrow() != m_alphabet[0].size())
                 msg.push_back("* alphabet length does not match number of rows in motif");
               Rcpp::StringVector alph_split;
               for (R_xlen_t i = 0; i < m_alphabet[0].size(); ++i) {
                 alph_split.push_back(m_alphabet[0][i]);
               }
               Rcpp::LogicalVector rownames_check = Rcpp::sort_unique(alph_split) == m_rownames;
               if (Rcpp::is_false(Rcpp::all(rownames_check)))
                 msg.push_back("* rownames must match alphabet and be in alphabetical order");
               break;
            }
  }

  return msg;

}

Rcpp::StringVector check_consensus(const Rcpp::StringVector &m_consensus,
    const Rcpp::NumericMatrix &m_motif, Rcpp::StringVector msg) {

  if (m_consensus.size() > 0) {

    if (m_consensus[0].size() != m_motif.ncol()) {

      msg.push_back("* consensus string must have the same number of letters as motif positions");

    } else {

      Rcpp::StringVector consensus_split;
      Rcpp::StringVector motif_colnames = Rcpp::colnames(m_motif);

      for (R_xlen_t i = 0; i < m_consensus[0].size(); ++i) {
        consensus_split.push_back(m_consensus[0][i]);
        if (consensus_split[i] != motif_colnames[i])
          msg.push_back("* consensus string must match colnames");
      }

    }

  }

  return msg;

}

Rcpp::StringVector check_gap(const Rcpp::RObject &gap, const R_xlen_t ncol,
    Rcpp::StringVector msg) {

  Rcpp::LogicalVector isgapped = gap.slot("isgapped");
  if (isgapped.length() != 1)
    msg.push_back("* isgapped must be a length one logical vector");

  Rcpp::NumericVector gaploc = gap.slot("gaploc");
  Rcpp::NumericVector mingap = gap.slot("mingap");
  Rcpp::NumericVector maxgap = gap.slot("maxgap");

  if (gaploc.length() != mingap.length() || gaploc.length() != maxgap.length())
    msg.push_back("* gaploc, mingap and maxgap should all be the same length");

  if (gaploc.length() > 1) {
    for (R_xlen_t i = 0; i < gaploc.length(); ++i) {
      if (gaploc[i] <= 0)
        msg.push_back("* position 0 gaps or less are not allowed");
      if (gaploc[i] >= ncol) {
        msg.push_back("* gap location values should not exceed motif size - 1");
      }
    }
  }

  return msg;

}

/* C++ ENTRY ---------------------------------------------------------------- */

// [[Rcpp::export(rng = false)]]
Rcpp::S4 universalmotif_cpp(
    Rcpp::NumericMatrix motif,
    Rcpp::String name = "new motif",
    Rcpp::StringVector altname = NA_STRING,
    Rcpp::StringVector family = NA_STRING,
    Rcpp::StringVector organism = NA_STRING,
    Rcpp::StringVector alphabet = "DNA",
    Rcpp::StringVector type = NA_STRING,
    Rcpp::NumericVector icscore = Rcpp::NumericVector::create(),
    Rcpp::NumericVector nsites = Rcpp::NumericVector::create(),
    double pseudocount = 1.0,
    Rcpp::NumericVector bkg = Rcpp::NumericVector::create(),
    Rcpp::NumericVector bkgsites = Rcpp::NumericVector::create(),
    Rcpp::StringVector consensus = NA_STRING,
    Rcpp::String strand = "+-",
    Rcpp::NumericVector pval = Rcpp::NumericVector::create(),
    Rcpp::NumericVector qval = Rcpp::NumericVector::create(),
    Rcpp::NumericVector eval = Rcpp::NumericVector::create(),
    Rcpp::StringVector extrainfo = NA_STRING,
    Rcpp::LogicalVector isgapped = NA_LOGICAL,
    Rcpp::NumericVector gaploc = Rcpp::NumericVector::create(),
    Rcpp::NumericVector mingap = Rcpp::NumericVector::create(),
    Rcpp::NumericVector maxgap = Rcpp::NumericVector::create()) {

  Rcpp::S4 x("universalmotif");

  Rcpp::NumericMatrix m_motif = Rcpp::clone(motif);
  Rcpp::NumericVector motif_colsums = Rcpp::colSums(m_motif);
  // sometimes nsites can slip through (?) as nan (not R_NaN)
  if (std::isnan(nsites[0])) nsites = Rcpp::NumericVector::create();

  if (Rcpp::StringVector::is_na(Rcpp::StringVector::create(name)[0]))
    name = "[no name]";
  // name
  x.slot("name") = name;

  // altname
  if (!Rcpp::StringVector::is_na(altname[0]) && altname.length() > 0)
    x.slot("altname") = altname;

  // family
  if (!Rcpp::StringVector::is_na(family[0]) && family.length() > 0)
    x.slot("family") = family;

  // organism
  if (!Rcpp::StringVector::is_na(organism[0]) && organism.length() > 0)
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
  if (!Rcpp::NumericVector::is_na(bkgsites[0]))
    x.slot("bkgsites") = bkgsites;

  // consensus (&m_motif)
  consensus = universalmotif_consensus(m_motif, alphabet, type, pseudocount, consensus);
  x.slot("consensus") = consensus;

  // motif
  x.slot("motif") = m_motif;

  // strand
  if (strand == "-+") strand = "+-";
  x.slot("strand") = strand;

  // pval
  if (!Rcpp::NumericVector::is_na(pval[0]))
    x.slot("pval") = pval;

  // qval
  if (!Rcpp::NumericVector::is_na(qval[0]))
    x.slot("qval") = qval;

  // eval
  if (!Rcpp::NumericVector::is_na(eval[0]))
    x.slot("eval") = eval;

  // extrainfo
  if (!Rcpp::StringVector::is_na(extrainfo[0]))
    x.slot("extrainfo") = extrainfo;

  // gapinfo
  Rcpp::S4 gap("universalmotif_gapped");

  if (!Rcpp::LogicalVector::is_na(isgapped[0]) && isgapped.length() == 1)
    gap.slot("isgapped") = isgapped;
  else
    gap.slot("isgapped") = false;

  if (!Rcpp::NumericVector::is_na(gaploc[0]) && gaploc.length() > 0)
    gap.slot("gaploc") = mingap;

  if (!Rcpp::NumericVector::is_na(mingap[0]) && mingap.length() > 0)
    gap.slot("mingap") = mingap;

  if (!Rcpp::NumericVector::is_na(maxgap[0]) && maxgap.length() > 0)
    gap.slot("maxgap") = maxgap;

  x.slot("gapinfo") = gap;

  return x;

}

// [[Rcpp::export(rng = false)]]
Rcpp::StringVector validObject_universalmotif(const Rcpp::S4 &motif,
    const bool throw_error = true) {

  Rcpp::StringVector msg;

  Rcpp::StringVector  m_name        = motif.slot("name");
  Rcpp::StringVector  m_altname     = motif.slot("altname");
  Rcpp::StringVector  m_family      = motif.slot("family");
  Rcpp::StringVector  m_organism    = motif.slot("organism");
  Rcpp::NumericMatrix m_motif       = motif.slot("motif");
  Rcpp::StringVector  m_alphabet    = motif.slot("alphabet");
  Rcpp::StringVector  m_type        = motif.slot("type");
  Rcpp::NumericVector m_icscore     = motif.slot("icscore");
  Rcpp::NumericVector m_nsites      = motif.slot("nsites");
  Rcpp::NumericVector m_pseudocount = motif.slot("pseudocount");
  Rcpp::NumericVector m_bkg         = motif.slot("bkg");
  Rcpp::NumericVector m_bkgsites    = motif.slot("bkgsites");
  Rcpp::StringVector  m_consensus   = motif.slot("consensus");
  Rcpp::StringVector  m_strand      = motif.slot("strand");
  Rcpp::NumericVector m_pval        = motif.slot("pval");
  Rcpp::NumericVector m_qval        = motif.slot("qval");
  Rcpp::NumericVector m_eval        = motif.slot("eval");
  Rcpp::RObject       m_gap         = motif.slot("gapinfo");

  // check gap
  msg = check_gap(m_gap, m_motif.ncol(), msg);

  // slot length checks
  msg = check_length(m_name, m_altname, m_family, m_organism, m_alphabet, m_type,
      m_icscore, m_nsites, m_pseudocount, m_bkgsites, m_consensus, m_strand,
      m_pval, m_qval, m_eval, msg);

  // character slot checks
  msg = check_char_slots(m_type, m_strand, msg);

  // bkg slot check
  msg = check_bkg(m_bkg, Rcpp::rownames(m_motif), msg);

  // check motif and type
  msg = check_motif_and_type(m_motif, m_type, m_nsites, msg);

  // check motif matches alphabet
  msg = check_alphabet(m_motif, m_alphabet, msg);

  // consensus check
  msg = check_consensus(m_consensus, m_motif, msg);

  if (msg.length() > 0) {

    msg = Rcpp::StringVector::create(all_checks_collapse(msg));

    if (throw_error) Rcpp::stop(Rcpp::as<std::string>(msg[0]));

  }

  return msg;

}

// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame summarise_motifs_cpp(const Rcpp::List &motifs) {

  R_xlen_t n = motifs.size();
  Rcpp::StringVector name(n);
  Rcpp::StringVector altname(n);
  Rcpp::StringVector family(n);
  Rcpp::StringVector organism(n);
  Rcpp::StringVector alphabet(n);
  Rcpp::NumericVector icscore(n);
  Rcpp::IntegerVector nsites(n);
  Rcpp::IntegerVector bkgsites(n);
  Rcpp::StringVector consensus(n);
  Rcpp::StringVector strand(n);
  Rcpp::NumericVector pval(n);
  Rcpp::NumericVector qval(n);
  Rcpp::NumericVector eval(n);

  for (R_xlen_t i = 0; i < n; ++i) {

    Rcpp::S4 mot = motifs[i];

    Rcpp::StringVector name_i = mot.slot("name");
    name[i] = name_i[0];

    Rcpp::StringVector altname_i = mot.slot("altname");
    altname[i] = altname_i.size() == 1 ? altname_i[0] : NA_STRING;

    Rcpp::StringVector family_i = mot.slot("family");
    family[i] = family_i.size() == 1 ? family_i[0] : NA_STRING;

    Rcpp::StringVector organism_i = mot.slot("organism");
    organism[i] = organism_i.size() == 1 ? organism_i[0] : NA_STRING;

    Rcpp::StringVector alphabet_i = mot.slot("alphabet");
    alphabet[i] = alphabet_i[0];

    Rcpp::NumericVector icscore_i = mot.slot("icscore");
    icscore[i] = icscore_i[0];

    Rcpp::IntegerVector nsites_i = mot.slot("nsites");
    nsites[i] = nsites_i.size() == 1 ? nsites_i[0] : NA_INTEGER;

    Rcpp::IntegerVector bkgsites_i = mot.slot("bkgsites");
    bkgsites[i] = bkgsites_i.size() == 1 ? bkgsites_i[0] : NA_INTEGER;

    Rcpp::StringVector consensus_i = mot.slot("consensus");
    consensus[i] = consensus_i.size() == 1 ? consensus_i[0] : NA_STRING;

    Rcpp::StringVector strand_i = mot.slot("strand");
    strand[i] = strand_i[0];

    Rcpp::NumericVector pval_i = mot.slot("pval");
    pval[i] = pval_i.size() == 1 ? pval_i[0] : NA_REAL;

    Rcpp::NumericVector qval_i = mot.slot("qval");
    qval[i] = qval_i.size() == 1 ? qval_i[0] : NA_REAL;

    Rcpp::NumericVector eval_i = mot.slot("eval");
    eval[i] = eval_i.size() == 1 ? eval_i[0] : NA_REAL;

  }

  return Rcpp::DataFrame::create(

      Rcpp::_["name"]      = name,
      Rcpp::_["altname"]   = altname,
      Rcpp::_["family"]    = family,
      Rcpp::_["organism"]  = organism,
      Rcpp::_["alphabet"]  = alphabet,
      Rcpp::_["icscore"]   = icscore,
      Rcpp::_["nsites"]    = nsites,
      Rcpp::_["bkgsites"]  = bkgsites,
      Rcpp::_["consensus"] = consensus,
      Rcpp::_["strand"]    = strand,
      Rcpp::_["pval"]      = pval,
      Rcpp::_["qval"]      = qval,
      Rcpp::_["eval"]      = eval,

      Rcpp::_["stringsAsFactors"] = false);

}

// [[Rcpp::export(rng = false)]]
Rcpp::List universalmotif_to_list(const Rcpp::S4 &motif) {

  return Rcpp::List::create(

      Rcpp::_["name"]        = motif.slot("name"),
      Rcpp::_["altname"]     = motif.slot("altname"),
      Rcpp::_["family"]      = motif.slot("family"),
      Rcpp::_["organism"]    = motif.slot("organism"),
      Rcpp::_["motif"]       = motif.slot("motif"),
      Rcpp::_["alphabet"]    = motif.slot("alphabet"),
      Rcpp::_["type"]        = motif.slot("type"),
      Rcpp::_["icscore"]     = motif.slot("icscore"),
      Rcpp::_["nsites"]      = motif.slot("nsites"),
      Rcpp::_["pseudocount"] = motif.slot("pseudocount"),
      Rcpp::_["bkg"]         = motif.slot("bkg"),
      Rcpp::_["bkgsites"]    = motif.slot("bkgsites"),
      Rcpp::_["consensus"]   = motif.slot("consensus"),
      Rcpp::_["strand"]      = motif.slot("strand"),
      Rcpp::_["pval"]        = motif.slot("pval"),
      Rcpp::_["qval"]        = motif.slot ("qval"),
      Rcpp::_["eval"]        = motif.slot("eval"),
      Rcpp::_["multifreq"]   = motif.slot("multifreq"),
      Rcpp::_["extrainfo"]   = motif.slot("extrainfo")

      );

}
