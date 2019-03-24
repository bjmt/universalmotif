#include <Rcpp.h>
using namespace Rcpp;

double score_seq(IntegerVector tmp_seq, NumericMatrix score_mat) {
  double score = 0;
  for (int i = 0; i < tmp_seq.length(); ++i) {
    score += score_mat(tmp_seq[i], i);
  }
  return score;
}

int score_seq_int(IntegerVector tmp_seq, IntegerMatrix score_mat) {
  int score = 0;
  for (int i = 0; i < tmp_seq.length(); ++i) {
    score += score_mat(tmp_seq[i], i);
  }
  return score;
}

// [[Rcpp::export]]
IntegerMatrix numeric_to_integer_matrix(NumericMatrix mat) {
  IntegerMatrix out(mat.nrow(), mat.ncol());
  for (int i = 0; i < mat.length(); ++i) {
    out[i] = mat[i] * 1000;
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector scan_seq_internal2(IntegerVector sequence, IntegerMatrix score_mat,
    int min_score) {

  // BUG FIX: can't deal with NAs generated from non-standard letters

  IntegerVector to_keep(sequence.length());

  int tmp_score;
  int max_step = sequence.size() - score_mat.ncol() + 1;

  for (int i = 0; i < max_step; ++i) {

    tmp_score = 0;
    for (int j = 0; j < score_mat.ncol(); ++j) {
      bool na_check = IntegerVector::is_na(sequence[i + j]);
      if (na_check)
        tmp_score += -999999;
      else
        tmp_score += score_mat(sequence[i + j], j);
    }
    if (tmp_score >= min_score) to_keep[i] = 1;

  }

  return to_keep;

}

// [[Rcpp::export]]
IntegerVector scan_seq_internal(IntegerVector sequence, IntegerMatrix score_mat,
    int min_score) {

  IntegerVector to_keep(sequence.length());

  int tmp_score;
  int max_step = sequence.size() - score_mat.ncol() + 1;

  for (int i = 0; i < max_step; ++i) {

    tmp_score = 0;
    for (int j = 0; j < score_mat.ncol(); ++j) {
      tmp_score += score_mat(sequence[i + j], j);
    }
    if (tmp_score >= min_score) to_keep[i] = 1;

  }

  return to_keep;

}

// [[Rcpp::export]]
IntegerVector LETTER_to_int(IntegerVector seqs, int k, IntegerVector letters) {

  IntegerVector out(seqs.length() / k, 0);
  int out_i, l_;
  int let_length = letters.length();

  if (k == 1) {

    out = seqs;

  } else {

    for (int i = 0; i < seqs.length(); ++i) {
      if (i % k == 0) {

        for (int l = 0; l < k; ++l) {

          l_ = pow(let_length, k - l - 1);
          out_i = seqs[i + l];
          out_i *= l_;
          out[i / k] += out_i;

        }

      }
    }

  }

  return out;

}

// [[Rcpp::export]]
IntegerVector string_to_factor(StringVector x, StringVector y) {

  StringVector lvls = sort_unique(y);
  IntegerVector out = match(x, lvls);

  out.attr("levels") = as<CharacterVector>(lvls);
  out.attr("class") = "factor";

  return out;

}

// [[Rcpp::export]]
IntegerVector res_to_index(IntegerVector x) {

  int x_len = x.length();

  for (int i = 0; i < x_len; ++i) {
    if (x[i] == 1) x[i] += i;
  }

  LogicalVector y = x != 0;
  x = x[y];

  return x;

}

// [[Rcpp::export]]
List parse_k_res_helper_1(IntegerVector seqs, IntegerVector to_keep,
    int mot_len, int k) {

  int n = to_keep.length();
  List out(n);

  for (int i = 0; i < n; ++i) {
    IntegerVector tmp(mot_len - k + 1);
    for (int j = 0; j < mot_len - k + 1; ++j) {
      tmp[j] = seqs[to_keep[i] - 1 + j];
    }
    out[i] = tmp;
  }

  return out;

}

// [[Rcpp::export]]
List parse_k_res_helper_2(StringVector sequence, IntegerVector to_keep,
    int mot_len) {

  int n = to_keep.length();
  List out(n);

  for (int i = 0; i < n; ++i) {
    StringVector tmp(mot_len);
    for (int j = 0; j < mot_len; ++j) {
      tmp[j] = sequence[to_keep[i] - 1 + j];
    }
    out[i] = tmp;
  }

  return out;

}

// [[Rcpp::export]]
List get_res_cpp(List to_keep, List seqs_aschar, List seq_ints,
    int mot_lens, double min_scores, double max_scores, String mot_names,
    StringVector seq_names, IntegerMatrix score_mats, String strand,
    IntegerVector seq_lens, int k) {

  int n = to_keep.length();
  IntegerVector n_rows(n);
  for (int i = 0; i < n; ++i) {
    IntegerVector tmp = to_keep(i);
    n_rows[i] = tmp.length();
  }
  int rows_all = sum(n_rows);

  StringVector col_motif(rows_all, mot_names);
  StringVector col_strand(rows_all, strand);
  NumericVector col_max_score(rows_all, max_scores);

  StringVector col_sequence(rows_all);
  NumericVector col_start(rows_all);
  NumericVector col_stop(rows_all);
  NumericVector col_score(rows_all);
  NumericVector col_score_pct(rows_all);
  StringVector col_match(rows_all);

  int row_offset = 0;
  for (int i = 0; i < n; ++i) {

    if (n_rows[i] == 0) continue;

    IntegerVector to_keep_i = to_keep[i];
    StringVector seqs_aschar_i = seqs_aschar[i];
    IntegerVector seq_ints_i = seq_ints[i];
    String seq_names_i = seq_names[i];
    List hits_i = parse_k_res_helper_1(seq_ints_i, to_keep_i, mot_lens, k);
    List matches_i = parse_k_res_helper_2(seqs_aschar_i, to_keep_i, mot_lens);

    for (int j = 0; j < n_rows[i]; ++j) {

      col_sequence[j + row_offset] = seq_names_i;

      if (strand == "+") {
        col_start[j + row_offset] = to_keep_i[j];
        col_stop[j + row_offset] = col_start[j + row_offset] + mot_lens + k - 2;
      } else if (strand == "-") {
        col_start[j + row_offset] = to_keep_i[j] + mot_lens + k - 2;
        col_stop[j + row_offset] = to_keep_i[j];
      }

      col_score[j + row_offset] = score_seq_int(hits_i[j], score_mats);
      col_score[j + row_offset] /= 1000.0;
      col_score_pct[j + row_offset] = col_score[j + row_offset] /
                                      col_max_score[j + row_offset] * 100;

      StringVector tmp = matches_i[j];
      col_match[j + row_offset] = collapse(tmp);

    }

    row_offset += n_rows[i];

  }

  List out = List::create(_["motif"] = col_motif,
      _["sequence"] = col_sequence, _["start"] = col_start,
      _["stop"] = col_stop, _["score"] = col_score,
      _["max.score"] = col_max_score, _["score.pct"] = col_score_pct,
      _["match"] = col_match, _["strand"] = col_strand,
      _["stringsAsFactors"] = false);

  return out;

}

// [[Rcpp::export]]
DataFrame res_list_to_df_cpp(List res) {

  int n = res.length();
  IntegerVector n_per(n);
  List tmp, tmp2;
  for (int i = 0; i < n; ++i) {
    tmp = res[i];
    tmp2 = tmp["motif"];
    n_per[i] = tmp2.length();
  }

  int n_all = sum(n_per);

  StringVector col_motif(n_all);
  StringVector col_sequence(n_all);
  IntegerVector col_start(n_all);
  IntegerVector col_stop(n_all);
  NumericVector col_score(n_all);
  NumericVector col_max_score(n_all);
  NumericVector col_score_pct(n_all);
  StringVector col_match(n_all);
  StringVector col_strand(n_all);

  int row_offset = 0;

  for (int i = 0; i < n; ++i) {

    List res_ = res[i];

    StringVector col_motif_ = res_["motif"];
    StringVector col_sequence_ = res_["sequence"];
    IntegerVector col_start_ = res_["start"];
    IntegerVector col_stop_ = res_["stop"];
    NumericVector col_score_ = res_["score"];
    NumericVector col_max_score_ = res_["max.score"];
    NumericVector col_score_pct_ = res_["score.pct"];
    StringVector col_match_ = res_["match"];
    StringVector col_strand_ = res_["strand"];

    for (int j = 0; j < n_per(i); ++j) {

      col_motif[j + row_offset] = col_motif_[j];
      col_sequence[j + row_offset] = col_sequence_[j];
      col_start[j + row_offset] = col_start_[j];
      col_stop[j + row_offset] = col_stop_[j];
      col_score[j + row_offset] = col_score_[j];
      col_max_score[j + row_offset] = col_max_score_[j];
      col_score_pct[j + row_offset] = col_score_pct_[j];
      col_match[j + row_offset] = col_match_[j];
      col_strand[j + row_offset] = col_strand_[j];

    }

    row_offset += n_per[i];

  }

  DataFrame out = DataFrame::create(_["motif"] = col_motif,
      _["sequence"] = col_sequence, _["start"] = col_start,
      _["stop"] = col_stop, _["score"] = col_score,
      _["max.score"] = col_max_score, _["score.pct"] = col_score_pct,
      _["match"] = col_match, _["strand"] = col_strand,
      _["stringsAsFactors"] = false);

  return out;

}
