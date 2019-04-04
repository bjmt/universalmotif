#include <Rcpp.h>
using namespace Rcpp;

int score_seq(IntegerVector tmp_seq, IntegerMatrix score_mat) {
  int score = 0;
  for (int i = 0; i < tmp_seq.length(); ++i) {
    score += score_mat(tmp_seq[i], i);
  }
  return score;
}

// [[Rcpp::export(rng = false)]]
IntegerMatrix numeric_to_integer_matrix(NumericMatrix mat) {
  IntegerMatrix out(mat.nrow(), mat.ncol());
  for (int i = 0; i < mat.length(); ++i) {
    out[i] = mat[i] * 1000;
  }
  return out;
}

// [[Rcpp::export(rng = false)]]
IntegerVector scan_seq_internal2(IntegerVector sequence, IntegerMatrix score_mat,
    int min_score) {

  // BUG FIX: Can't deal with NAs generated from non-standard letters.

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

// [[Rcpp::export(rng = false)]]
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

// [[Rcpp::export(rng = false)]]
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

// [[Rcpp::export(rng = false)]]
IntegerVector string_to_factor(StringVector x, StringVector y) {

  StringVector lvls = sort_unique(y);
  IntegerVector out = match(x, lvls);

  out.attr("levels") = as<CharacterVector>(lvls);
  out.attr("class") = "factor";

  return out;

}

// [[Rcpp::export(rng = false)]]
IntegerVector res_to_index(IntegerVector x) {

  int x_len = x.length();

  for (int i = 0; i < x_len; ++i) {
    if (x[i] == 1) x[i] += i;
  }

  LogicalVector y = x != 0;
  x = x[y];

  return x;

}

// [[Rcpp::export(rng = false)]]
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

// [[Rcpp::export(rng = false)]]
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

StringVector create_col_sequence(StringVector seq_names, IntegerVector n_rows,
    int n, int rows_all) {

  int row_offset = 0;
  StringVector out(rows_all);

  for (int i = 0; i < n; ++i) {
    if (n_rows[i] == 0) continue;
    String seq_names_i = seq_names[i];
    for (int j = 0; j < n_rows[i]; ++j) {
      out[j + row_offset] = seq_names_i;
    }
    row_offset += n_rows[i];
  }

  return out;

}

IntegerVector create_col_start(List to_keep, IntegerVector n_rows, int n,
    int rows_all, int k, int mot_lens, String strand) {

  int row_offset = 0;
  IntegerVector out(rows_all);

  if (strand == "+") {

    for (int i = 0; i < n; ++i) {
      if (n_rows[i] == 0) continue;
      IntegerVector to_keep_i = to_keep[i];
      for (int j = 0; j < n_rows[i]; ++j) {
        out[j + row_offset] = to_keep_i[j];
      }
      row_offset += n_rows[i];
    }

  } else if (strand == "-") {

    for (int i = 0; i < n; ++i) {
      if (n_rows[i] == 0) continue;
      IntegerVector to_keep_i = to_keep[i];
      for (int j = 0; j < n_rows[i]; ++j) {
        out[j + row_offset] = to_keep_i[j] + mot_lens + k - 2;
      }
      row_offset += n_rows[i];
    }

  } else {

    stop("strand must be one of +, -");

  }

  return out;

}

IntegerVector create_col_stop(List to_keep, IntegerVector n_rows, int n,
    int rows_all, int k, int mot_lens, String strand) {

  int row_offset = 0;
  IntegerVector out(rows_all);

  if (strand == "+") {

    for (int i = 0; i < n; ++i) {
      if (n_rows[i] == 0) continue;
      IntegerVector to_keep_i = to_keep[i];
      for (int j = 0; j < n_rows[i]; ++j) {
        out[j + row_offset] = to_keep_i[j] + mot_lens + k - 2;
      }
      row_offset += n_rows[i];
    }

  } else if (strand == "-") {

    for (int i = 0; i < n; ++i) {
      if (n_rows[i] == 0) continue;
      IntegerVector to_keep_i = to_keep[i];
      for (int j = 0; j < n_rows[i]; ++j) {
        out[j + row_offset] = to_keep_i[j];
      }
      row_offset += n_rows[i];
    }

  } else {

    stop("strand must be one of +, -");

  }

  return out;

}

NumericVector create_col_score(List to_keep, IntegerVector n_rows, int n,
    int rows_all, List seq_ints, IntegerMatrix score_mats, int mot_lens,
    int k) {

  int row_offset = 0;
  NumericVector out(rows_all);

  for (int i = 0; i < n; ++i) {
    if (n_rows[i] == 0) continue;
    IntegerVector seq_ints_i = seq_ints[i];
    IntegerVector to_keep_i = to_keep[i];
    List hits_i = parse_k_res_helper_1(seq_ints_i, to_keep_i, mot_lens, k);
    for (int j = 0; j < n_rows[i]; ++j) {
      out[j + row_offset] = score_seq(hits_i[j], score_mats);
    }
    for (int j = 0; j < n_rows[i]; ++j) {
      out[j + row_offset] /= 1000.0;
    }
    row_offset += n_rows[i];
  }

  return out;

}

NumericVector create_col_score_pct(NumericVector col_score, double min_scores,
    double max_scores, int rows_all) {

  double min_scores_abs = abs(min_scores);
  double total_score = min_scores_abs + abs(max_scores);

  NumericVector out = 100.0 * (col_score + min_scores_abs) / total_score;

  return out;

}

StringVector create_col_match(List to_keep, List seqs_aschar, int mot_lens,
    IntegerVector n_rows, int n, int rows_all) {

  int row_offset = 0;
  StringVector out(rows_all);

  for (int i = 0; i < n; ++i) {
    if (n_rows[i] == 0) continue;
    IntegerVector to_keep_i = to_keep[i];
    StringVector seqs_aschar_i = seqs_aschar[i];
    List matches_i = parse_k_res_helper_2(seqs_aschar_i, to_keep_i, mot_lens);
    for (int j = 0; j < n_rows[i]; ++j) {
      StringVector matches_i_j = matches_i[j];
      out[j + row_offset] = collapse(matches_i_j);
    }
    row_offset += n_rows[i];
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
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

  StringVector  col_motif(rows_all, mot_names);
  StringVector  col_strand(rows_all, strand);
  NumericVector col_max_score(rows_all, max_scores);

  StringVector  col_sequence  = create_col_sequence(seq_names, n_rows, n,
                                rows_all);
  IntegerVector col_start     = create_col_start(to_keep, n_rows, n, rows_all, k,
                                mot_lens, strand);
  IntegerVector col_stop      = create_col_stop(to_keep, n_rows, n, rows_all, k,
                                mot_lens, strand);
  NumericVector col_score     = create_col_score(to_keep, n_rows, n, rows_all,
                                seq_ints, score_mats, mot_lens, k);
  NumericVector col_score_pct = create_col_score_pct(col_score, min_scores,
                                max_scores, rows_all);
  StringVector  col_match     = create_col_match(to_keep, seqs_aschar, mot_lens,
                                n_rows, n, rows_all);

  List out = List::create(

      _["motif"]     = col_motif,
      _["sequence"]  = col_sequence,
      _["start"]     = col_start,
      _["stop"]      = col_stop,
      _["score"]     = col_score,
      _["max.score"] = col_max_score,
      _["score.pct"] = col_score_pct,
      _["match"]     = col_match,
      _["strand"]    = col_strand

      );

  return out;

}

IntegerVector join_int_vecs(List x, String list_name, IntegerVector n_per,
    int n, int out_len) {

  int row_offset = 0;
  List tmp;
  IntegerVector tmp_vec;
  IntegerVector out(out_len);

  for (int i = 0; i < n; ++i) {
    tmp = x[i];
    tmp_vec = tmp[list_name];
    for (int j = 0; j < n_per[i]; ++j) {
      out[j + row_offset] = tmp_vec[j];
    }
    row_offset += n_per[i];
  }

  return out;

}

NumericVector join_num_vecs(List x, String list_name, IntegerVector n_per,
    int n, int out_len) {

  int row_offset = 0;
  List tmp;
  NumericVector tmp_vec;
  NumericVector out(out_len);

  for (int i = 0; i < n; ++i) {
    tmp = x[i];
    tmp_vec = tmp[list_name];
    for (int j = 0; j < n_per[i]; ++j) {
      out[j + row_offset] = tmp_vec[j];
    }
    row_offset += n_per[i];
  }

  return out;

}

StringVector join_str_vecs(List x, String list_name, IntegerVector n_per,
    int n, int out_len) {

  int row_offset = 0;
  List tmp;
  StringVector tmp_vec;
  StringVector out(out_len);

  for (int i = 0; i < n; ++i) {
    tmp = x[i];
    tmp_vec = tmp[list_name];
    for (int j = 0; j < n_per[i]; ++j) {
      out[j + row_offset] = tmp_vec[j];
    }
    row_offset += n_per[i];
  }

  return out;

}

// [[Rcpp::export(rng = false)]]
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

  StringVector  col_motif     = join_str_vecs(res, "motif", n_per, n, n_all);
  StringVector  col_sequence  = join_str_vecs(res, "sequence", n_per, n, n_all);
  IntegerVector col_start     = join_int_vecs(res, "start", n_per, n, n_all);
  IntegerVector col_stop      = join_int_vecs(res, "stop", n_per, n, n_all);
  NumericVector col_score     = join_num_vecs(res, "score", n_per, n, n_all);
  NumericVector col_max_score = join_num_vecs(res, "max.score", n_per, n, n_all);
  NumericVector col_score_pct = join_num_vecs(res, "score.pct", n_per, n, n_all);
  StringVector  col_match     = join_str_vecs(res, "match", n_per, n, n_all);
  StringVector  col_strand    = join_str_vecs(res, "strand", n_per, n, n_all);

  DataFrame out = DataFrame::create(

      _["motif"]     = col_motif,
      _["sequence"]  = col_sequence,
      _["start"]     = col_start,
      _["stop"]      = col_stop,
      _["score"]     = col_score,
      _["max.score"] = col_max_score,
      _["score.pct"] = col_score_pct,
      _["match"]     = col_match,
      _["strand"]    = col_strand,

      _["stringsAsFactors"] = false);

  return out;

}
