#include <Rcpp.h>
#include <RcppThread.h>
#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include "types.h"
#include "score-grid.h"

/* Scanning with gapped motifs.
 *
 * A gapped motif is a series of pieces with a gap of variable size between
 * each pair of neighbouring pieces, so every combination of gap sizes amounts
 * to a separate motif with blank columns in its gaps. Rather than scanning
 * each of those, every piece is scored once per position and the combinations
 * are built by adding up the scores of the pieces, dropping a partial
 * combination as soon as the pieces still to come can no longer reach the
 * threshold. Blank columns score zero, so the scores are those of scanning
 * each combination in full. Letters outside the alphabet score -999999 each
 * wherever they fall in a window, gaps included, as they do when scanning
 * ungapped motifs.
 *
 * Each motif is given in the orientation it is scanned in, which for the
 * reverse strand means the reverse complement of the scoring matrix and the
 * pieces and gaps in reverse order. Combinations of gap sizes are numbered
 * in the motif's own orientation with the first gap varying fastest, and
 * `strides` holds the step in that number for each gap, in scanning order.
 */

namespace {

const int NA_SCORE = -999999;
const std::size_t GAP_CHUNK = 1 << 16;

typedef std::array<int, 256> gap_lookup_t;

struct GapLayout {
  list_mat_t pieces;      // pieces[j][column][letter]
  vec_int_t widths;       // width of each piece
  vec_int_t lo, hi;       // smallest and largest size of each gap
  vec_int_t strides;      // combination step of each gap
  std::vector<std::int64_t> best_rest;  // best total score of pieces j..n-1
  vec_int_t min_rest;     // narrowest span from the start of piece j on
  int max_span;
  int min_score;
};

struct GapHit {
  int motif;
  int sequence;
  int start;  // 0-based
  int width;
  int score;
  int combo;
};

/* The search for the combinations of one motif in one chunk of one
 * sequence. Piece scores are held for the positions from `offset` on. */
struct GapFinder {
  const GapLayout &lay;
  const list_int_t &ps;
  const vec_int_t &na_prefix;  // empty if the sequence has no unknown letters
  std::size_t offset;
  std::size_t seqlen;
  int motif;
  int sequence;
  std::vector<GapHit> &out;

  void extend(std::size_t start, std::size_t j, std::size_t pos,
      std::int64_t partial, int combo) {
    const std::size_t n = lay.widths.size();
    if (j + 1 == n) {
      if (partial >= lay.min_score) {
        int score = partial < INT_MIN ? INT_MIN : static_cast<int>(partial);
        out.push_back({motif, sequence, static_cast<int>(start),
                       static_cast<int>(pos + lay.widths[j] - start),
                       score, combo});
      }
      return;
    }
    if (partial + lay.best_rest[j + 1] < lay.min_score) return;
    std::size_t gap_start = pos + lay.widths[j];
    for (int g = lay.lo[j]; g <= lay.hi[j]; ++g) {
      std::size_t next = gap_start + g;
      if (next + lay.min_rest[j + 1] > seqlen) break;
      std::int64_t s = partial + ps[j + 1][next - offset];
      if (!na_prefix.empty()) {
        int nas = na_prefix[next - offset] - na_prefix[gap_start - offset];
        s += static_cast<std::int64_t>(nas) * NA_SCORE;
      }
      extend(start, j + 1, next, s, combo + (g - lay.lo[j]) * lay.strides[j]);
    }
  }
};

void scan_gapped_seq(const GapLayout &lay, const std::string &seq,
    const gap_lookup_t &lookup, bool has_na, int motif, int sequence,
    std::vector<GapHit> &out) {

  const std::size_t L = seq.size();
  const std::size_t n = lay.widths.size();
  const std::size_t narrowest = lay.min_rest[0];
  if (L < narrowest) return;
  const std::size_t last_start = L - narrowest;  // inclusive

  vec_int_t codes;
  vec_int_t na_prefix;
  list_int_t ps(n);

  for (std::size_t c = 0; c <= last_start; c += GAP_CHUNK) {
    std::size_t e = std::min(c + GAP_CHUNK, last_start + 1);  // starts [c, e)
    std::size_t buf_end = std::min(L, e - 1 + lay.max_span);  // letters [c, buf_end)
    std::size_t nbuf = buf_end - c;

    codes.resize(nbuf);
    for (std::size_t p = 0; p < nbuf; ++p) {
      codes[p] = lookup[static_cast<unsigned char>(seq[c + p])];
    }
    if (has_na) {
      na_prefix.assign(nbuf + 1, 0);
      for (std::size_t p = 0; p < nbuf; ++p) {
        na_prefix[p + 1] = na_prefix[p] + (codes[p] < 0 ? 1 : 0);
      }
    } else {
      na_prefix.clear();
    }

    // Each piece's score wherever it fits in the chunk's letters.
    for (std::size_t j = 0; j < n; ++j) {
      const std::size_t w = lay.widths[j];
      const list_int_t &piece = lay.pieces[j];
      std::size_t npos = nbuf >= w ? nbuf - w + 1 : 0;
      ps[j].resize(npos);
      for (std::size_t p = 0; p < npos; ++p) {
        int tmp = 0;
        for (std::size_t k = 0; k < w; ++k) {
          int idx = codes[p + k];
          tmp += idx < 0 ? NA_SCORE : piece[k][idx];
        }
        ps[j][p] = tmp;
      }
    }

    GapFinder find{lay, ps, na_prefix, c, L, motif, sequence, out};
    for (std::size_t s = c; s < e; ++s) {
      std::int64_t first = ps[0][s - c];
      find.extend(s, 0, s, first, 0);
    }
  }

}

}  // namespace

// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame scan_gapped_cpp(const Rcpp::List &score_mats,
    const std::vector<std::string> &seq_vecs, const std::string &alph,
    const std::vector<double> &min_scores, const Rcpp::List &piece_widths,
    const Rcpp::List &gap_min, const Rcpp::List &gap_max,
    const Rcpp::List &gap_strides, const int &nthreads,
    const bool &warnNA = true) {

  gap_lookup_t lookup;
  lookup.fill(-1);
  for (std::size_t a = 0; a < alph.size(); ++a) {
    lookup[static_cast<unsigned char>(alph[a])] = static_cast<int>(a);
  }

  std::vector<char> has_na(seq_vecs.size(), 0);
  for (std::size_t i = 0; i < seq_vecs.size(); ++i) {
    for (std::size_t j = 0; j < seq_vecs[i].size(); ++j) {
      if (lookup[static_cast<unsigned char>(seq_vecs[i][j])] < 0) {
        has_na[i] = 1;
        break;
      }
    }
  }
  if (warnNA && std::find(has_na.begin(), has_na.end(), 1) != has_na.end()) {
    Rcpp::warning("Non-standard letters detected. These were ignored.");
  }

  const std::size_t nmot = score_mats.size();
  std::vector<GapLayout> layouts(nmot);
  for (std::size_t i = 0; i < nmot; ++i) {
    GapLayout &lay = layouts[i];
    Rcpp::NumericMatrix mat = score_mats[i];
    lay.widths = Rcpp::as<vec_int_t>(piece_widths[i]);
    lay.lo = Rcpp::as<vec_int_t>(gap_min[i]);
    lay.hi = Rcpp::as<vec_int_t>(gap_max[i]);
    lay.strides = Rcpp::as<vec_int_t>(gap_strides[i]);
    const std::size_t n = lay.widths.size();
    if (n == 0 || lay.lo.size() + 1 != n || lay.hi.size() + 1 != n ||
        lay.strides.size() + 1 != n)
      Rcpp::stop("malformed gapped motif layout");

    double val = score_grid_ceiling(min_scores[i]);
    if (val >= static_cast<double>(INT_MAX)) lay.min_score = INT_MAX;
    else if (val <= static_cast<double>(INT_MIN)) lay.min_score = INT_MIN;
    else lay.min_score = static_cast<int>(val);

    // Quantise as the ungapped scanner does (* 1000, truncated toward zero),
    // so that adding up the pieces gives the same integer scores.
    R_xlen_t col = 0;
    std::vector<std::int64_t> best(n, 0);
    lay.pieces.resize(n);
    for (std::size_t j = 0; j < n; ++j) {
      lay.pieces[j].reserve(lay.widths[j]);
      for (int k = 0; k < lay.widths[j]; ++k, ++col) {
        vec_int_t column(mat.nrow());
        int colmax = INT_MIN;
        for (R_xlen_t r = 0; r < mat.nrow(); ++r) {
          column[r] = static_cast<int>(mat(r, col) * 1000.0);
          colmax = std::max(colmax, column[r]);
        }
        best[j] += colmax;
        lay.pieces[j].push_back(std::move(column));
      }
    }
    if (col != mat.ncol()) Rcpp::stop("gapped motif pieces do not match its width");

    lay.best_rest.assign(n + 1, 0);
    lay.min_rest.assign(n + 1, 0);
    for (std::size_t j = n; j-- > 0; ) {
      lay.best_rest[j] = lay.best_rest[j + 1] + best[j];
      lay.min_rest[j] = lay.min_rest[j + 1] + lay.widths[j] +
        (j + 1 < n ? lay.lo[j] : 0);
    }
    lay.max_span = 0;
    for (std::size_t j = 0; j < n; ++j) {
      lay.max_span += lay.widths[j] + (j + 1 < n ? lay.hi[j] : 0);
    }
  }

  std::vector<std::vector<GapHit>> hit_bufs(nmot);
  RcppThread::parallelFor(0, nmot,
      [&hit_bufs, &layouts, &seq_vecs, &lookup, &has_na] (std::size_t i) {
        if (layouts[i].min_score == INT_MAX) return;
        for (std::size_t j = 0; j < seq_vecs.size(); ++j) {
          scan_gapped_seq(layouts[i], seq_vecs[j], lookup, has_na[j],
              int(i) + 1, int(j) + 1, hit_bufs[i]);
        }
      }, nthreads);

  std::size_t total = 0;
  for (auto &b : hit_bufs) total += b.size();
  Rcpp::IntegerVector motif(total), sequence(total), start(total), stop(total),
    combo(total);
  Rcpp::NumericVector score(total);
  Rcpp::CharacterVector match(total);
  std::size_t h = 0;
  std::string m;
  for (std::size_t i = 0; i < nmot; ++i) {
    const GapLayout &lay = layouts[i];
    for (auto &hit : hit_bufs[i]) {
      motif[h] = hit.motif;
      sequence[h] = hit.sequence;
      start[h] = hit.start + 1;
      stop[h] = hit.start + hit.width;
      score[h] = static_cast<double>(hit.score) / 1000;
      combo[h] = hit.combo;
      // The window's letters, with the gaps of this combination as dots.
      m = seq_vecs[hit.sequence - 1].substr(hit.start, hit.width);
      std::size_t at = 0;
      for (std::size_t j = 0; j + 1 < lay.widths.size(); ++j) {
        at += lay.widths[j];
        int size = lay.hi[j] - lay.lo[j] + 1;
        int g = lay.lo[j] + (hit.combo / lay.strides[j]) % size;
        for (int k = 0; k < g; ++k) m[at + k] = '.';
        at += g;
      }
      match[h] = m;
      ++h;
    }
  }

  return Rcpp::DataFrame::create(
        Rcpp::_["motif"] = motif,
        Rcpp::_["motif.i"] = motif,
        Rcpp::_["sequence"] = sequence,
        Rcpp::_["sequence.i"] = sequence,
        Rcpp::_["start"] = start,
        Rcpp::_["stop"] = stop,
        Rcpp::_["score"] = score,
        Rcpp::_["match"] = match,
        Rcpp::_["combo"] = combo,
        Rcpp::_["stringsAsFactors"] = false
      );

}
