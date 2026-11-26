#' Name motifs by the database motifs they resemble most.
#'
#' `name_motifs()` finds the motifs of a database, such as JASPAR, that each
#' motif resembles most, to name the factors behind motifs found by
#' discovery.
#'
#' Each motif is compared with every database motif by
#' [compare_motifs_lite()] with `complete = TRUE`, which scores the whole of
#' the query as TOMTOM does and picks the same best match as TOMTOM for most
#' motifs. Matches are ranked by P-value. Their q-values are
#' Benjamini-Hochberg adjusted for the number of database motifs each motif
#' was compared with, so that a motif's matches mean the same whatever other
#' motifs are named alongside it. Up to `top` matches with a q-value of at
#' most `qvalue` are kept for each motif. Gapped motifs are compared with
#' their gaps, as by [compare_motifs_lite()].
#'
#' @param motifs See [convert_motifs()] for accepted motif formats. The
#'   motifs to name, DNA or RNA only.
#' @param database The motifs to name them by, such as those of
#'   [read_jaspar()] or [read_meme()] for a JASPAR file. Their `altname`,
#'   usually the factor's name, is given beside their `name`.
#' @param top `integer(1)`. The most matches to keep for each motif.
#'   Default `3`.
#' @param qvalue `numeric(1)`. The largest q-value of a match kept. Default
#'   `0.05`.
#' @param min.overlap `integer(1)`. The fewest columns two motifs must
#'   overlap to be compared, lowered to the width of the narrower motif when
#'   needed. Default `5`.
#' @param nthreads `numeric(1)`. Number of threads. `nthreads = 0` uses all
#'   available threads.
#'
#' @return A `data.frame` with one row per match, ordered by motif and then
#'   rank, with columns `motif` and `motif.i` (the motif's name and position
#'   in `motifs`), `rank` (counting from 1), `match` and `altname` (the
#'   database motif's names), `score`, `pvalue` and `qvalue`, `strand` and
#'   `offset` (as [compare_motifs_lite()] gives them, with the motif as the
#'   query) and `consensus` and `match.consensus` (the consensus sequences of
#'   the aligned columns). Motifs without a match passing `qvalue` have no
#'   rows.
#'
#' @references
#'
#' Benjamini Y, Hochberg Y (1995). "Controlling the false discovery rate: a
#' practical and powerful approach to multiple testing." *Journal of the
#' Royal Statistical Society: Series B*, **57**(1), 289-300.
#'
#' Gupta S, Stamatoyannopoulos JA, Bailey TL, Noble WS (2007). "Quantifying
#' similarity between motifs." *Genome Biology*, **8**(2), R24.
#'
#' @examples
#' db <- list(create_motif("TGACTCA", name = "MA0099", altname = "FOS::JUN"),
#'            create_motif("AGATAA", name = "MA0035", altname = "GATA1"),
#'            create_motif("CACGTG", name = "MA0004", altname = "Arnt"))
#' found <- list(create_motif("ATGAGTCAT", name = "pos_1"),
#'               create_motif("CTTATCT", name = "pos_2"))
#' name_motifs(found, db, top = 1, qvalue = 1)
#'
#' @seealso [compare_motifs_lite()], [compare_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @family lite motif functions
#' @export
name_motifs <- function(motifs, database, top = 3, qvalue = 0.05,
                        min.overlap = 5, nthreads = 1) {

  if (missing(motifs) || missing(database))
    stop("`motifs` and `database` are required", call. = FALSE)
  if (!is.numeric(top) || length(top) != 1L || is.na(top) || top < 1)
    stop("`top` must be a single number of at least 1", call. = FALSE)
  if (!is.numeric(qvalue) || length(qvalue) != 1L || is.na(qvalue) ||
      qvalue <= 0 || qvalue > 1)
    stop("`qvalue` must be a single numeric in (0, 1]", call. = FALSE)

  as_list <- function(x) {
    if (is.list(x) && !length(x)) return(list())
    x <- convert_motifs(x)
    if (!is.list(x)) list(x) else x
  }
  motifs <- as_list(motifs)
  database <- as_list(database)
  if (!length(motifs) || !length(database))
    stop("give at least one motif and one database motif", call. = FALSE)

  ## The motifs are renamed so that names shared between them and the
  ## database cannot get in the way; they are compared by position.
  nq <- length(motifs)
  rename <- function(x, prefix)
    Map(function(m, i) { m@name <- paste0(prefix, i); m }, x, seq_along(x))
  cmp <- compare_motifs_lite(c(rename(motifs, "q"), rename(database, "d")),
                             compare.to = seq_len(nq), qvalue = 1,
                             complete = TRUE, min.overlap = min.overlap,
                             nthreads = nthreads)
  cmp <- cmp[cmp$target.i > nq, , drop = FALSE]
  cmp$db.i <- cmp$target.i - nq

  ## Benjamini-Hochberg over every database motif, within each motif.
  cmp <- cmp[order(cmp$subject.i, cmp$pvalue, -cmp$score, cmp$db.i), ,
             drop = FALSE]
  n.db <- length(database)
  cmp$qvalue <- stats::ave(cmp$pvalue, cmp$subject.i, FUN = function(p) {
    q <- p * n.db / seq_along(p)
    pmin(rev(cummin(rev(q))), 1)
  })
  cmp <- cmp[cmp$qvalue <= qvalue, , drop = FALSE]
  cmp$rank <- stats::ave(cmp$subject.i, cmp$subject.i, FUN = seq_along)
  cmp <- cmp[cmp$rank <= top, , drop = FALSE]

  altname <- vapply(database, function(m) {
    a <- m@altname
    if (length(a) && !is.na(a[1])) a[1] else NA_character_
  }, character(1))
  out <- data.frame(
    motif           = vapply(motifs, function(m) m@name, character(1))[cmp$subject.i],
    motif.i         = as.integer(cmp$subject.i),
    rank            = as.integer(cmp$rank),
    match           = vapply(database, function(m) m@name, character(1))[cmp$db.i],
    altname         = altname[cmp$db.i],
    score           = cmp$score,
    pvalue          = cmp$pvalue,
    qvalue          = cmp$qvalue,
    strand          = cmp$strand,
    offset          = as.integer(cmp$offset),
    consensus       = cmp$subject.consensus,
    match.consensus = cmp$target.consensus,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out

}
