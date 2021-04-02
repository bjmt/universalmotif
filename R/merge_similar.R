#' Identify and merge similar motifs within a collection of motifs.
#'
#' Given a list of motifs, [merge_similar()] will identify similar motifs with
#' [compare_motifs()], and merge similar ones with [merge_motifs()].
#'
#' @param threshold `numeric(1)` The minimum (for similarity metrics) or maximum (for
#'   distance metrics) threshold score for merging.
#' @param threshold.type `character(1)` Type of score used for thresholding.
#'   Currently unused.
#' @param method `character(1)` One of PCC, EUCL, SW, KL, BHAT, HELL,
#'   SEUCL, MAN, WEUCL, WPCC. See [compare_motifs()]. (The ALLR and ALLR_LL
#'   methods cannot be used for distance matrix construction.)
#' @param score.strat.compare `character(1)` The `score.strat` parameter used
#'   by [compare_motifs()]. For clustering purposes, the `"sum"` option cannot
#'   be used.
#' @param score.strat.merge `character(1)` The `score.strat` parameter used
#'   by [merge_motifs()]. As discussed in [merge_motifs()], the `"sum"` option
#'   is recommended over `"a.mean"` to maximize the overlap between motifs.
#'
#' @return See [convert_motifs()] for available output formats.
#'
#' @details
#' See [compare_motifs()] for more info on comparison parameters, and
#' [merge_motifs()] for more info on motif merging.
#'
#' @examples
#' \dontrun{
#' library(MotifDb)
#' motifs <- filter_motifs(MotifDb, family = "bHLH")[1:50]
#' length(motifs)
#' motifs <- merge_similar(motifs)
#' length(motifs)
#' }
#'
#' @seealso [compare_motifs()], [merge_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @inheritParams compare_motifs
#' @export
merge_similar <- function(motifs,
  threshold = 0.95, threshold.type = "score.abs", method = "PCC",
  use.type = "PPM", min.overlap = 6, min.mean.ic = 0,
  tryRC = TRUE, relative_entropy = FALSE, normalise.scores = FALSE,
  min.position.ic = 0, score.strat.compare = "a.mean",
  score.strat.merge = "sum", nthreads = 1) {

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  if (length(motifs) == 1) {
    return(.internal_convert(motifs[[1]], unique(CLASS_IN)))
  }

  # For now, just merge based on score -- add p-value merging later.

  if (score.strat.compare == "sum")
    stop("`score.strat.compare` cannot be \"sum\".", call. = FALSE)
  if (method %in% c("ALLR", "ALLR_LL"))
    stop(wmsg("`method` cannot be \"ALLR\", \"ALLR_LL\" as the resulting similarity",
      " matrix cannot be converted to a distance matrix."), call. = FALSE)
  
  comp.mat <- compare_motifs(motifs, method = method,
    use.type = use.type, min.overlap = min.overlap, min.mean.ic = min.mean.ic,
    tryRC = tryRC, relative_entropy = relative_entropy,
    normalise.scores = normalise.scores, min.position.ic = min.position.ic,
    score.strat = score.strat.compare, nthreads = nthreads)

  if (anyNA(comp.mat))
    stop(wmsg("Clustering is not possible when NA values are present in the ",
        "comparison matrix. Lower `min.mean.ic` and/or `min.position.ic` until ",
        "no more NA values are generated, or omit the problematic motifs."),
      call. = FALSE)

  if (method %in% c("PCC", "SW", "WPCC")) {
    threshold <- comp.mat[1] - threshold
    comp.mat <- comp.mat[1] - comp.mat
  }

  comp.clust <- hclust(as.dist(comp.mat))
  comp.groups <- cutree(comp.clust, h = threshold)

  final.mots <- unname(tapply(motifs, comp.groups,
      function(x) merge_motifs(sort_by_ic(x),
      method = method, use.type = use.type, min.overlap = min.overlap,
      min.mean.ic = min.mean.ic, tryRC = tryRC,
      relative_entropy = relative_entropy, normalise.scores = normalise.scores,
      min.position.ic = min.position.ic, score.strat = score.strat.merge)))

  .internal_convert(final.mots, unique(CLASS_IN))

}

sort_by_ic <- function(x) {
  # Should help with getting the best overall alignments.
  ICs <- vapply(x, function(x) x@icscore, numeric(1))
  x[order(ICs, decreasing = TRUE)]
}
