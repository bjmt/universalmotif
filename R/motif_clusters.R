#' Find clusters of motifs in sequences.
#'
#' @param motifs See `convert_motifs()` for acceptable motif formats.
#' @param min.dist.inner `integer(1)` Min distance between motifs.
#' @param min.dist.outer `integer(1)` Min distance between motifs (including
#'    the size of the motifs).
#' @param max.dist.inner `integer(1)` Max distance between motifs.
#' @param max.dist.outer `integer(1)` Max distance between motifs (including
#'    the size of the motifs).
#' @param min.size `integer(1)` Min number of motifs per cluster.
#' @param max.size `integer(1)` Max number of motifs per cluster.
#' @param allow.dup `logical(1)` Allow motifs to occur more than once within a
#'    cluster.
#' @param allow.disorder `logical(1)` Ignore the order of motifs within a
#'    cluster between occurrences of the cluster.
#' @param allow.overlap `logical(1)` Allow motifs within a cluster to overlap.
#' @param scan.threshold `numeric(1)` See [scan_sequences()].
#' @param scan.threshold.type `character(1)` See [scan_sequences()].
#' @param scan.use.freq `integer(1)` See [scan_sequences()].
#'
#' @return `list`
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @seealso [scan_sequences()], [enrich_motifs()]
#' @inheritParams scan_sequences
#' @noRd
motif_clusters <- function(motifs, sequences,
                           min.dist.inner = 0, min.dist.outer = 20,
                           max.dist.inner = 20, max.dist.outer = 100,
                           min.size = 2, max.size = 5, allow.dup = TRUE,
                           allow.disorder = TRUE, allow.overlap = FALSE,
                           RC = FALSE, scan.threshold = 0.001,
                           scan.threshold.type = "pvalue",
                           scan.use.freq = 1, motif_pvalue.k = 8,
                           nthreads = 1) {

  scan.res <- scan_sequences(motifs, sequences, threshold = scan.threshold,
                             threshold.type = scan.threshold.type, RC = RC,
                             use.freq = scan.use.freq, nthreads = nthreads,
                             motif_pvalue.k = motif_pvalue.k)



}
