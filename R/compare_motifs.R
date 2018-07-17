#' Compare motifs.
#'
#' TFBSTools implementation of motif comparison, rewritten in C++. See
#' \code{\link[TFBSTools]{PWMSimilarity}} for details.
#'
#' @param motifs List of motifs. If not a \linkS4class{universalmotif} object,
#'    they will be converted. See \code{\link{convert_motifs}} for supported
#'    classes.
#' @param method One of 'Euclidean', 'Pearson', and 'KL'.
#' @param k Numeric. If '1', compare regular motif matrices. Otherwise, compare
#'    k-freq multifreq slot.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#' @param progress_bar Logical. Show progress bar.
#'
#' @return Distance matrix for Euclidean or KL; similarity matrix for Pearson.
#'    Each distance/similarity score is between 0 and 1.
#'
#' @details
#'    The implementations of this function are the exact same as that of the
#'    \code{\link[TFBSTools]{PWMSimilarity}} function, except that this function
#'    allows for more than two motifs to be compared as well as providing
#'    significant performance gains.
#'
#'    Each possible motif pairs are compared to generate the final matrix. This
#'    is done by first aligning the two motifs, then performing the
#'    distance/similarity calculation for each position pairs. If the two
#'    motifs are not the same length, then the calculation is performed
#'    repeatedly by moving the smaller motif along the larger motif. Afterwards,
#'    either the smallest distance or the largest similarity is reported.
#'
#' @examples
#'    motif1 <- create_motif()
#'    motif2 <- create_motif()
#'    motif1vs2 <- compare_motifs(list(motif1, motif2))
#'    # to get a dist object:
#'    as.dist(motif1vs2)
#'
#' @references
#'    \insertRef{tfbstools}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso \code{\link{convert_motifs}}, \code{\link[TFBSTools]{PWMSimilarity}},
#'    \code{\link{motif_tree}}
#' @export
compare_motifs <- function(motifs, method = "Euclidean", k = 1,
                           BPPARAM = SerialParam(),
                           progress_bar = FALSE) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)

  mot.names <- vapply(motifs, function(x) x["name"], character(1))

  pb_prev <- BPPARAM$progressbar
  if (progress_bar) BPPARAM$progressbar <- TRUE

  if (k > 1) {

    for (i in seq_along(motifs)) {
      if (is.null(motifs[[i]]["multifreq"][[as.character(k)]])) {
        stop("The following motif is missing ", k, "-letter freqs: ",
             motifs[[i]]["name"])
      }
    }

    .compare_k <- function(x) {
      y <- vector("numeric", length = length(motifs))
      for (j in seq_along(motifs)) {
        y[j] <- motif_simil_internal(motifs[[x]]["multifreq"][[as.character(k)]],
                                     motifs[[j]]["multifreq"][[as.character(k)]],
                                     method)
      }
      data.frame(y)
    }

    comparisons <- bplapply(seq_along(motifs), .compare_k, BPPARAM = BPPARAM)

  } else {

    .compare <- function(x) {
      y <- vector("numeric", length = length(motifs))
      for (j in seq_along(motifs)) {
        y[j] <- motif_simil_internal(motifs[[x]]["motif"], motifs[[j]]["motif"],
                                     method)
      }
      data.frame(y)
    }

    comparisons <- bplapply(seq_along(motifs), .compare, BPPARAM = BPPARAM)
  }

  comparisons <- unlist(comparisons)
  comparisons <- matrix(comparisons, ncol = length(motifs))

  colnames(comparisons) <- mot.names
  rownames(comparisons) <- mot.names

  BPPARAM$progressbar <- pb_prev

  comparisons

}
