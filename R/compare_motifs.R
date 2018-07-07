#' Compare motifs.
#'
#' TFBSTools implementation of motif comparison. Rewritten in C++.
#'
#' @param motifs List of motifs.
#' @param method One of Euclidean, Pearson, and KL.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Distance matrix for Euclidean and KL; similarity matrix for Pearson.
#'
#' @references
#'    \insertRef{tfbstools}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
compare_motifs <- function(motifs, method = "Euclidean", BPPARAM = bpparam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)

  mot.names <- vapply(motifs, function(x) x["name"], character(1))
  .compare <- function(x) {
    y <- vector("numeric", length = length(motifs))
    for (j in seq_along(motifs)) {
      y[j] <- motif_simil_internal(motifs[[x]]["motif"], motifs[[j]]["motif"],
                                   method)
    }
    data.frame(y)
  }

  comparisons <- bplapply(seq_along(motifs), .compare, BPPARAM = BPPARAM)
  comparisons <- unlist(comparisons)
  comparisons <- matrix(comparisons, ncol = length(motifs))

  colnames(comparisons) <- mot.names
  rownames(comparisons) <- mot.names

  comparisons

}
