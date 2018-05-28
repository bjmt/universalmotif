#' Motif similarity matrix.
#'
#' TFBSTools implementation of motif similarity. 
#'
#' @param motifs List of motifs.
#' @param method One of Euclidean, Pearson, KL.
#' @param BPPARAM Param for bplapply.
#
#'
#' @return Similarity matrix.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.simil <- motif_simil(jaspar)
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_simil <- function(motifs, method = "Pearson", BPPARAM = bpparam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)

  motif_types <- vapply(motifs, function(x) x["type"], character(1))
  if (any(motif_types != "PWM")) {
    motifs[motif_types != "PWM"] <- lapply(motifs[motif_types != "PWM"],
                                           function(x) convert_type(x, "PWM",
                                                          BPPARAM = BPPARAM))
  }

  motif_names <- vapply(motifs, function(x) x["name"], character(1))

  .similarities <- function(x) {
    y <- vector("numeric", length = length(motifs))
    for (j in seq_along(motifs)) {
      y[j] <- PWMSimilarity(motifs[[x]]["motif"], motifs[[j]]["motif"],
                            method = method)
    }
    data.frame(y)
  }
  similarities <- bplapply(seq_along(motifs), .similarities)
  similarities <- unlist(similarities)
  similarities <- matrix(similarities, ncol = length(motifs))

  colnames(similarities) <- motif_names
  rownames(similarities) <- motif_names

  similarities

}
