#' Motif similarity matrix.
#'
#' TFBSTools implementation of motif similarity. Supports parallel execution
#' via foreach package.
#'
#' @param motifs List of motifs.
#' @param method One of Euclidean, Pearson, KL.
#'
#' @return Similarity matrix.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_simil <- function(motifs, method = "Pearson") {

  motifs <- convert_motifs(motifs)

  motif_types <- vapply(motifs, function(x) x["type"], character(1))
  if (any(motif_types != "PWM")) {
    motifs[motif_types != "PWM"] <- lapply(motifs[motif_types != "PWM"],
                                           function(x) convert_type(x, "PWM"))
  }

  motif_names <- vapply(motifs, function(x) x["name"], character(1))

  similarities <- foreach(i = seq_along(motifs), .combine = "cbind") %:%
    foreach(j = seq_along(motifs), .combine = "c") %dopar% {
      PWMSimilarity(motifs[[i]]["motif"], motifs[[j]]["motif"],
                    method = method)
    }

  colnames(similarities) <- motif_names
  rownames(similarities) <- motif_names

  similarities

}
