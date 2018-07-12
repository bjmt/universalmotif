#' Generate binding sites from a motif.
#'
#' Given probabilities for a sequence as represented by a motif, generate
#' random sequences with the same length as the motif.
#'
#' @param motif Motif object.
#' @param n Numeric. Number of sites to generate.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return \linkS4class{XStringSet} object.
#'
#' @examples
#' motif <- create_motif()
#' sites <- sample_sites(motif)
#'
#' @seealso \code{\link{create_sequences}}, \code{\link{create_motif}},
#'    \code{\link{add_multifreq}}
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
sample_sites <- function(motif, n = 100, BPPARAM = SerialParam()) {

  motif <- convert_motifs(motif, BPPARAM = BPPARAM)
  motif <- convert_type(motif, "PPM", BPPARAM = BPPARAM)

  mot.mat <- motif["motif"]
  alph <- motif["alphabet"]

  .get_sites <- function(x) {
    site <- apply(mot.mat, 2, function(x) sample(rownames(mot.mat), 1, TRUE, x))
    site <- paste(site, collapse = "")
    site
  }

  sites <- bplapply(seq_len(n), .get_sites, BPPARAM = BPPARAM)
  sites <- unlist(sites)

  if (alph == "DNA") {
    sites <- DNAStringSet(sites)
  } else if (alph == "RNA") {
    sites <- RNAStringSet(sites)
  } else if (alph == "AA") {
    sites <- AAStringSet(sites)
  } else sites <- BStringSet(sites)

  sites

}
