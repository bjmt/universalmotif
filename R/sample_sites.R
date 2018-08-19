#' Generate binding sites from a motif.
#'
#' Given probabilities for a sequence as represented by a motif, generate
#' random sequences with the same length as the motif.
#'
#' @param motif See \code{\link{convert_motifs}} for acceptable formats.
#' @param n \code{numeric(1)} Number of sites to generate.
#' @param use.freq \code{numeric(1)} If one, use regular motif matrix. Otherwise,
#'    use respective \code{multifreq} matrix.
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
sample_sites <- function(motif, n = 100, use.freq = 1, BPPARAM = SerialParam()) {

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list(n = args$n), 1, FALSE, "numeric")
  s4_check <- check_fun_params(list(BPPARAM = args$BPPARAM), numeric(), FALSE, "S4")
  all_checks <- c(num_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motif <- convert_motifs(motif, BPPARAM = BPPARAM)
  motif <- convert_type(motif, "PPM", BPPARAM = BPPARAM)

  if (use.freq == 1) {
    mot.mat <- motif["motif"]
  } else {
    mot.mat <- motif["multifreq"][[as.character(use.freq)]]
  }

  alph <- motif["alphabet"]

  .get_sites <- function(x) {
    site <- apply(mot.mat, 2, function(x) sample(rownames(mot.mat), 1, TRUE, x))
    sites <- collapse_cpp(site)
    sites
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
