#' Generate binding sites from a motif.
#'
#' Given probabilities for a sequence as represented by a motif, generate
#' random sequences with the same length as the motif.
#'
#' @param motif See [convert_motifs()] for acceptable formats.
#' @param n `numeric(1)` Number of sites to generate.
#' @param use.freq `numeric(1)` If one, use regular motif matrix. Otherwise,
#'    use respective `multifreq` matrix.
#'
#' @return [Biostrings::XStringSet-class] object.
#'
#' @examples
#' motif <- create_motif()
#' sites <- sample_sites(motif)
#'
#' @seealso [create_sequences()], [create_motif()], [add_multifreq()]
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
sample_sites <- function(motif, n = 100, use.freq = 1) {

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list(n = args$n), 1, FALSE, "numeric")
  all_checks <- c(num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motif <- convert_motifs(motif)
  motif <- convert_type(motif, "PPM")

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

  sites <- lapply(seq_len(n), .get_sites)
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
