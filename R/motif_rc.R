#' Get the reverse complement of a motif.
#'
#' @param motifs List of motifs or a single motif.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return List of motifs or single motif object.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.rc <- motif_rc(jaspar)
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_rc <- function(motifs, BPPARAM = bpparam()) {

  if (class(motifs) == "list") {
    motifs <- bplapply(motifs, motif_rc, BPPARAM = BPPARAM)
    return(motifs)
  }

  CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)

  if (length(motifs@multifreq) > 1) {
    warning("'multifreq' information is lost")
  }

  motifs <- universalmotif(name = motifs["name"], altname = motifs["altname"],
                           family = motifs["family"],
                           motif = matrix(rev(as.numeric(motifs["motif"])),
                                          nrow = 4, byrow = FALSE),
                           alphabet = motifs["alphabet"], type = motifs["type"],
                           organism = motifs["organism"],
                           nsites = motifs["nsites"],
                           pseudocount = motifs["pseudocount"],
                           bkg = motifs["bkg"], bkgsites = motifs["bkgsites"],
                           strand = motifs["strand"], pval = motifs["pval"],
                           qval = motifs["qval"], eval = motifs["eval"],
                           extrainfo = motifs["extrainfo"])

  motifs <- .internal_convert(motifs, CLASS_IN, BPPARAM = BPPARAM)
  motifs

}
