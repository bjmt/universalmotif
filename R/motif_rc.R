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

  hmmfirst <- motifs["hmmfirst"]
  if (length(hmmfirst) > 1) {
    hmmfirst.dims <- dimnames(hmmfirst)
    hmmfirst <- matrix(rev(as.numeric(hmmfirst)), nrow = 16, byrow = FALSE)
    dimnames(hmmfirst) <- hmmfirst.dims
  }
  hmmsecond <- motifs["hmmsecond"]
  if (length(hmmsecond) > 1) {
    hmmsecond.dims <- dimnames(hmmsecond)
    hmmsecond <- matrix(rev(as.numeric(hmmsecond)), nrow = 64, byrow = FALSE)
    dimnames(hmmsecond) <- hmmsecond.dims
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
                           extrainfo = motifs["extrainfo"],
                           hmmfirst = hmmfirst,
                           hmmsecond = hmmsecond)

  motifs <- .internal_convert(motifs, CLASS_IN, BPPARAM = BPPARAM)
  motifs

}
