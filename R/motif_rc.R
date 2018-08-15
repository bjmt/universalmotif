#' Get the reverse complement of a motif.
#'
#' For any motif, change the 'motif' slot to it's reverse complement. If the
#' 'multifreq' slot is filled, then it is also applied. No other slots are
#' affected.
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
motif_rc <- function(motifs, BPPARAM = SerialParam()) {

  if (is.list(motifs)) {
    motifs <- bplapply(motifs, motif_rc, BPPARAM = BPPARAM)
    return(motifs)
  }

  CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)

  multifreq <- motifs@multifreq
  if (length(multifreq) > 0) {
    for (i in seq_along(multifreq)) {
      rown.prev <- rownames(multifreq[[i]])
      multifreq[[i]] <- matrix(rev(as.numeric(multifreq[[i]])),
                               nrow = nrow(multifreq[[i]]))
      rownames(multifreq[[i]]) <- rown.prev
      colnames(multifreq[[i]]) <- seq_len(ncol(multifreq[[i]]))
    }
  }

  motifs <- universalmotif_cpp(name = motifs["name"], altname = motifs["altname"],
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
  if (length(multifreq) > 0) motif@multifreq <- multifreq

  msg <- validObject_universalmotif(motifs)
  if (length(msg) > 0) stop(msg)

  motifs <- .internal_convert(motifs, CLASS_IN, BPPARAM = BPPARAM)
  motifs

}
