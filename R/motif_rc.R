#' Get the reverse complement of a DNA or RNA motif.
#'
#' For any motif, change the `motif` slot to it's reverse complement. If the
#' `multifreq` slot is filled, then it is also applied. No other slots are
#' affected.
#'
#' @param motifs See [convert_motifs()] for acceptable formats
#' @param ignore.alphabet `logical(1)` If `TRUE`, then [motif_rc()] throws
#'    an error when it detects a non-DNA/RNA motif. If `FALSE`, it will
#'    proceed regardless.
#'
#' @return See [convert_motifs()] for available output formats.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.rc <- motif_rc(jaspar)
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
motif_rc <- function(motifs, ignore.alphabet = FALSE) {

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
  if (is.list(motifs)) was.list <- TRUE else was.list <- FALSE
  if (!is.list(motifs)) motifs <- list(motifs)

  motifs <- lapply(motifs, motif_rc_internal, ignore.alphabet = ignore.alphabet)

  motifs <- .internal_convert(motifs, unique(CLASS_IN))
  if (length(motifs) == 1 && !was.list)
    motifs <- motifs[[1]]

  motifs

}

motif_rc_internal <- function(motifs, ignore.alphabet = FALSE) {

  if (!motifs@alphabet %in% c("DNA", "RNA") && !ignore.alphabet)
    stop("Only DNA/RNA motifs are allowed")

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

  motifs <- universalmotif_cpp(name = motifs@name, altname = motifs@altname,
                           family = motifs@family,
                           motif = matrix(rev(as.numeric(motifs@motif)),
                                          nrow = 4, byrow = FALSE),
                           alphabet = motifs@alphabet, type = motifs@type,
                           organism = motifs@organism,
                           nsites = motifs@nsites,
                           pseudocount = motifs@pseudocount,
                           bkg = motifs@bkg, bkgsites = motifs@bkgsites,
                           strand = motifs@strand, pval = motifs@pval,
                           qval = motifs@qval, eval = motifs@eval,
                           extrainfo = motifs@extrainfo)

  if (length(multifreq) > 0) motif@multifreq <- multifreq

  validObject_universalmotif(motifs)
  motifs

}
