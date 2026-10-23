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

  ## Each k-mer row takes the values of its reverse complement, in the
  ## reverse column order.
  multifreq <- motifs@multifreq
  if (length(multifreq) > 0) {
    comp <- if (motifs@alphabet == "RNA") c("ACGU", "UGCA") else c("ACGT", "TGCA")
    for (i in seq_along(multifreq)) {
      mat <- multifreq[[i]]
      rc.names <- vapply(rownames(mat), function(x)
        collapse_cpp(rev(safeExplode(chartr(comp[1], comp[2], x)))), character(1))
      mat <- mat[rc.names, rev(seq_len(ncol(mat))), drop = FALSE]
      dimnames(mat) <- dimnames(multifreq[[i]])
      multifreq[[i]] <- mat
    }
  }

  ## A gap after column g of a motif of width n falls after column n - g of
  ## its reverse complement.
  gapinfo <- motifs@gapinfo
  if (length(gapinfo@gaploc) > 0) {
    o <- rev(seq_along(gapinfo@gaploc))
    gapinfo@gaploc <- ncol(motifs@motif) - gapinfo@gaploc[o]
    gapinfo@mingap <- gapinfo@mingap[o]
    gapinfo@maxgap <- gapinfo@maxgap[o]
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

  if (length(multifreq) > 0) motifs@multifreq <- multifreq
  motifs@gapinfo <- gapinfo

  validObject_universalmotif(motifs)
  motifs

}
