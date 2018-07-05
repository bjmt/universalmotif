#' Switch between DNA and RNA alphabets.
#'
#' @param motifs Motif object or list of.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return The DNA/RNA version of the motif(s).
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
switch_alph <- function(motifs, BPPARAM = bpparam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  if (!is.list(motifs)) motifs <- list(motifs)

  .switch_alph <- function(motif) {
    old.alph <- motif["alphabet"]
    if (!old.alph %in% c("DNA", "RNA")) {
      stop("motif '", motif["name"], "' has alphabet '", old.alph,
           "'; will not convert")
    }
    if (old.alph == "DNA") new.alph <- "RNA" else new.alph <- "DNA"
    if (new.alph == "RNA") {
      motif@alphabet <- "RNA"
      rownames(motif@motif) <- RNA_BASES
      colnames(motif@motif) <- gsub("T", "U", colnames(motif@motif))
      motif@consensus <- gsub("T", "U", motif@consensus)
      if (length(motif@multifreq) > 0) {
        for (i in names(motif@multifreq)) {
          rownames(motif@multifreq[[i]]) <- gsub("T", "U", rownames(motif@multifreq[[i]]))
        }
      }
    } else {
      motif@alphabet <- "DNA"
      rownames(motif@motif) <- DNA_BASES
      colnames(motif@motif) <- gsub("U", "T", colnames(motif@motif))
      motif@consensus <- gsub("U", "T", motif@consensus)
      if (length(motif@multifreq) > 0) {
        for (i in names(motif@multifreq)) {
          rownames(motif@multifreq[[i]]) <- gsub("U", "T", rownames(motif@multifreq[[i]]))
        }
      }
    }
    motif
  }

  motifs <- bplapply(motifs, .switch_alph, BPPARAM = BPPARAM)

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs

}
