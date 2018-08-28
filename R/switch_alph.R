#' Switch between DNA and RNA alphabets.
#'
#' Convert a motif from DNA to RNA, or RNA to DNA.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#'
#' @return The DNA/RNA version of the motifs.
#'
#' @examples
#' DNA.motif <- create_motif()
#' RNA.motif <- switch_alph(DNA.motif)
#'
#' @seealso [create_motif()]
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
switch_alph <- function(motifs) {

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
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

  motifs <- lapply(motifs, .switch_alph)

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs <- .internal_convert(motifs, CLASS_IN)
  motifs

}
