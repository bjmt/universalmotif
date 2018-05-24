#' Convert motif class.
#'
#' @param motifs Single motif object or list.
#' @param class Desired motif class. E.g. 'motifStack-pfm'.
#'
#' @return Single motif object or list.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
convert_motifs <- function(motifs, class = "universalmotif") {

  if (class(motifs) == "list") {
    motifs <- lapply(motifs, .convert_motifs, class = class)
  } else {
    motifs <- .convert_motifs(motifs, class = class)
  }

  motifs

}

# justification for not using a generic:
#   (maybe I could just create a generic for 'ANY')
#   - multiple packages export classes with identical names
#   - MotIV package does not export the pwm2 class despite using it

.convert_motifs <- function(motifs, class) {

  out_class <- strsplit(class, "-")[[1]][2]
  out_class_pkg <- strsplit(class, "-")[[1]][1]
  in_class <- class(motifs)[1]
  in_class_pkg <- attributes(class(motifs))$package

  ## convert in_class to universalmotif:

  if (in_class_pkg == "MotifDb" && in_class == "MotifList") {
    motifs_out <- list()
    motifdb_fun <- function(x) {
      universalmotif(name = x@elementMetadata@listData$providerName,
                     altname = x@elementMetadata@listData$geneSymbol,
                     family = x@elementMetadata@listData$tfFamily,
                     organism = x@elementMetadata@listData$organism,
                     motif = x@listData[[1]], alphabet = "DNA",
                     type = "PPM")
    }
    for (i in seq_len(length(motifs))) {
      motifs_out[[i]] <- motifdb_fun(motifs[i])
    }
    motifs <- lapply(motifs_out, function(x) .convert_motifs(x, class = class))
    return(motifs)
  }

  ## convert universalmotif to out_class:

  if (out_class_pkg == "MotIV" && out_class == "pwm2") {
    motifs <- convert_type(motifs, "PPM")
    motifs <- MotIV::makePWM(motifs["motif"], alphabet = motifs["alphabet"])
  }

  motifs

}
