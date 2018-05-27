#' Create a motif.
#'
#' @param consensus Character. Consensus string.
#' @param matrix Matrix.
#' @param sequences XStringSet.
#' @param name Character.
#' @param pseudoweight Numeric.
#' @param alphabet Character.
#' @param bkg Numeric.
#' @param nsites Numeric.
#'
#' @return Motif object.
#'
#' @examples
#' my.motifs <- create_motif("TATAA")
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("create_motif", function(consensus, matrix, sequences,
                                    name = "motif", pseudoweight = 0.8,
                                    alphabet, bkg, nsites)
           standardGeneric("create_motif"))

#' Convert motif class.
#'
#' @param motifs Single motif object or list.
#' @param class Desired motif class. E.g. 'motifStack-pfm'.
#'
#' @return Single motif object or list.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   pacakge = "universalmotif"))
#' jaspar.motifstack.pfm <- convert_motifs(jaspar, "motifStack-pfm")
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("convert_motifs", function(motifs,
                                      class = "universalmotif-universalmotif")
           standardGeneric("convert_motifs"))
