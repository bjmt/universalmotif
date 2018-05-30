#' Create a motif.
#'
#' @param input One of character vector, matrix (PCM, PPM, or PWM), or 
#'              XStringSet.
#' @param alphabet Character.
#' @param type Character.
#' @param name Character.
#' @param pseudoweight Numeric.
#' @param bkg Numeric.
#' @param nsites Numeric.
#' @param altname Character.
#' @param family Character.
#' @param organism Character.
#' @param bkgsites Numeric.
#' @param strand Character.
#' @param pval Numeric.
#' @param qval Numeric.
#' @param eval Numeric.
#' @param extrainfo Character.
#'
#' @return Motif object.
#'
#' @examples
#' ## create motifs from a single string
#' DNA.motif <- create_motif("TATAAWW")
#' RNA.motif <- create_motif("UUUCCG")
#' AA.motif <- create_motif("AVLK", alphabet = "AA")
#' custom.motif <- create_motif("QWER", alphabet = "custom")
#' custom.motif <- create_motif("QWER", alphabet = "QWERASDF")
#'
#' ## create motifs from multiple strings of equal length
#' 
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("create_motif", function(input, alphabet, type,
                                    name = "motif", pseudoweight = 0.8,
                                    bkg, nsites, altname, family,
                                    organism, bkgsites, strand, pval, qval,
                                    eval, extrainfo)
           standardGeneric("create_motif"))

#' Convert motif class.
#'
#' @param motifs Single motif object or list.
#' @param class Desired motif class. E.g. 'motifStack-pfm'.
#' @param BPPARAM Param for bplapply.
#'
#' @return Single motif object or list.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.motifstack.pfm <- convert_motifs(jaspar, "motifStack-pfm")
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("convert_motifs", function(motifs,
                                      class = "universalmotif-universalmotif",
                                      BPPARAM = bpparam())
           standardGeneric("convert_motifs"))
