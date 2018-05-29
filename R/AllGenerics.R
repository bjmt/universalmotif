#' Create a motif.
#'
#' @param input One of character vector, matrix, or XStringSet.
#' @param name Character.
#' @param pseudoweight Numeric.
#' @param alphabet Character.
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
#' my.motifs <- create_motif("TATAA")
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("create_motif", function(input,
                                    name = "motif", pseudoweight = 0.8,
                                    alphabet, bkg, nsites, altname, family,
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
