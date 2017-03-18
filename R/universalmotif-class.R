######################################################################
## Benjamin Tremblay
##
## S4 class for storing motifs
##
######################################################################

#' [UNDER CONSTRUCTION] universalmotif: Motif class.
#'
#' Placeholder description.
#'
#' @slot name Character. Motif name.
#' @slot motif Matrix. Contains motif.
#' @slot alphabet Character. Describes alphabet.
#' @slot letters Character. Alphabet letters.
#' @slot type Character. Describes how the motif is stores; i.e., as a PWM,
#'       etc.
#' @slot nsites Integer. Number of times motif was found in original set.
#' @slot penrich Numeric. Motif enrichment p-value, if present.
#' @slot pseudoweights Numeric vector. Motif pseudoweights.
#' @slot eval Numeric. Motif e-value, if present.
#' @slot pdetect Numeric. Motif detection p-value, if present.
#' @slot tpos Numeric. Average position of motif in target sequences.
#' @slot tstd Numeric. Standard deviation of position of motif in target
#'       sequences.
#' @slot bpos Numeric. Average position of motif in background sequences.
#' @slot bstd Numeric. Standard deviation of average position in background
#'       sequences.
#' @slot strandbias Numeric. Log ratio of + strand occurences to - strand
#'       occurrences.
#' @slot tpercent Numeric. Percent of target sequences with motif.
#' @slot bpercent Numeric. Percent of background sequences with motif.
#' @slot consensus Character. Motif consensus sequences.
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @export
setClass("universalmotif", slots = list(name = "character", motif = "matrix",
                              alphabet = "character", letters = "character",
                              type = "character", nsites = "integer",
                              penrich = "numeric", pseudoweights = "numeric",
                              eval = "numeric", pdetect = "numeric",
                              tpos = "numeric", tstd = "numeric",
                              bpos = "numeric", bstd = "numeric",
                              strandbias = "numeric", tpercent = "numeric",
                              bpercent = "numeric", consensus = "character"))

# idea: make universalmotif a virtual class with fewer slots
