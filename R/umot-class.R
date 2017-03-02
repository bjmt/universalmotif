######################################################################
## Benjamin Tremblay
##
## S4 class for storing motifs
##
######################################################################

#' [UNDER CONSTRUCTION] umot: Motif class.
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
#' @slot pval Double. Motif p-value, if present.
#' @slot pseudoweights Numeric vector. Motif pseudoweights.
#' @slot eval Double. Motif e-value, if present.
#'
#' @export
setClass("umot", slots = list(name = "character", motif = "matrix",
                              alphabet = "character", letters = "character",
                              type = "character", nsites = "integer",
                              pval = "numeric", pseudoweights = "numeric",
                              eval = "numeric"))
