######################################################################
## Benjamin Tremblay
##
## Generics
##
######################################################################

#' Accessor function for class universalmotif.
#'
#' Easy access to universalmotif slots.
#'
#' @param object Motif object of class \linkS4class{universalmotif}.
#'
#' @return List of motif info.
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @export
setGeneric("motif_slots", function(object, ...) standardGeneric("motif_slots"))

#' [UNDER CONSTRUCTION] Convert motifs between various classes.
#'
#' Convert motif objects between various supported classes; some of these
#' require their respective packages being installed.
#'
#' @param motif Motif object or list of motif objects..
#' @param out_class Character. Desired motif class to convert to.
#'
#' @return Motif object in the desired class, or list of motifs..
#'
#' @examples
#'   motifs <- system.file("extdata", "minimal.meme", package = "universalmotif")
#'   motifs <- read_meme(motifs, out_class = "matrix-2")
#'   motifs <- convert_motifs(motifs, out_class = "matrix-1")
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @export
setGeneric("convert_motifs", function(motif, out_class, ...) 
           standardGeneric("convert_motifs"))

#' [UNDER CONSTRUCTION] Convert \linkS4class{universalmotif} type.
#'
#' Convert between PCM, PPM and PWM.
#'
#' @param motif Motif object.
#' @param out_type Character. Either PCM, PPM, or PWM.
#'
#' @return Motif object of class \linkS4class{universalmotif}.
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @export
setGeneric("convert_type", function(motif, out_type, ...)
           standardGeneric("convert_type"))
