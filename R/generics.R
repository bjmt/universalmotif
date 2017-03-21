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
setGeneric("motif_slots", function(object, slots) standardGeneric("motif_slots"))

#' @export
setGeneric("motif_slots<-", function(object, slot, value)
           standardGeneric("motif_slots<-"))

#' [UNDER CONSTRUCTION] Convert motifs between various classes.
#'
#' Convert motif objects between various supported classes; some of these
#' require their respective packages being installed.
#'
#' @param motif Motif object or list of motif objects.
#' @param out_class Character. Desired motif class to convert to.
#'
#' @return Motif object in the desired class, or list of motifs.
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
#' Convert between PCM, PPM, PWM and ICM. Cannot convert back from ICM.
#'
#' @param motif Motif object.
#' @param out_type Character. Either PCM, PPM, PWM or ICM.
#' @param pseudoweight Numeric. Leaving this to NULL will cause the function
#'    to use the pseudoweight found within the motif object; set this to 0
#'    to prevent information loss between conversions.
#' @param background Numeric. Leaving this to NULL will cause the function
#'    to use the background frequencies found within the motif object. NOTE:
#'    when converting to 'ICM', the function will use a uniform
#'    background frequency of 0.25 if NULL (otherwise this can result in IC
#'    higher than 2!). The function can be forced to use different
#'    background frequencies if given.
#' @param IC_floor Logical. Whether or not to apply a floor of 0 to the IC
#'    results.
#'
#' @return Motif object of class \linkS4class{universalmotif}.
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @export
setGeneric("convert_type", function(motif, out_type, pseudoweight = NULL,
                                    background = NULL, IC_floor = FALSE)
           standardGeneric("convert_type"))

#' [UNDER CONSTRUCTION] Create motifs.
#'
#' Manually generate motifs based on a consensus string, matrix, or set of
#' sequences of matching lengths.
#'
#' @param consensus Character. Consensus string.
#' @param matrix Matrix.
#' @param sequences DNAStringSet or RNAStringSet. Set of sequences of
#'    identical length.
#' @param out_type Character. Either PCM, PPM, PWM, or ICM.
#' @param out_class Character.
#' @param pseudoweight Numeric.
#' @param alphabet Character. Either DNA or RNA.
#' @param background Numeric. Must be of length 4.
#' @param nsites Numeric. If left as-is, will assume nsites as 100 when
#'    needed.
#'
#' @return Motif object of the desired class.
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @export
setGeneric("create_motif", function(consensus, matrix, sequences,
                                     name, out_type = "PPM",
                                     out_class = "universalmotif",
                                     pseudoweight = 0.8, alphabet = "DNA",
                                     background = c(0.25, 0.25, 0.25, 0.25),
                                     nsites = numeric(0))
           standardGeneric("create_motif"))
