######################################################################
## Benjamin Tremblay
##
## Generics
##
######################################################################

#' Accessor function for \linkS4class{universalmotif}.
#'
#' Easy access to \linkS4class{universalmotif} slots.
#'
#' @param object Motif object.
#' @param slots Chraracter. One or more motif slots to be accessed. See
#'    \linkS4class{universalmotif} for a description of the slots.
#'
#' @examples
#'    motifs <- system.file("extdata", "minimal.meme",
#'                          package = "universalmotif")
#'    motifs <- read_motifs(motifs)
#'    motif <- motifs[[1]]
#'    ## expose the 'name' slot:
#'    motif_slots(motif, "name")
#'    ## expose all filled slots as a list:
#'    motif_slots(motif)
#'    ## expose multiple slots:
#'    motif_slots(motif, c("icscore", "pseudoweight"))
#'
#' @return List of motif info.
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @export
setGeneric("motif_slots", function(object, slots) standardGeneric("motif_slots"))

#' Replacement function for \linkS4class{universalmotif}.
#'
#' Replace individual slots for \linkS4class{universalmotif}.
#'
#' @param object Motif object.
#' @param slot Character. Single slot to replace. See
#'    \linkS4class{universalmotif} for a description of the slots.
#' @param value Information to replace the slot with.
#'
#' @return Motif object.
#'
#' @examples
#'    motifs <- system.file("extdata", "minimal.meme",
#'                          package = "universalmotif")
#'    motifs <- read_motifs(motifs)
#'    motif <- motifs[[1]]
#'    ## change the name:
#'    motif_slots(motif, "name") <- "New_name"
#'    ## the 'motif' object is of 'type' PPM; changing it manually to PCM
#'    ## will cause an error to occur
#'    \dontrun{
#'       motif_slots(motif, "type") <- "PCM"
#'    }
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @seealso \code{\link{convert_type}} to convert motif 'type'
#' @export
setGeneric("motif_slots<-", function(object, slot, value)
           standardGeneric("motif_slots<-"))

#' [UNDER CONSTRUCTION] Convert motifs between various classes.
#'
#' Convert motif objects between various supported classes.
#'
#' @param motif Motif object or list of motif objects.
#' @param out_class Character. Desired motif class to convert to.
#' @param ... Additional parameters specified in the method for the
#'    class being acted upon.
#'
#' @return Motif object in the desired class, or list of motifs.
#'
#' @examples
#'   motifs <- system.file("extdata", "minimal.meme",
#'                         package = "universalmotif")
#'   motifs <- read_motifs(motifs)
#'   motifs <- convert_motifs(motifs, out_class = "TFBSTools-PFMatrix")
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @export
setGeneric("convert_motifs", function(motif, out_class = "universalmotif", ...) 
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
#' @examples
#'    motifs <- system.file("extdata", "minimal.meme",
#'                          package = "universalmotif")
#'    motifs <- read_motifs(motifs)
#'    motifs <- convert_type(motifs, out_type = "PWM")
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
#' @param sequences XStringSet. Set of sequences of identical length.
#' @param name Character. Required.
#' @param out_type Character. Either PCM, PPM, PWM, or ICM.
#' @param out_class Character.
#' @param pseudoweight Numeric.
#' @param alphabet Character. Currently either DNA or RNA.
#' @param background Numeric. Must be of length 4.
#' @param nsites Numeric. If left as-is, will assume nsites as 100 when
#'    needed.
#'
#' @return Motif object of the desired class.
#'
#' @examples
#'    motif <- create_motif(consensus = "AAVBWWWT", name = "Test_motif")
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
