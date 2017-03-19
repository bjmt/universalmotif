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
#' @param object Motif object of class universalmotif.
#'
#' @return List of motif info.
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @export
setGeneric("motif_slots", function(object) standardGeneric("motif_slots"))
