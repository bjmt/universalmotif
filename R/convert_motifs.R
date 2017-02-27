######################################################################
## Benjamin Tremblay
##
## Convert motifs between different classes
##
######################################################################

# IDEAS:
  
  # convert all motifs to a new S4 class which holds every single piece of info
  # then convert out to desired class

#' @title [UNDER CONSTRUCTION] Convert motifs between various classes.
#'
#' @description
#' Convert motif objects between various supported classes; some of these
#' require their respective packages being installed.
#'
#' @param motif_list List of motif objects.
#' @param out_class Character. Desired motif class to convert to.
#'
#' @return List of motif objects in the desired class.
#'
#' @examples
#'   motifs <- system.file("extdata", "minimal.meme", package = "universalmotif")
#'   motifs <- read_meme(motifs, show_warnings = FALSE, out_class = "matrix-2")
#'   motifs <- convert_motifs(motifs, out_class = "matrix-1")
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @export
convert_motifs <- function(motif_list, out_class = "matrix=1") {
  print("hello")
}

######################################################################


