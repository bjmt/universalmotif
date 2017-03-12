######################################################################
## Benjamin Tremblay
##
## Read UNIPROBE motifs from a text file.
##
######################################################################

#' [UNDER CONSTRUCTION] Load UNIPROBE motifs from a text file.
#'
#' Support for uniprobe-style motifs. Each motif has a header and a number
#' matrix.
#'
#' @param motif_file Character.
#' @param verbose Logical.
#' @param mot_length_cutoff Integer.
#' @param out_class Character.
#'
#' @return a list of motif objects of the desired class.
#'
#' @examples
#'    motifs <- system.file("extdata", "uniprobe.txt",
#'                          package = "universalmotif")
#'    rmotifs <- read_uniprobe(motifs)
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @include utils.R
#' @export
read_uniprobe <- function(motif_file, verbose = FALSE,
                          mot_length_cutoff = NULL, out_class = "matrix-2") {

  # check args
  check_logi_args(as.list(environment())[2])
  check_filter_args(as.list(environment())[3])
  check_out_class(out_class)

  # read file
  uniprobe_raw <- readLines(con <- file(motif_file)); close(con)
  if (length(uniprobe_raw) == 0) stop("could not read file, or file is empty",
                                      call. = FALSE)
  names(uniprobe_raw) <- seq_along(uniprobe_raw)

}
