######################################################################
## Benjamin Tremblay
##
## Main motif reading function
##
######################################################################

#' [UNDER CONSTRUCTION] Load motifs from a text file.
#'
#' Motifs can be autodetected or explicitly named. Each file must only contain
#' one format; mixed motif files are not supported.
#'
#' @param motif_file Character.
#' @param verbose Logical.
#' @param motif_format Character.
#' @param out_class Character.
#' @param ... Additional arguments to be passed to the motif format-specific
#'    function.
#'
#' @return A list of motifs.
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @include utils.R
#' @export
read_motifs <- function(motif_file, verbose = FALSE,
                        motif_format = "autodetect", out_class = "matrix-2",
                        ...) {

  # check args
  check_logi_args(as.list(environment())[2])  # utils.R
  check_out_class(out_class)  # utils.R
  if (length(motif_format) != 1) stop("only one 'motif_format' can be used",
                                      call. = FALSE)
  if (!motif_format %in% c("autodetect", "meme", "jaspar", "homer")) {
    stop("'", motif_format, "' is not a supported motif format", call. = FALSE)
  }

  # do autodetection, if requested
  if (motif_format == "autodetect") {
    if (verbose) cat("Attempting to detect motif format..\n")
    motif_format <- find_format(motif_file, verbose)
  } 

  available_formats <- list(c("meme", "read_meme"), c("jaspar", "read_jaspar"),
                            c("homer", "read_homer"))

  # get function name
  final_format <- available_formats[vapply(available_formats, function(x)
                            motif_format == x[1], logical(1))]

  if (verbose) cat("Using '", final_format[[1]][2], "' to parse motifs..\n\n",
                   sep = "")

  # send args to correct function
  motifs <- do.call(final_format[[1]][2], c(list("motif_file" = motif_file,
                                                 "verbose" = verbose,
                                                 "out_class" = out_class),
                                            ...))

  return(motifs)

}

######################################################################
######################################################################

find_format <- function(motif_file, verbose) {

  # read file
  con <- file(motif_file)
  motifs_raw <- readLines(con)
  close(con)
  if (length(motifs_raw) == 0) stop("could not read file, or file is empty",
                                    call. = FALSE)
  names(motifs_raw) <- seq_along(motifs_raw)

  motif_format <- "meme"

  if (verbose) cat("Detected as: '", motif_format, "'\n", sep = "")

  return("meme")

}
