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
#' @param use_warning Logical.
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
                        use_warning = FALSE, ...) {

  # check args
  check_logi_args(as.list(environment())[c(2, 5)])  # utils.R
  check_out_class(out_class)  # utils.R
  if (length(motif_format) != 1) stop("only one 'motif_format' can be used",
                                      call. = FALSE)

  available_formats <- list(c("MEME", "read_meme"), c("Jaspar", "read_jaspar"),
                            c("Homer", "read_homer"))

  motif_args <- c(list("motif_file" = motif_file, "verbose" = verbose,
                       "out_class" = out_class), ...)

  # handling autodetection
  if (motif_format == "autodetect") {

    if (verbose) cat("Attempting to detect motif format..\n")
    motif_format <- find_format(motif_file, verbose)
    final_format <- available_formats[vapply(available_formats, function(x)
                                      any(motif_format == x[1]), logical(1))]

    for (i in final_format) {

      if (verbose) cat("Detected as: '", i[1], "'\n", sep = "")
      if (verbose) cat("Using '", i[2], "' to parse motifs..\n",
                     sep = "")

      motifs <- extract_motifs(i[2], motif_args = motif_args,
                               verbose = verbose)

      if (!is.null(motifs)) return(motifs)

    }

  } else {

    final_format <- available_formats[vapply(available_formats, function(x)
                                      motif_format == x[1], logical(1))]
    if (length(final_format) == 0) {
      stop("\"", motif_format, "\" is not a supported motif format", call. = FALSE)
    }

    motifs <- do.call(final_format[[1]][2], args = motif_args) 
    return(motifs)

  }

  # inform user of autodetection failure
  if (!use_warning) stop("autodetection failed", call. = FALSE)
  
  warning("autodetection failed", call. = FALSE)

  return(invisible(NULL))

}

######################################################################
######################################################################

find_format <- function(motif_file, verbose, start_detection = NULL) {

  # read file
  con <- file(motif_file)
  motifs_raw <- readLines(con)
  close(con)
  if (length(motifs_raw) == 0) stop("could not read file, or file is empty",
                                    call. = FALSE)
  names(motifs_raw) <- seq_along(motifs_raw)

  ##########

  # autodetection algorithm goes here

  motif_format <- c("Jaspar", "MEME")

  ##########
  
  return(motif_format)

}

extract_motifs <- function(final_format, motif_args, verbose) {
  motifs <- try(do.call(final_format, args = motif_args), silent = TRUE)
  if (class(motifs) == "try-error") {
    if (verbose) cat("Parsing of motifs with '", final_format, "' failed!",
                     " Continuing autodetection..\n", sep = "") 
    return(NULL)
  } else return(motifs)

}
