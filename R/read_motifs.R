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
                            c("Homer", "read_homer"),
                            c("Transfac", "read_transfac"))

  motif_args <- c(list("motif_file" = motif_file, "verbose" = verbose,
                       "out_class" = out_class), ...)

  # handling autodetection
  if (motif_format == "autodetect") {

    if (verbose) cat("Attempting to detect motif format..\n")
    motif_format <- find_format(motif_file, verbose, available_formats)
    final_format <- available_formats[vapply(available_formats, function(x)
                                      any(motif_format == x[1]), logical(1))]

    if (is.null(final_format)) {
      if (use_warning) {
        warning("autodetection failed", call. = FALSE)
        return(invisible(NULL))
      } else stop("autodetection failed", call. = FALSE)
    }

    for (i in final_format) {

      if (verbose) cat("Detected as: '", i[1], "'\n", sep = "")
      if (verbose) cat("Using '", i[2], "' to parse motifs..\n",
                     sep = "")

      motifs <- extract_motifs(i[2], motif_args = motif_args, verbose = verbose)

      if (!is.null(motifs)) return(motifs)

    }

  } else {

    # if no autodetection, then just move file along to correct function

    final_format <- available_formats[vapply(available_formats, function(x)
                                      grepl(paste0("^", motif_format, "$"), x[1],
                                            ignore.case = TRUE), logical(1))]
    if (length(final_format) == 0) {
      stop("\"", motif_format, "\" is not a supported motif format",
           call. = FALSE)
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

find_format <- function(motif_file, verbose, available_formats) {

  # read file
  con <- file(motif_file)
  motifs_raw <- readLines(con)
  close(con)
  if (length(motifs_raw) == 0) stop("could not read file, or file is empty",
                                    call. = FALSE)
  names(motifs_raw) <- seq_along(motifs_raw)

  final_format <- vector("character", length(available_formats))

  ##########

  # 1 MEME

  if (any(grepl("MEME version", motifs_raw,
                fixed = TRUE))) final_format[1] <- "MEME"

  # 2 Jaspar

  if (any(grepl("\\[", motifs_raw, fixed = TRUE)) &&
      any(grepl("\\]", motifs_raw, fixed = TRUE))) final_format[2] <- "Jaspar" 

  # 3 Transfac

  if (any(grepl("^XX$", motifs_raw))) final_format[3] <- "Transfac"

  # 4 Homer

  arrow_test <- motifs_raw[vapply(motifs_raw,
                                  function(x) grepl(">", x, fixed = TRUE),
                                  logical(1))]

  if (all(!is.na(vapply(arrow_test, function(x)
                        strsplit(x, split = "\\s+")[[1]][3], character(1))))) {
    final_format[4] <- "Homer"
  } 

  ##########

  final_format <- final_format[final_format != ""]

  if (length(final_format) == 0) return(NULL)

  return(final_format)

}

extract_motifs <- function(final_format, motif_args, verbose) {
  motifs <- try(do.call(final_format, args = motif_args), silent = TRUE)
  if (class(motifs) == "try-error") {
    if (verbose) cat("Parsing of motifs with '", final_format, "' failed!",
                     " Continuing autodetection..\n", sep = "") 
    return(NULL)
  } else return(motifs)

}
