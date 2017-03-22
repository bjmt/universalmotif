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
#' @examples
#'    motifs <- system.file("extdata", "minimal.meme",
#'                          package = "universalmotif")
#'    motifs <- read_motifs(motifs)
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
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

  # motif format support and related functions
  available_formats <- list(

    # the sequence for these matters; functions will be called in this order
        
        # - meme is early on due to how unique the autodetection condition is
        # - transfac must be before homer since it often gets detected as homer

           c("MEME",       "read_meme"),
         c("Jaspar",       "read_jaspar"), 
       c("Transfac",       "read_transfac"),
          c("Homer",       "read_homer"),
       c("uniprobe",       "read_uniprobe")

    )

  # make life easier for later: 
  motif_args <- c(list("motif_file" = motif_file, "verbose" = verbose,
                       "out_class" = out_class), ...)

  # handling autodetection
  if (motif_format == "autodetect") {

    if (verbose) cat("Attempting to detect motif format..\n")
    motif_format <- find_format(motif_file, verbose, length(available_formats))
    final_format <- available_formats[vapply(available_formats, function(x)
                                      any(motif_format == x[1]), logical(1))]

    if (is.null(final_format)) {
      if (use_warning) {
        warning("autodetection failed", call. = FALSE)
        return(invisible(NULL))
      } else stop("autodetection failed", call. = FALSE)
    }

    # find_format can return a list of possibilities, which are tried
    # in the order they are listed in available_formats
    for (i in final_format) {

      if (verbose) cat("Detected as: ", i[1], "\n", 
                       "Using '", i[2], "' to parse motifs..\n",
                       sep = "")

      motifs <- extract_motifs(i[2], motif_args = motif_args, verbose = verbose)

      if (!is.null(motifs)) {
        if (verbose) cat("Done.\n")  
        return(motifs)
      }

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

find_format <- function(motif_file, verbose, available_formats_len) {

  # read file
  con <- file(motif_file)
  motifs_raw <- readLines(con)
  close(con)
  if (length(motifs_raw) == 0) stop("could not read file, or file is empty",
                                    call. = FALSE)
  # names(motifs_raw) <- seq_along(motifs_raw)

  final_format <- vector("character", available_formats_len)

  ##########
  # autodetection methods are pretty quick and dirty (at least for now)

  # 1 MEME
  # benchmarks: ~160 000 motifs takes a couple mins to parse

  if (any(grepl("MEME version", motifs_raw,
                fixed = TRUE))) final_format[1] <- "MEME"

  # 2 Jaspar

  if (any(grepl("[", motifs_raw, fixed = TRUE)) &&
      any(grepl("]", motifs_raw, fixed = TRUE))) final_format[2] <- "Jaspar" 

  # 3 Transfac

  if (any(grepl("^XX$", motifs_raw))) final_format[3] <- "Transfac"

  # 4 Homer

  arrow_test <- motifs_raw[vapply(motifs_raw, function(x) grepl("^>", x,),
                                  logical(1))]

  if (length(arrow_test) !=0) {
    if (all(!is.na(vapply(arrow_test, function(x)
                          strsplit(x, split = "\\s+")[[1]][3],
                          character(1))))) final_format[4] <- "Homer"
  }

  # 5 uniprobe

  if (any(grepl("A:", motifs_raw, fixed = TRUE)) &&
      any(grepl("C:", motifs_raw, fixed = TRUE)) &&
      any(grepl("G:", motifs_raw, fixed = TRUE)) &&
      any(grepl("T:", motifs_raw, fixed = TRUE))) final_format[5] <- "uniprobe"

  ##########

  final_format <- final_format[final_format != ""]

  if (length(final_format) == 0) return(NULL)

  return(final_format)

}

extract_motifs <- function(final_format, motif_args, verbose) {

  motifs <- tryCatch(

    do.call(final_format, args = motif_args),

    error = function(cond) {
      if (verbose) cat("Parsing of motifs with '", final_format, "' failed!\n",
                       "Continuing autodetection..\n", sep = "") 
      return(NULL)
    }

  )

  return(motifs)

}
