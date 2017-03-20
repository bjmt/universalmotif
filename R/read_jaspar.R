######################################################################
## Benjamin Tremblay
##
## Import Jaspar formatted motifs
##
######################################################################

#' @title [UNDER CONSTRUCTION] Load Jaspar motifs from a text file.
#'
#' @description
#' Motifs can be a mix of raw and Jaspar formatted motifs, as described
#' in \url{http://jaspar.genereg.net/html/TEMPLATES/help.html}.
#'
#' @param motif_file Character. Text file containing motifs.
#' @param verbose Logical.
#' @param mot_length_cutoff Integer.
#' @param out_class Character.
#'
#' @return A list of motifs.
#'
#' @examples
#'    motifs <- system.file("extdata", "example.jaspar",
#'                          package = "universalmotif")
#'    rmotifs <- read_jaspar(motifs)
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @include utils.R universalmotif-class.R universalmotif-methods.R
#' @export
read_jaspar <- function(motif_file, verbose = FALSE,
                        mot_length_cutoff = NULL, out_class = "matrix-2") {

  # check args
  check_logi_args(as.list(environment())[2])  # utils.R
  check_filter_args(as.list(environment())[3])  # utils.R
  check_out_class(out_class)

  # read file
  con <- file(motif_file)
  jaspar_raw <- readLines(con)
  close(con)
  if (length(jaspar_raw) == 0) stop("could not read file, or file is empty",
                                    call. = FALSE)
  names(jaspar_raw) <- seq_along(jaspar_raw)

  # get motif_info
  motif_info <- jaspar_raw[vapply(jaspar_raw, function(x)
                                  grepl(">", x, fixed = TRUE), logical(1))]

  # warning/verbose
  if (length(motif_info) == 0) stop("could not find any motifs", call. = FALSE)
  if (verbose) cat("Found", length(motif_info), "motifs.\n")

  # motif indices
  beg_mots <- as.integer(names(motif_info)) + 1
  end_mots <- c(beg_mots[2:length(beg_mots)] - 2, length(jaspar_raw))

  # load motifs
  motifs <- mapply(jasp_load, beg_mots, end_mots,
                   MoreArgs = list(jaspar_raw = jaspar_raw),
                   SIMPLIFY = FALSE)

  # parse motif names
  mot_names <- vapply(motif_info, function(x) strsplit(x, "\\s+")[[1]][2],
                      character(1))
  names(motifs) <- mot_names

  # filter
  motifs <- jasp_filter(motifs, mot_length_cutoff)

  if (is.null(motifs)) return(invisible(NULL))

  motifs <- mapply(jaspar_to_umot, motifs, mot_names, SIMPLIFY = FALSE)
  
  return(motifs)

}

######################################################################
######################################################################

jasp_load <- function(beg_mot, end_mot, jaspar_raw) {

  motif <- jaspar_raw[beg_mot:end_mot]
  motif <- motif[motif != ""]
  if (length(motif) != 4) stop("motifs must have 4 rows", call. = FALSE)

  if (all(mapply(grepl, c("A", "C", "G", "T"), motif))) {

    motif <- mapply(sub, x = motif, 
                    pattern = c("A", "C", "G", "T"),
                    MoreArgs = list(replacement = ""))

    if (all(grepl("[", motif, fixed = TRUE)) &&
        all(grepl("]", motif,  fixed = TRUE))) {
      motif <- vapply(motif, function(x) sub("[", "", x, fixed = TRUE),
                      character(1))
      motif <- vapply(motif, function(x) sub("]", "", x, fixed = TRUE),
                      character(1))
    } else if (all(grepl("|", motif, fixed = TRUE))) {
      motif <- vapply(motif, function(x) sub("|", replacement = "", x,
                                             fixed = TRUE), character(1))
    } else stop("motifs are not formatted as Jaspar; see `?read_jaspar`",
                call. = FALSE)

  }

  motif <- as.matrix(read.table(text = motif))
  rownames(motif) <- c("A", "C", "G", "T")
  colnames(motif) <- NULL
  
  return(motif)

}

jasp_filter <- function(motifs, mot_length_cutoff) {

  if (!is.null(mot_length_cutoff)) {
    motifs <- motifs[vapply(motifs, function(x) ncol(x) >= mot_length_cutoff,
                            logical(1))]
    if (length(motifs) == 0) {
      warning("All motifs were filtered out.\n", call. = FALSE)
      return(NULL)
    }
  }

  return(motifs)

}

jaspar_to_umot <- function(motif, name) {
  
  motif <- new("universalmotif", motif = motif, name = name, type = "PCM")
  
  return(motif)

}
