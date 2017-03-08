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
#'
#' @return A list of motifs.
#'
#' @examples
#'    motifs <- system.file("extdata", "example.jaspar",
#'                          package = "universalmotif")
#'    rmotifs <- read_jaspar(motifs)
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @include utils.R
#' @export
read_jaspar <- function(motif_file, verbose = FALSE,
                        mot_length_cutoff = NULL) {

  # check args
  check_logi_args(as.list(environment())[2])  # utils.R
  check_filter_args(as.list(environment())[3])  # utils.R

  # read file
  con <- file(motif_file)
  jaspar_raw <- readLines(con)
  close(con)
  if (length(jaspar_raw) == 0) stop("could not read file, or file is empty",
                                    call. = FALSE)
  names(jaspar_raw) <- seq_along(jaspar_raw)

  # get motif_info
  motif_info <- jaspar_raw[vapply(jaspar_raw, function(x) grepl(">", x),
                                  logical(1))]

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

    if (all(grepl("\\[", motif)) && all(grepl("\\]", motif))) {
      motif <- vapply(motif, function(x) sub("\\[", "", x), character(1))
      motif <- vapply(motif, function(x) sub("\\]", "", x), character(1))
    } else if (all(grepl("|", motif))) {
      motif <- vapply(motif, function(x) sub("|", replacement = "", x),
                      character(1))
    } else stop("motifs are not formatted as Jaspar; see `?read_jaspar`",
                call. = FALSE)

  }

  motif <- as.matrix(read.table(text = motif))
  rownames(motif) <- c("A", "C", "G", "T")
  colnames(motif) <- rep("", ncol(motif))
  
  return(motif)

}

jasp_filter <- function(motifs, mot_length_cutoff) {

  if (!is.null(mot_length_cutoff)) {
    motifs <- motifs[vapply(motifs, function(x) ncol(x) >= mot_length_cutoff,
                            logical(1))]
    if (length(motifs) == 0) stop("all motifs filtered out", call. = FALSE)
  }

  return(motifs)

}
