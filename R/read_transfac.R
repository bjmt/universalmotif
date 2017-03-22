######################################################################
## Benjamin Tremblay
##
## Read TRANSFAC motifs from a text file.
##
######################################################################

#' [UNDER CONSTRUCTION] Load TRANSFAC motifs from a text file.
#'
#' Support for TRANSFAC and TRANSFAC-like motifs. Each motif must be terminated
#' with a line containing 'XX' or '//'.
#'
#' @param motif_file Character.
#' @param verbose Logical.
#' @param mot_length_cutoff Integer.
#' @param out_class Character.
#'
#' @return A list of motif objects of the specified class.
#'
#' @examples
#'    motifs <- system.file("extdata", "transfac.txt",
#'                          package = "universalmotif")
#'    motifs <- read_transfac(motifs)
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @include utils.R universalmotif-class.R universalmotif-methods.R
#' @export
read_transfac <- function(motif_file, verbose = FALSE,
                          mot_length_cutoff = NULL, out_class = "matrix-2") {
  
  # check args
  check_logi_args(as.list(environment())[2])  # utils.R
  check_filter_args(as.list(environment())[3])  # utils.R
  check_out_class(out_class)

  # read file
  transfac_raw <- readLines(con <- file(motif_file)); close(con)
  transfac_raw <- transfac_raw[transfac_raw != ""]
  if (length(transfac_raw) == 0) stop("could not read file, or file is empty",
                                      call. = FALSE)
  names(transfac_raw) <- seq_along(transfac_raw)

  # get motif indices with grepl for speed
  xx_transfac <- transfac_raw[vapply(transfac_raw, function(x) grepl("^XX$", x),
                                     logical(1))]
  # be safe:
  # xx_transfac <- xx_transfac[vapply(xx_transfac, function(x)
  #                                   strsplit(x, split = "\\s+")[[1]][1] == "XX",
  #                                   logical(1))]
  fs_transfac <- transfac_raw[vapply(transfac_raw, function(x) grepl("//", x),
                                     logical(1))]
  end_mots <- transfac_indices(xx_transfac, fs_transfac)
  if (length(end_mots) == 0) stop("could not find any motifs", call. = FALSE)
  if (verbose) cat("Found", length(end_mots), "motifs.\n")
  beg_mots <- c(1, end_mots[-length(end_mots)] + 1)

  mot_names <- mapply(transfac_name, beg_mots, end_mots,
                      MoreArgs = list(transfac_raw = transfac_raw),
                      SIMPLIFY = FALSE)

  motifs <- mapply(transfac_load, beg_mots, end_mots,
                   MoreArgs = list(transfac_raw = transfac_raw),
                   SIMPLIFY = FALSE)

  mot_names <- mapply(function(x, y) ifelse(is.null(y), x, y),
                      seq_along(mot_names), mot_names)
  names(motifs) <- mot_names

  motifs <- mapply(transfac_to_umot, motifs, mot_names, SIMPLIFY = FALSE)

  return(motifs)

}

######################################################################
######################################################################

transfac_indices <- function(xx_transfac, fs_transfac) {

  xx_t <- as.integer(names(xx_transfac))
  names(xx_t) <- rep("xx_t", length(xx_t))
  fs_t <- as.integer(names(fs_transfac)) - 1
  names(fs_t) <- rep("fs_t", length(fs_t))

  final_t <- c(xx_t, fs_t)
  final_t <- final_t[!duplicated(final_t)]
  final_t <- ifelse(names(final_t) == "fs_t", final_t + 1, final_t)

  return(sort(final_t))

}

transfac_name <- function(beg_mot, end_mot, transfac_raw) {

  motif <- transfac_raw[beg_mot:end_mot]

  # get title, if any
  mot_name <- motif[vapply(motif, function(x)
                                  strsplit(x, split = "\\s+")[[1]][1] == "NA",
                                  logical(1))]
  if (length(mot_name) == 0) {
    mot_name <- motif[vapply(motif, function(x)
                                    strsplit(x, split = "\\s+")[[1]][1] == "DE",
                                    logical(1))]
  }
  if (length(mot_name) == 0) {
    mot_name <- motif[vapply(motif, function(x)
                             strsplit(x, split = "\\s+")[[1]][1] == "ID",
                             logical(1))]
  }

  if (length(mot_name) == 0) return(NULL) else {
    return(strsplit(mot_name, split = "\\s+")[[1]][2])
  }

}

transfac_load <- function(beg_mot, end_mot, transfac_raw) {

  motif <- transfac_raw[beg_mot:end_mot]

  # get letters, if present
  transfac_letters <- motif[vapply(motif, function(x)
                                   strsplit(x, split = "\\s+")[[1]][1] == "P0",
                                   logical(1))]

  beg_mot01 <- which(vapply(motif, function(x)
                            strsplit(x, split = "\\s+")[[1]][1] == "01",
                            logical(1)))
  if (length(beg_mot01) == 0) {
    beg_mot01 <- which(vapply(motif, function(x)
                              strsplit(x, split = "\\s+")[[1]][1] == "1",
                              logical(1)))
  }
  
  if (length(transfac_letters) != 0) {
    transfac_letters <- strsplit(transfac_letters, split = "\\s+")[[1]]
    transfac_letters <- transfac_letters[-1]
    motif <- read.table(text = motif[beg_mot01:(length(motif) - 1)],
                        row.names = 1)
    motif <- motif[, 1:length(transfac_letters)]
    motif <- as.matrix(motif)
    colnames(motif) <- transfac_letters
  } else {
    motif <- read.table(text = motif[beg_mot01:(length(motif) - 1)],
                        row.names = 1)
    motif <- motif[, 1:4]
    motif <- as.matrix(motif)
    colnames(motif) <- c("A", "C", "G", "T")
  }

  rownames(motif) <- NULL
  return(motif)

}

transfac_to_umot <- function(motif, name) {

  motif <- new("universalmotif", motif = t(motif), name = name, type = "PCM")

}
