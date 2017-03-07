######################################################################
## Benjamin Tremblay
##
## Read MEME motifs from a file connection into a list
##
######################################################################

#' @title [UNDER CONSTRUCTION] Load MEME motifs from a text file.
#'
#' @description
#' \code{read_meme} can load MEME motifs from a file connection to a list of
#' motifs in the class specified by the user. Both full MEME format and 
#' minimal MEME format files can be read; the latter requires at least a MEME
#' version in the preamble. \code{read_meme} was written for MEME version 4 and
#' above; reading a file with a lower version of MEME will throw a warning, as
#' this may result in incorrect parsing.
#'
#' NOTE: parsing full meme format files is quite slower than minimal meme.
#'
#' @param motif_file Character. Text file containing MEME format motifs.
#' @param verbose Logical. \code{read_meme} will return information regarding any 
#'   detected information in the file preamble.
#' @param use_alt_title Logical. Set this to TRUE if you wish to use the
#'   alternate name in the resulting motif object.
#' @param mot_length_cutoff Integer. Setting this option will cause \code{
#'   read_meme} to skip over motifs which are shorter than specified.
#' @param source_sites_cutoff Integer. Setting this option will cause \code{
#'   read_meme} to skip over motifs which are found in fewer numbers than
#'   requested.
#' @param e_val_cutoff Double. \code{read_meme} will skip over motifs with an
#'   E value higher than specified.
#' @param out_class Character. Specify the class the resulting motif objects
#'   will be stored as. These classes include:
#'   \itemize{
#'     \item \code{matrix-1}: Matrix with alphabet letters as rows.
#'     \item \code{matrix-2}: Matrix with alphabet letters as columns.
#'   }
#' @param motif_type Character. MEME motifs can be provided as either letter-
#'   probability matrices or log-odds matrices. The load the former, set this
#'   option to \code{lpm}; for the latter, \code{lom}. 
#'
#' @return A list of motif objects of the specified class.
#'
#' @examples
#'   motifs <- system.file("extdata", "minimal.meme", package = "universalmotif")
#'   rmotifs <- read_meme(motifs, out_class = "matrix-2")
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @include utils.R
#' @export
read_meme <- function(motif_file, verbose = FALSE,
                      use_alt_title = FALSE, mot_length_cutoff = NULL,
                      source_sites_cutoff = NULL, e_val_cutoff = NULL,
                      out_class = "matrix-2", motif_type = "lpm") {

  # check args
  check_logi_args(as.list(environment())[2:3])  # utils.R
  check_filter_args(as.list(environment())[4:6])  # utils.R
  if (!motif_type %in% c("lpm", "lom")) {
    stop("'motif_type' must be \"lpm\" or \"lom\"",
         call. = FALSE)
  }
  if (!out_class %in% c("matrix-1", "matrix-2")) {
    stop("please see `?read_meme` for a list of available 'out_class' options",
         call. = FALSE)
  }

  # read file
  con <- file(motif_file)
  meme_raw <- readLines(con)
  close(con)
  if (length(meme_raw) == 0) stop("could not read file, or file is empty",
                                  call. = FALSE)
  names(meme_raw) <- seq_along(meme_raw)

  # inspect preamble
  version <- meme_ver(meme_raw)
  alphabet <- meme_alph(meme_raw)
  alph_type <- parse_alph(alphabet)
  strands <- meme_strands(meme_raw)
  background <- meme_bkg(meme_raw, alph_type)

  # verbose call
  if (verbose) {
    cat("\nReading MEME preamble:\n\n")
    cat(paste0("\t", version, "\n"))
    if (!is.null(alphabet)) cat("\t", alphabet[[1]], "\n", sep = "")
    if (!is.null(strands)) cat("\tstrands:", strands, "\n")
    if (!is.null(background)) cat("\tBackground letter frequencies\n")
    if (!is.null(background)) cat("\t", background[[1]], "\n\n", sep = "")
  }

  # get motif names 
  mot_names <- get_names(meme_raw, use_alt_title)
  mot_names <- mot_names[!is.na(mot_names)]  # to be honest I don't quite get 
                                             # why there are NA values here

  # get motif matrix indices
  posmotifs <- pos_mots(meme_raw, motif_type)
  nposmotifs <- as.list(names(posmotifs))

  if (length(posmotifs) == 0) stop("no motifs detected", call. = FALSE)
  
  # get motif info
  info_mots <- lapply(posmotifs, get_info, motif_type =  motif_type)

  # detect and read motifs
  motifs <- mapply(load_mots, nposmotifs, info_mots,
                   MoreArgs = list(meme_raw = meme_raw),
                   SIMPLIFY = FALSE)

  # another verbose call
  if (verbose) cat("Found", length(posmotifs), "motif(s) of type:",
                          alph_type[[1]], "\n\n")

  if (!identical(length(mot_names), length(motifs))) {
    stop("the number of motif names do not match the number of motif matrices",
         call. = FALSE)
  }

  # convert motifs to desired class
  motifs <- mapply(convert_mots, mot_names, info_mots, motifs,
                   MoreArgs = list(out_class = out_class,
                                   alph_type = alph_type),
                   SIMPLIFY = FALSE)

  names(motifs) <- mot_names

  # final step: get rid of motifs which do not match filter options
  if (!is.null(c(mot_length_cutoff, source_sites_cutoff, e_val_cutoff))) {
    motifs <- filter_meme(motifs, info_mots, mot_length_cutoff,
                          source_sites_cutoff, e_val_cutoff, verbose)
  }

  if (length(motifs) == 0) return(NULL)
  return(motifs)

}

######################################################################
######################################################################

# preamble functions: uses for loops to not waste time looking at every line

meme_ver <- function(meme_raw) {
  for (i in seq_along(meme_raw)) {
    if (grepl("MEME version", meme_raw[i])) {
      ver <- strsplit(meme_raw[i], split = "\\s+")[[1]][3]
      ver <- strsplit(ver, split = "\\.")[[1]][1]
      if (as.integer(ver) < 4) {
        warning("MEME version less than 4 detected; this may cause parsing issues",
                call. = FALSE)
      }
      return(meme_raw[i])
    }
  }
  stop("no MEME version detected", call. = FALSE)
}

meme_alph <- function(meme_raw) {
  for (i in seq_along(meme_raw)) {
    if (grepl("ALPHABET=", meme_raw[i])) {
      alph <- strsplit(meme_raw[i], split = "\\s+")[[1]][2]
      return(list(meme_raw[i], i, alph))
    }
  }
  return(NULL)
}

parse_alph <- function(alphabet) {
  alph_string <- alphabet[[3]]
  alph_parsed <- strsplit(alph_string, split = "")[[1]]
  if (length(alph_parsed) == 4) {
    if (all(c("A", "C", "G") %in% alph_parsed)) {
      if ("T" %in% alph_parsed) return(list(alph = "DNA", len = 4,
                                            letters = alph_parsed))
      if ("U" %in% alph_parsed) return(list(alph = "RNA", len = 4,
                                            letters = alph_parsed))
      return(list(alph = "custom", len = 4, letters = alph_parsed))
    }  # TODO: add support for amino acid alphabets
  }
  warning("non-standard alphabet detected; this may cause issues with other packages",
          call. = FALSE)
  return(list(alph = "custom", len = nchar(alph_string), letters = alph_parsed))
}

meme_strands <- function(meme_raw) {
  for (i in seq_along(meme_raw)) {
    if (grepl("strands:", meme_raw[i])) {
      strands <- strsplit(meme_raw[i], split = "\\s+")[[1]]
      if (length(strands) == 3 && all(c("+", "-") %in% strands)) return(strands[2:3])
      if (length(strands) == 2 && any(c("+", "-") %in% strands)) return(strands[2])
      warning("could not parse strand information", call. = FALSE)
    }
  }
  return(NULL)
}

meme_bkg <- function(meme_raw, alph_type) {
  for (i in seq_along(meme_raw)) {
    if (grepl("Background letter frequencies", meme_raw[i])) {
      frequencies <- strsplit(meme_raw[i + 1], split = "\\s+")[[1]]
      frequencies <- frequencies[seq(from = 2, to = (alph_type[[2]] * 2), by = 2)]
      f_test <- sum(as.double(frequencies))
      if (f_test < 0.99 || f_test > 1.01) {
        warning("background letter frequencies do not add up to 1",
                call. = FALSE)
      }
      return(list(meme_raw[i + 1], i))
    }
  }
  return(NULL)
}

######################################################################
######################################################################

# get motif names

get_names <- function(meme_raw, use_alt_title) {
  meme_raw <- meme_raw[vapply(meme_raw, function(x) {
                              all(strsplit(x,
                                  split = "\\s+")[[1]][1] == "MOTIF")},
                              logical(1))]
  if (!use_alt_title) {
    return(vapply(meme_raw, function(x) strsplit(x, split = "\\s+")[[1]][2],
                  character(1)))
  } else {
    return(vapply(meme_raw, function(x) strsplit(x, split = "\\s+")[[1]][3],
                  character(1)))
  }
  stop("no motif names detected", call. = FALSE)
}

# find motifs

pos_mots <- function(meme_raw, motif_type) {
  if (motif_type == "lpm") mtype <- c("letter-probability", "matrix:")
  if (motif_type == "lom") mtype <- c("log-odds", "matrix:")
  if (!exists("mtype")) stop("'motif_type' must be either \"lpm\" or \"lom\"",
                             call. = FALSE)
  meme_raw <- meme_raw[vapply(meme_raw, function(x) {
                              all(strsplit(as.character(x),
                              split = "\\s+")[[1]][1:2] == mtype)},
                              logical(1))]
  return(meme_raw[!is.na(meme_raw)])
}

# load motif info (such as source sites, etc)

get_info <- function(posmotifs, motif_type) {
  posmotifs <- strsplit(posmotifs, split = "\\s+")[[1]]
  # alength=, w=, nsites, E=
  names(posmotifs) <- seq_along(posmotifs)
  alength <- posmotifs[vapply(posmotifs, function(x) identical(x, "alength="),
                              logical(1))]
  if (length(alength) > 0) alength <- posmotifs[as.integer(names(alength)) + 1]
  mlength <- posmotifs[vapply(posmotifs, function(x) identical(x, "w="),
                              logical(1))]
  if (length(mlength) > 0) mlength <- posmotifs[as.integer(names(mlength)) + 1]
  if (motif_type == "lpm") nsites <- "nsites=" else nsites <- "n="
  nsites <- posmotifs[vapply(posmotifs, function(x) identical(x, nsites),
                             logical(1))]
  if (length(nsites) > 0) nsites <- posmotifs[as.integer(names(nsites)) + 1]
  e_val <- posmotifs[vapply(posmotifs, function(x) identical(x, "E="),
                            logical(1))]
  if (length(nsites) > 0) e_val <- posmotifs[as.integer(names(e_val)) + 1]
  return(c("alength" = alength, "w" = mlength, "nsites" = nsites, "E" = e_val))
}

# load motifs as matrices

load_mots <- function(posmotifs, info_mots, meme_raw) {
  posmot <- as.integer(posmotifs)
  if (!is.na(info_mots["w.6"])) {
    mot <- meme_raw[(posmot + 1):(posmot + as.integer(info_mots["w.6"]))]
  } else {
    mot <- vector(length = (length(meme_raw) - posmot))  # must be a better way
    mot[seq_along(mot)] <- NA
    for (i in seq_along(mot)) {
      if (meme_raw[posmot + i] == "") break
      mot[i] <- meme_raw[posmot + i]
    }
    mot <- mot[!is.na(mot)]
  }
  mot_mat <- read.table(text = mot)
  if (!is.na(info_mots["alength.4"])) {
    if (ncol(mot_mat) != as.integer(info_mots["alength.4"])) {
      warning("motif 'alength' and actual number of columns do not match",
              call. = FALSE)
    }
  }
  return(as.matrix(mot_mat))
}

######################################################################
######################################################################

# filter out unwanted motifs

filter_meme <- function(motifs, info_mots, mot_length_cutoff,
                        source_sites_cutoff, e_val_cutoff, verbose) {
  if (length(motifs) != length(info_mots)) {
    stop("please check that the file is properly formatted", call. = FALSE)
  }
  index <- names(info_mots)
  names(index) <- seq_along(info_mots)
  if (!is.null(mot_length_cutoff)) {
    info_mots <- info_mots[vapply(info_mots, function(x) {
                                    as.integer(x["w.6"]) >= mot_length_cutoff},
                                  logical(1))]
  }
  if (!is.null(source_sites_cutoff)) {
    info_mots <- info_mots[vapply(info_mots, function(x) {
                                    as.integer(x["nsites.8"]) >= source_sites_cutoff},
                                  logical(1))]
  }
  if (!is.null(e_val_cutoff)) {
    info_mots <- info_mots[vapply(info_mots, function(x) {
                                    as.double(x["E.10"]) <= e_val_cutoff},
                                  logical(1))]
  }
  index <- index[vapply(index, function(x) any(x == names(info_mots)),
                        logical(1))]
  index <- as.integer(names(index))
  if (verbose) cat("Number of filtered out motifs:",
                   (length(motifs) - length(index)), "\n\n")
  return(motifs[index])
}


######################################################################
######################################################################

# final conversion to desired class; this is done by convert_motifs.R
# (with the exception of 'matrix-1' and 'matrix-2')

convert_mots <- function(mot_names, info_mots, motifs, out_class,
                         alph_type) {
  if (out_class == "matrix-2") {
    colnames(motifs) <- alph_type[[3]]
    return(motifs)
  }
  if (out_class == "matrix-1") {
    motifs <- t(motifs)
    rownames(motifs) <- alph_type[[3]]
    return(motifs)
  }
}
