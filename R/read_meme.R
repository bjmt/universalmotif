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
#' @param motif_file Character. Text file containing MEME format motifs.
#' @param verbose Logical. \code{read_meme} will return information regarding any 
#'   detected information in the file preamble.
#' @param show_warnings Logical. Set this option to FALSE if you using this function
#'   in a loop and are sure the warnings can be ignored.
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
#'     \item \code{matrix-1} Matrix with alphabet letters as rows.
#'     \item \code{matrix-2} Matrix with alphabet letters as columns.
#'   }
#' @param motif_type Character. MEME motifs can be provided as either letter-
#'   probability matrices or log-odds matrices. The load the former, set this
#'   option to \code{lpm}; for the latter, \code{lom}.
#'
#' @return A list of motif objects of the specified class.
#'
#' @examples
#'   motifs <- system.file("extdata", "minimal.meme", package = "universalmotif")
#'   rmotifs <- read_meme(motifs, show_warnings = FALSE, out_class = "matrix-2")
#'
#' @author Benjamin Tremblay <b2trembl@uwaterloo.ca>
#' @export
read_meme <- function(motif_file, verbose = FALSE, show_warnings = TRUE,
                      use_alt_title = FALSE, mot_length_cutoff = NULL,
                      source_sites_cutoff = NULL, e_val_cutoff = NULL,
                      out_class = "matrix-2", motif_type = "lpm") {

  # read file
  meme_raw <- readLines(motif_file)
  if (length(meme_raw) == 0) stop("Could not read file, or file is empty.")
  names(meme_raw) <- seq_along(meme_raw)

  # inspect preamble
  version <- meme_ver(meme_raw, show_warnings)
  alphabet <- meme_alph(meme_raw)
  alph_type <- parse_alph(alphabet, show_warnings)
  strands <- meme_strands(meme_raw, show_warnings)
  background <- meme_bkg(meme_raw, alph_type, show_warnings)

  # verbose call
  if (verbose) {
    cat("\nReading MEME preamble:\n\n")
    cat(paste0("\t", version, "\n"))
    if (!is.null(alphabet)) cat(paste0("\t", alphabet[[1]], "\n"))
    if (!is.null(strands)) cat(paste0("\tstrands:", strands, "\n"))
    if (!is.null(background)) cat("\tBackground letter frequencies\n")
    if (!is.null(background)) cat(paste0("\t", background[[1]], "\n\n"))
  }

  # get motif names 
  mot_names <- get_names(meme_raw, use_alt_title)
  mot_names <- mot_names[!is.na(mot_names)]  # to be honest I don't quite get 
                                             # why there are NA values here

  # get motif matrix indices
  posmotifs <- pos_mots(meme_raw, motif_type)
  nposmotifs <- as.list(names(posmotifs))

  if (length(posmotifs) == 0) stop("No motifs detected.")
  
  # get motif info
  info_mots <- lapply(posmotifs, get_info, motif_type =  motif_type)

  # detect and read motifs
  motifs <- mapply(load_mots, nposmotifs, info_mots,
                   MoreArgs = list(meme_raw = meme_raw,
                                   show_warnings = show_warnings),
                   SIMPLIFY = FALSE)

  # another verbose call
  if (verbose) cat(paste("Found", length(posmotifs), "motif(s) of type:",
                          alph_type[[1]], "\n\n"))

  if (!identical(length(mot_names), length(motifs))) {
    stop("The number of motif names do not match the number of motif matrices.")
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
                          source_sites_cutoff, e_val_cutoff)
  }

  return(motifs)

}

######################################################################

# preamble functions: uses for loops to not waste time looking at every line

#' @keywords internal
meme_ver <- function(meme_raw, show_warnings) {
  for (i in seq_along(meme_raw)) {
    if (grepl("MEME version", meme_raw[i])) {
      ver <- strsplit(meme_raw[i], split = "\\s+")[[1]][3]
      ver <- strsplit(ver, split = "\\.")[[1]][1]
      if (as.integer(ver) < 4 && show_warnings) {
        warning("MEME version less than 4 detected; this may cause parsing issues.")
      }
      return(meme_raw[i])
    }
  }
  stop("No MEME version detected.")
}

#' @keywords internal
meme_alph <- function(meme_raw) {
  for (i in seq_along(meme_raw)) {
    if (grepl("ALPHABET=", meme_raw[i])) {
      alph <- strsplit(meme_raw[i], split = "\\s+")[[1]][2]
      return(list(meme_raw[i], i, alph))
    }
  }
  return(NULL)
}

#' @keywords internal
parse_alph <- function(alphabet, show_warnings) {
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
  if (show_warnings) {
    warning("Non-standard alphabet detected; this may cause issues with other packages.")
  }
  return(list(alph = "custom", len = nchar(alph_string), letters = alph_parsed))
}

#' @keywords internal
meme_strands <- function(meme_raw, show_warnings) {
  for (i in seq_along(meme_raw)) {
    if (grepl("strands:", meme_raw[i])) {
      strands <- strsplit(meme_raw[i], split = "\\s+")[[1]]
      if (length(strands) == 3 && all(c("+", "-") %in% strands)) return(strands[2:3])
      if (length(strands) == 2 && any(c("+", "-") %in% strands)) return(strands[2])
      if (show_warnings) warning("Could not parse strand information.")
    }
  }
  return(NULL)
}

#' @keywords internal
meme_bkg <- function(meme_raw, alph_type, show_warnings) {
  for (i in seq_along(meme_raw)) {
    if (grepl("Background letter frequencies", meme_raw[i])) {
      if (show_warnings) {
        frequencies <- strsplit(meme_raw[i + 1], split = "\\s+")[[1]]
        frequencies <- frequencies[seq(from = 2, to = alph_type[[2]], by = 2)]
        f_test <- sum(as.double(frequencies))
        if (f_test < 0.99 || f_test > 1.01) {
          warning("Background letter frequencies do not add up to 1.")
        }
      }
      return(list(meme_raw[i + 1], i))
    }
  }
  return(NULL)
}

######################################################################

# get motif names

#' @keywords internal
get_names <- function(meme_raw, use_alt_title) {
  meme_raw <- meme_raw[sapply(meme_raw, function(x) {
                              all(strsplit(x,
                                  split = " ")[[1]][1] == "MOTIF")})]
  if (!use_alt_title) {
    return(sapply(meme_raw, function(x) strsplit(x, split = "\\s+")[[1]][2]))
  } else {
    return(sapply(meme_raw, function(x) strsplit(x, split = "\\s+")[[1]][3]))
  }
  stop("No motif names detected.")
}

# find motifs

#' @keywords internal
pos_mots <- function(meme_raw, motif_type) {
  if (motif_type == "lpm") mtype <- c("letter-probability", "matrix:")
  if (motif_type == "lom") mtype <- c("log-odds", "matrix:")
  if (!exists("mtype")) stop("motif_type must be either lpm or lom")
  meme_raw <- meme_raw[sapply(meme_raw,function(x) {
                              all(strsplit(as.character(x),
                              split = "\\s+")[[1]][1:2] == mtype)})]
  return(meme_raw[!is.na(meme_raw)])
}

# load motif info (such as source sites, etc)

#' @keywords internal
get_info <- function(posmotifs, motif_type) {
  posmotifs <- strsplit(posmotifs, split = "\\s+")[[1]]
  # alength=, w=, nsites, E=
  names(posmotifs) <- seq_along(posmotifs)
  alength <- posmotifs[sapply(posmotifs, function(x) identical(x, "alength="))]
  if (length(alength) > 0) alength <- posmotifs[as.integer(names(alength)) + 1]
  mlength <- posmotifs[sapply(posmotifs, function(x) identical(x, "w="))]
  if (length(mlength) > 0) mlength <- posmotifs[as.integer(names(mlength)) + 1]
  if (motif_type == "lpm") nsites <- "nsites=" else nsites <- "n="
  nsites <- posmotifs[sapply(posmotifs, function(x) identical(x, nsites))]
  if (length(nsites) > 0) nsites <- posmotifs[as.integer(names(nsites)) + 1]
  e_val <- posmotifs[sapply(posmotifs, function(x) identical(x, "E="))]
  if (length(nsites) > 0) e_val <- posmotifs[as.integer(names(e_val)) + 1]
  return(c("alength" = alength, "w" = mlength, "nsites" = nsites, "E" = e_val))
}

# load motifs as matrices

#' @keywords internal
load_mots <- function(posmotifs, info_mots, meme_raw, show_warnings) {
  posmot <- as.integer(posmotifs)
  if (!is.na(info_mots["w.6"])) {
    mot <- meme_raw[(posmot + 1):(posmot + as.integer(info_mots["w.6"]))]  # ..seems dangerous
  } else {
    mot <- vector(length = (length(meme_raw) - posmot))
    mot[seq_along(mot)] <- NA
    for (i in seq_along(mot)) {
      if (meme_raw[posmot + i] == "") break
      mot[i] <- meme_raw[posmot + i]
    }
    mot <- mot[!is.na(mot)]
  }
  mot_mat <- read.table(textConnection(mot))
  if (!is.na(info_mots["alength.4"])) {
    if (ncol(mot_mat) != as.integer(info_mots["alength.4"]) && show_warnings) {
      warning("Motif alength and actual number of columns do not match.")
    }
  }
  return(as.matrix(mot_mat))
}

######################################################################

# filter out unwanted motifs

#' @keywords internal
filter_meme <- function(motifs, info_mots, mot_length_cutoff,
                        source_sites_cutoff, e_val_cutoff) {
  if (length(motifs) != length(info_mots)) {
    stop("Please check that the file is properly formatted.")
  }
  index <- names(info_mots)
  names(index) <- seq_along(info_mots)
  if (!is.null(mot_length_cutoff)) {
    info_mots <- info_mots[sapply(info_mots, function(x) {
                                    as.integer(x["w.6"]) > mot_length_cutoff})]
  }
  if (!is.null(source_sites_cutoff)) {
    info_mots <- info_mots[sapply(info_mots, function(x) {
                                    as.integer(x["nsites.8"]) > source_sites_cutoff})]
  }
  if (!is.null(e_val_cutoff)) {
    info_mots <- info_mots[sapply(info_mots, function(x) {
                                    as.double(x["E.10"]) < e_val_cutoff})]
  }
  index <- index[sapply(index, function(x) any(x == names(info_mots)))]
  index <- as.integer(names(index))
  return(motifs[index])
}


######################################################################

# final conversion to desired class; this is done by convert_motifs.R
# (with the exception of 'matrix-1' and 'matrix-2')

#' @keywords internal
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
