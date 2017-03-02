######################################################################
## Benjamin Tremblay
##
## Import Homer-formatted motifs into R
##
######################################################################

#' @title [UNDER CONSTRUCTION] Load Homer motifs from a text file.
#'
#' @description
#' PLACEHOLDER TEXT.
#'
#' @param motif_file Character.
#' @param verbose Logical.
#' @param show_warnings Logical.
#' @param mot_length_cutoff Integer.
#' @param bkg_cutoff Integer.
#' @param target_cutoff Integer.
#' @param log_odds_thresh Double.
#' @param p_enrich_cutoff Double.
#' @param p_detect_cutoff Double.
#' @param tpos_cutoff Integer.
#' @param tstd_cutoff Double.
#' @param bpos_cutoff Integer.
#' @param bstd_cutoff Double.
#' @param strand_bias_cutoff Double.
#' @param multiplicity_cutoff Double.
#'
#' @return A list of motifs.
#'
#' @examples
#'   motifs <- system.file("extdata", "example.homer", package = "universalmotif")
#'   rmotifs <- read_homer(motifs)
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @export
read_homer <- function(motif_file, verbose = FALSE, show_warnings = TRUE,
                       mot_length_cutoff = NULL, bkg_cutoff = NULL,
                       target_cutoff = NULL, log_odds_thresh = NULL,
                       p_enrich_cutoff = NULL, p_detect_cutoff = NULL,
                       tpos_cutoff = NULL, tstd_cutoff = NULL,
                       bpos_cutoff = NULL, bstd_cutoff = NULL,
                       strand_bias_cutoff = NULL, multiplicity_cutoff = NULL) {

# log odds threshold is pretty important; can be generated manually w/ seq2profile.pl
# formula:
#
# Score for GGATGT
# score = log(pG1/0.25) + log(pG2/0.25) + log(pA3/0.25) + log(pT4/0.25) + log(pG5/0.25) + log(pT6/0.25)

  homer_raw <- readLines(motif_file)
  if (length(homer_raw) == 0) stop("could not read file, or file is empty")
  names(homer_raw) <- seq_along(homer_raw)

  # get motif info
  motif_info <- homer_raw[vapply(homer_raw, (function(x) grepl(">", x)),
                                 logical(1))]

  if (length(motif_info) == 0) stop("could not find any motifs")
  if (verbose) cat("Found", length(motif_info), "motifs.\n")

  # get motif matrix indices
  beg_mots <- vector(length = length(motif_info))
  end_mots <- vector(length = length(motif_info))
  for (i in seq_along(motif_info)) {
    beg_mots[i] <- as.integer(names(motif_info[i])) + 1
    end_mots[i] <- as.integer(names(motif_info[i + 1])) - 1
  }
  end_mots[length(motif_info)] <- length(homer_raw)

  # load motifs
  motifs <- mapply(hom_load, beg_mots, end_mots,
                   MoreArgs = list(homer_raw = homer_raw),
                   SIMPLIFY = FALSE)

  return(motifs)

}

######################################################################
######################################################################

hom_load <- function(beg_mot, end_mot, homer_raw) {
  as.matrix(read.table(textConnection(homer_raw[beg_mot:end_mot])))
}
