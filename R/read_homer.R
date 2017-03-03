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
#' @param bkg_cutoff Double. Percent. Motifs with more than num are cut.
#' @param target_cutoff Double. Percent. Motifs with more than num are kept.
#' @param log_odds_thresh Double.
#' @param p_enrich_cutoff Double.
#' @param p_detect_cutoff Double.
#' @param tpos_cutoff Double.
#' @param tstd_cutoff Double.
#' @param bpos_cutoff Double.
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
#' @include utils.R
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

  con <- file(motif_file)
  homer_raw <- readLines(con)
  close(con)
  if (length(homer_raw) == 0) stop("could not read file, or file is empty")
  names(homer_raw) <- seq_along(homer_raw)

  # get motif info
  motif_info <- homer_raw[vapply(homer_raw, function(x) grepl(">", x),
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

  # parse motif info
  info <- hom_info(motif_info)
  names(motifs) <- info[[1]]
  info_occ <- lapply(info[[4]], hom_occ)

  # warning checks
  if (show_warnings && any(is.na(info[[2]]))) {
    warning("motifs have missing log odds detection threshold values")
  }

  return(motifs)

}

######################################################################
######################################################################

hom_load <- function(beg_mot, end_mot, homer_raw) {
  con <- textConnection(homer_raw[beg_mot:end_mot])
  x <- as.matrix(read.table(con))
  close(con)
  colnames(x) <- c("A", "C", "G", "T")
  return(x)
}

hom_info <- function(motif_info) {
  info <- lapply(motif_info, function(x) {
                   con <- textConnection(x)
                   y <- scan(con, what = "", quiet = TRUE)
                   close(con)
                   return(y)})
  mot_names <- vapply(info, function(x) x[2], character(1))
  logodds <- vapply(info, function(x) as.double(x[3]), double(1))
  penr <- vapply(info, function(x) as.double(x[4]), double(1))
  occinfo <- vapply(info, function(x) x[6], character(1))
  motstats <- vapply(info, function(x) x[7], character(1))
  info <- list(mot_names, logodds, penr, occinfo, motstats)
  return(info)
}

# T:17311.0(44.36%),B:2181.5(5.80%),P:1e-10317
hom_occ <- function(occinfo) {
  occ <- strsplit(occinfo, split = ",")[[1]]
  tnum <- occ[vapply(occ, function(x) grepl("T", x), logical(1))]
  bnum <- occ[vapply(occ, function(x) grepl("B", x), logical(1))]
  pnum <- occ[vapply(occ, function(x) grepl("P", x), logical(1))]
  occ <- c("T" = tnum, "B" = bnum, "P" = pnum)

  occ[c("T", "B")] <- vapply(occ[c("T", "B")], function(x) {
                   x <- strsplit(x, split = "\\(")[[1]][2]
                   x <- strsplit(x, split = "%")[[1]][1]
                   return(x)}, character(1))
  occ["P"] <- strsplit(occ[3], split = ":")[[1]][2]
  occ <- vapply(occ, as.double, double(1))
  return(occ)
}
