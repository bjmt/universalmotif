#######################################################################
#                          Benjamin Tremblay                          #
#                             2017-02-18                              #
#                  import meme format motifs into r                   #
#######################################################################

######################
#  public functions  #
######################

#-----------------------------------------------------------
read_meme <- function(motif_file, verbose = TRUE, show_warnings = TRUE,
                      use_alt_title = FALSE, mot_length_cutoff = NULL,
                      source_sites_cutoff = NULL, e_val_cutoff = NULL,
                      out_format = "by_col", motif_type = "pwm") {
  meme_raw <- readLines(motif_file)
  if (length(meme_raw) == 0) stop("Empty file.")
  meme_raw <- meme_raw[sapply(meme_raw, function(x) !identical(x, ""))]
  meme_raw <- meme_raw[sapply(meme_raw, function(x) !grepl("-----", x))]
  meme_raw <- meme_raw[sapply(meme_raw, function(x) !grepl("\\*\\*\\*\\*\\*",
                                                           x))]
  if (verbose) cat("\nReading MEME preamble:\n\n")
  version <- meme_ver(meme_raw, show_warnings)
  if (verbose) cat(paste0("\t", version[[1]], "\n"))
  alphabet <- meme_alph(meme_raw)
  if (verbose && !is.null(alphabet)) cat(paste0("\t", alphabet[[1]], "\n"))
  alph_type <- parse_alph(alphabet)
  background <- meme_bkg(meme_raw)
  if (verbose && !is.null(background)) cat("\tBackground letter frequencies\n")
  if (verbose && !is.null(background)) cat(paste0("\t", background[[1]],
                                                  "\n\n"))
  posmotifs <- pos_mots(meme_raw)
  if (nrow(posmotifs) == 0) stop("No motifs detected.")
  if (verbose) cat(paste0("Found ", nrow(posmotifs), " motif(s) of type: ",
                          alph_type[[1]], "\n\n"))
  motifs <- load_mots(meme_raw, posmotifs, show_warnings, use_alt_title,
                      out_format, alphabet, alph_type)
  return(motifs)
}
#-----------------------------------------------------------

########################
#  internal functions  #
########################

#-----------------------------------------------------------
meme_ver <- function(meme_raw, show_warnings) {
  for (i in seq_along(meme_raw)) {
    if (grepl("MEME version", meme_raw[i])) {
      ver <- strsplit(meme_raw[i], split = " ")[[1]][3]
      ver <- strsplit(ver, split = "\\.")[[1]][1]
      if (as.integer(ver) < 4 && show_warnings) {
        warning("MEME version less than 4 detected; please be careful.")
      }
      return(list(meme_raw[i], i))
    }
  }
  stop("No MEME version detected.")
}
#-----------------------------------------------------------
# TODO: strand info
#-----------------------------------------------------------
meme_alph <- function(meme_raw) {
  for (i in seq_along(meme_raw)) {
    if (grepl("ALPHABET=", meme_raw[i])) {
      alph <- strsplit(meme_raw[i], split = " ")[[1]][2]
      return(list(meme_raw[i], i, alph))
    }
  }
  return(NULL)
}
#-----------------------------------------------------------
parse_alph <- function(alphabet) {
  alph_string <- alphabet[[3]]
  alph_parsed <- strsplit(alph_string, split = "")[[1]]
  if (length(alph_parsed) == 4) {
    if ("A" %in% alph_parsed &&
        "C" %in% alph_parsed &&
        "G" %in% alph_parsed) {
      if ("T" %in% alph_parsed) return(list(alph = "DNA", len = 4,
                                            letters = alph_parsed))
      if ("U" %in% alph_parsed) return(list(alph = "RNA", len = 4,
                                            letters = alph_parsed))
      return(list(alph = "custom", len = nchar(alph_string),
                  letters = alph_parsed))
    }
  }
  warning("Non-standard alphabet detected; this may cause issues downstream.")
  return(list(alph = "custom", len = nchar(alph_string)))
}
#-----------------------------------------------------------
meme_bkg <- function(meme_raw) {
  for (i in seq_along(meme_raw)) {
    if (grepl("Background letter frequencies", meme_raw[i])) {
      return(list(meme_raw[i + 1], i))
    }
  }
  return(NULL)
}  # TODO: check that bkg frequencies add up to 1
#-----------------------------------------------------------
pos_mots <- function(meme_raw) {
  meme_raw2 <- data.frame(cbind(seq_along(meme_raw), meme_raw))
  meme_raw3 <- meme_raw2[sapply(meme_raw,
                                function(x) {
                                  any(strsplit(as.character(x),
                                      split = " ")[[1]][1] == "MOTIF")
                                }), ]
  return(meme_raw3)
}  # TODO: make this compatible with full MEME format (this is will take a lot
   # of code reworking)
#-----------------------------------------------------------
load_mots <- function(meme_raw, posmotifs, show_warnings, use_alt_title,
                      out_format, alphabet, alph_type) {
  allmots <- as.list(seq_len(nrow(posmotifs)))
  for (i in seq_len(nrow(posmotifs))) {
    if (i == nrow(posmotifs)) j <- length(meme_raw) else {
      j <- as.integer(row.names(posmotifs[i + 1, ])) - 1
    }  # TODO: get rid of this garbage row.names code
    if (grepl("letter-probability matrix:",
              meme_raw[as.integer(row.names(posmotifs[i, ])) + 1])) {
      lpm <- meme_raw[(as.integer(row.names(posmotifs[i, ])) + 2):j]
      lpm <- read.table(textConnection(lpm))
      test1 <- rowSums(lpm)
      if (all(test1 > 1.01) || all(test1 < 0.99) && show_warnings) {
        warning(paste(posmotifs[i, 2], "has positions which do not sum to 1."))
      }
      lpm <- format_mots(lpm, out_format, alphabet, alph_type)
      allmots[[i]] <- lpm
    }
  }
  if (!use_alt_title) {
    names(allmots) <- sapply(posmotifs[, 2],
                             function(x) strsplit(as.character(x),
                                                  split = " ")[[1]][2])
  } else if (use_alt_title) {
    names(allmots) <- sapply(posmotifs[, 2],
                             function(x) strsplit(as.character(x),
                                                  split = " ")[[1]][3])
  }  # add check for NA names?
  return(allmots)
}
#-----------------------------------------------------------
format_mots <- function(lpm, out_format, alphabet, alph_type) {
  if (out_format == "by_col") {
    lpm <- as.matrix(lpm)
    colnames(lpm) <- alph_type[[3]]
  }
  if (out_format == "by_row") {
    lpm <- as.matrix(lpm)
    colnames(lpm) <- alph_type[[3]]
    lpm <- t(lpm)
  }
  if (out_format == "seqLogo_pwm") {
    if (requireNamespace("seqLogo", quietly = TRUE)) {
      lpm <- seqLogo::makePWM(t(lpm))
    } else stop("seqLogo package is not installed.")
  }
    return(lpm)  # TODO: allow for RNA letters
}  # TODO: add more out format options
# modifying seqLogo pwm-class:
  # dimnames(lpm@pwm)[[1]]  <- c("A", "C", "G", "U")
  # lpm@alphabet <- "RNA"
  # also need to adjust @consensus
  # ..that doesn't actually seem to change the seqLogo graph output..
  # look into how the motifStack package is doing things
#-----------------------------------------------------------
# TODO: add filter function
