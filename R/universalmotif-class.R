######################################################################
## Benjamin Tremblay
##
## S4 class for storing motifs
##
######################################################################

#' [UNDER CONSTRUCTION] universalmotif: Motif class.
#'
#' Placeholder description.
#'
#' @slot name Character. Motif name.
#' @slot motif Matrix. Contains motif.
#' @slot alphabet Character. Describes alphabet. Can be: DNA, RNA, AA, custom.
#' @slot type Character. Can be: PCM, PPM, or PWM.
#' @slot icscore Numeric. Total information content for all positions.
#' @slot nsites Numeric. Total number of sites containing motif.
#' @slot pseudoweight Numeric. Amount of smoothing to apply.
#' @slot bkg Numeric. Background letter frequencies.
#' @slot consensus Character. Motif consensus sequence.
#' @slot strand Character. '+' or '-'.
#' @slot pval Numeric.
#' @slot eval Numeric.
#' @slot extrachar Character.
#' @slot extranum Numeric.
#' @slot bkgsites Numeric.
#'
#' @return A motif object.
#'
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
#' @export
universalmotif <- setClass("universalmotif",
                           slots = list(name = "character", motif = "matrix",
                                        alphabet = "character",
                                        type = "character", icscore = "numeric",
                                        nsites = "numeric", 
                                        pseudoweight = "numeric",
                                        bkg = "numeric", 
                                        consensus = "character",
                                        strand = "character", pval = "numeric",
                                        extrachar = "character", 
                                        extranum = "numeric",
                                        eval = "numeric", bkgsites = "numeric"))

setValidity("universalmotif",
            function(object) {

              msg <- NULL
              valid <- TRUE

              if (length(motif_slots(object, "name")) != 1) {
                valid <- FALSE
                msg <- c(msg, "motif name must be a character vector of length 1")
              }

              if (!motif_slots(object, "alphabet") %in% c("DNA", "RNA", "AA",
                                                          "custom") ||
                  length(motif_slots(object, "alphabet")) != 1) {
                valid <- FALSE
                msg <- c(msg, "motif alphabet must be DNA, RNA, AA, or custom")
              }

              if (!motif_slots(object, "type") %in% c("PCM", "PPM", "PWM",
                                                      "ICM") ||
                  length(motif_slots(object, "type")) != 1) {
                valid <- FALSE
                msg <- c(msg, "motif type must be PCM, PPM, PWM or ICM")
              }

              # if (motif_slots(object, "type") == "PCM") {
              #   nsitestest <- unique(colSums(object@motif))
              #   if (nsitestest != motif_slots(object, "nsites")) {
              #     valid <- FALSE
              #     msg <- c(msg, "PCM colSums and motif nsites must be equal")
              #   }
              #   if (any(motif_slots(object, "motif") < 0)) {
              #     valid <- FALSE
              #     msg <- c(msg, "PCM values must be positive")
              #   }
                # if (!any(all.equal(colSums(motif_slots(object, "motif"))), 1,
                #          tolerance = 0.01)) {
                #   valid <- FALSE
                #   msg <- c(msg, "PCM colSums must be equal to 1")
                # }
              # }

              if (!any(motif_slots(object, "strand") %in% c("+", "-"))) {
                valid <- FALSE
                msg <- c(msg, "motif strand can only be '+', '-'")
              }

              if (length(motif_slots(object, "bkg")) > 0) {
                if (length(motif_slots(object, "bkg")) !=
                    nrow(motif_slots(object, "motif"))) {
                  valid <- FALSE
                  msg <- c(msg, "motif bkg frequencies and matrix nrow must match")
                }
              }

              if (length(motif_slots(object, "pseudoweight")) != 1) {
                valid <- FALSE
                msg <- c(msg, "pseudoweight must be a numeric of length 1")
              }

              if (length(motif_slots(object, "nsites")) > 1) {
                valid <- FALSE
                msg <- c(msg, "nsites must be a numeric of length 0 or 1")
              }

              if (length(motif_slots(object, "bkgsites")) > 1) {
                valid <- FALSE
                msg <- c(msg, "bkgsites must be a numeric of length 0 or 1")
              }

              if (length(motif_slots(object, "pval")) > 1) {
                valid <- FALSE
                msg <- c(msg, "pval must be a numeric of length 0 or 1")
              }

              if (length(motif_slots(object, "eval")) > 1) {
                valid <- FALSE
                msg <- c(msg, "eval must be a numeric of length 0 or 1")
              }

              if (object@type == "PCM") {
                test1 <- colSums(object@motif)
                if (length(unique(test1)) > 1) {
                  valid <- FALSE
                  msg <- c(msg, "motif of type PCM must have equal colSums")
                } else {
                  nsitestest <- unique(colSums(object@motif))
                  if (nsitestest != motif_slots(object, "nsites")) {
                    valid <- FALSE
                    msg <- c(msg, "PCM colSums and motif nsites must be equal")
                  }
                }
              }

              if (object@type == "PPM") {
                test2 <- colSums(object@motif)
                if (any(test2 > 1.01) || any(test2 < 0.99)) {
                  valid <- FALSE
                  msg <- c(msg, "motif of type PPM must have colSums of 1")
                }
                if (any(motif_slots(object, "motif") < 0)) {
                  valid <- FALSE
                  msg <- c(msg, "PPM values must be positive")
                }
              }

              if (valid) TRUE else msg

            })

# position frequency matrix (PFM): counts for hits at that position
# aka position count matrix (PCM)

  # WARNING: from what I can tell, PFMs are sometimes actually referring to
  # PPM

# position probability matrix (PPM): counts at that position divided by total
# aka position-specific probability matrix (PSPM)
  # the probability of a particular combo can be calculated by multiplying
    # pseudocounts: used when calculated PPM from a small dataset; otherwise
    # certain positions may have a probability of 0, regardless of other
    # positions

# position weight matrix (PWM): for each position of the PPM, apply the
# formulat: log2(position/bkgprob)
# aka position-specific weight matrix (PSWM), aka position-specific
# scoring matrix (PSSM)
    # without pseudocounts, there will be positions with -infinity
  # the probability of a particular combo can be calculated by addition;
  # if positive, more likely functional; if 0, random; if less than 0,
  # more likely a random site

# good ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
  # suggests a pseudocount of 0.8

# PSEUDOCOUNTS:
  
  # typically, 'pseudocounts' are added when converting from PFM to PPM.
  # the formula is thus modified:

    # position prob = (count + [pseudocount / 4]) / (num of seqs + pseudocount]

# pseudocounts can also be applied to PPMs:
  # prob * 100 = PFM
  # then recalculate PPM with pseudoweights formula

# workflow:
# multiple alignment --> PFM --> PPM --> PWM

# information content (IC) at a given position:
  # IC = log2(number of alph letters) - entropy(position)
  # for DNA/RNA: IC = 2 + entropy(position)
  # entropy(given position) = A[posprob*log2(posprob)] + C[..] + G[..] + T[..]

# PWMEnrich::motifIC

# seqLogo:::pwm2ic
  # ic[i] = 2 + sum( if x>0 (x * log2(x)) else 0 FOR ALL LETTERS)
  # e.g. for a = 0.25, c = 0.25, g = 0.25, t = 0.25:
  # 2 + 0.25*log2(0.25) + 0.25*log2(0.25) +..

# to get the bits for each letter of a position: prob * IC

# however the pwm2ic assumes a bkg of 0.25*4; this can be changed to get
# relative entropy:
  # 0.25*log2(0.25/lett prob) + etc..

# homer scores:
# Score for GGATGT

# score = log(pG1/0.25) + log(pG2/0.25) + log(pA3/0.25) + log(pT4/0.25) + 
#         log(pG5/0.25) + log(pT6/0.25)
