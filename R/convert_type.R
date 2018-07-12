#' Convert \linkS4class{universalmotif} type.
#'
#' Switch between position count matrix (PCM), position probability matrix
#' (PPM), position weight matrix (PWM), and information count matrix (ICM)
#' types.
#'
#' @param motifs \linkS4class{universalmotif} objects.
#' @param type Character. 'PCM', 'PPM', 'PWM', or 'ICM'.
#' @param pseudocount Numeric. Correction to be applied to prevent \code{-Inf}
#'                    from apearing in PWM matrices. If missing,
#'                    the pseudocount stored in the
#'                    \linkS4class{universalmotif} 'pseudocount' slot will be
#'                    used, which defaults to 0.8 (the suggested value from
#'                    \insertCite{pseudo;textual}{universalmotif}).
#' @param nsize_correction Logical. If true, the ICM
#'                          at each position will be corrected to account
#'                          for small sample sizes. Only used if 
#'                          \code{relative_entropy = FALSE}.
#' @param relative_entropy Logical. If true, the ICM will be
#'                         calculated as relative entropy. (See details.)
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return \linkS4class{universalmotif} objects.
#'
#' @details
#'    Position count matrix (PCM), also known as position frequency matrix
#'    (PFM). For n sequences from which the motif was built, each position is
#'    represented by the numbers of each letter at that position. In theory
#'    all positions should have sums equal to n, but not all databases are
#'    this consistent. If converting from another type to PCM, column sums
#'    will be equal to the 'nsites' slot; if empty, 100 is used.
#'
#'    Position probability matrix (PPM), also known as position frequency
#'    matrix (PFM). At each position, the probability of individual letters
#'    is calculated by dividing the count for that letter by the total sum of
#'    counts at that position (\code{letter_count / position_total}). 
#'    As a result, each position will sum to 1. Letters with counts of 0 will
#'    thus have a probability of 0, which can be undesirable when searching for
#'    motifs in a set of sequences. To avoid this a pseudocount can be added 
#'    (\code{(letter_count + pseudocount) / (position_total + pseudocount)}).
#'
#'    Position weight matrix (PWM; \insertCite{pwm;textual}{universalmotif}),
#'    also known as position-specific weight
#'    matrix (PSWM), position-specific scoring matrix (PSSM), or
#'    log-odds matrix. At each position, each letter is represented by it's
#'    log-likelihood (\code{log2(letter_probability / background_probility)}),
#'    which is normalized using the background letter frequencies. A PWM matrix
#'    is constructed from a PPM; if any position has 0-probability letters to
#'    which pseudocounts were not added, then the final log-likelihodd of these
#'    letters will be \code{-Inf}.
#'
#'    Information content matrix (ICM; \insertCite{icm;textual}{universalmotif}).
#'    An ICM is a PPM where each letter probability is multiplied by the total
#'    information content at that position. The information content of each
#'    position is determined as: \code{totalIC - Hi}, where the total information
#'    totalIC
#'    
#'    \code{totalIC <- log2(alphabet_length)}, and the Shannon entropy 
#'    \insertCite{shannon}{universalmotif} for a specific 
#'    position (Hi)
#'
#'    \code{Hi <- -sum(sapply(alphabet_frequencies,
#'                            function(x) x * log(2))}.
#'
#'    As a result, the total sum or height of each position is representative of
#'    it's sequence conservation, measured in the unit 'bits', which is a unit
#'    of energy (\insertCite{bits;textual}{universalmotif}; see
#'    \url{https://fr-s-schneider.ncifcrf.gov/logorecommendations.html}
#'    for more information). However not all programs will calculate
#'    information content the same; some will 'correct' the total information
#'    content at each position using a correction factor as described by
#'    \insertCite{correction;textual}{universalmotif}. This correction can
#'    applied by setting \code{nsize_correction = TRUE}, however it will only
#'    be applied if the 'nsites' slot is not empty. This is done using
#'    \code{TFBSTools:::schneider_correction}
#'    \insertCite{tfbstools}{universalmotif}. As such, converting from an ICM to
#'    which some form of correction has been applied will result in a
#'    PCM/PPM/PWM with slight inaccuracies.
#'
#'    Another method of calculating information content is calculating the
#'    relative entropy, also known as Kullback-Leibler divergence 
#'    \insertCite{kl}{universalmotif}. This accounts for background
#'    frequencies, which
#'    can be useful for genomes with a heavy imbalance in letter frequencies.
#'    For each position, the individual letter frequencies are calculated as
#'    \code{letter_freq * log2(letter_freq / bkg_freq)}. When calculating 
#'    information content using Shannon entropy, the maximum content for
#'    each position will always be \code{log2(alphabet_length)}; this does
#'    not hold for information content calculated as relative entropy.
#'    Please note that conversion from ICM assumes the information content
#'    was _not_ calculated as relative entropy.
#'
#' @examples
#' jaspar.pcm <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                       package = "universalmotif"))
#'
#' # pseudocounts default to 0.8
#' jaspar.pwm <- convert_type(jaspar.pcm, type = "PPM")
#' 
#' # setting pseudocounts to 0 will prevent any correction from being
#' # applied to PPM/PWM matrices
#' jaspar.pwm <- convert_type(jaspar.pcm, type = "PWM", pseudocount = 0)
#'
#' @references
#'    \insertRef{pseudo}{universalmotif}
#'
#'    \insertRef{pwm}{universalmotif}
#'
#'    \insertRef{icm}{universalmotif}
#'
#'    \insertRef{bits}{universalmotif}
#'
#'    \insertRef{correction}{universalmotif}
#'
#'    \insertRef{tfbstools}{universalmotif}
#'
#'    \insertRef{kl}{universalmotif}
#'
#'    \insertRef{shannon}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
convert_type <- function(motifs, type, pseudocount,
                         nsize_correction = FALSE, 
                         relative_entropy = FALSE,
                         BPPARAM = SerialParam()) {

  if (!missing(pseudocount)) {
    if (!is.numeric(pseudocount)) stop("pseudocount must be a numeric vector")
    if (length(pseudocount) > 1) stop("pseudocount must a length one vector")
  }

  if (!any(is.logical(nsize_correction), is.logical(relative_entropy))) {
    stop("'nsize_correction' and 'relative_entropy' must be TRUE or FALSE")
  }

  if (!type %in% c("PCM", "PPM", "PWM", "ICM")) {
    stop("unrecognized 'type'")
  }

  motif <- motifs
  if (class(motif) == "list") {
    margs <- list(type = type) 
    if (!missing(pseudocount)) margs <- c(margs, list(pseudocount = pseudocount))
    motif <- bplapply(motif, function(x) do.call(convert_type,
                                                 c(list(motifs = x), 
                                                   margs)),
                      BPPARAM = BPPARAM)
    return(motif)
  }

  nalph <- nrow(motif["motif"])
  nbkg <- length(motif["bkg"])
  if (nalph != nbkg) stop("length of bkg must match alphabet length")
    
  CLASS_IN <- .internal_convert(motif, BPPARAM = BPPARAM)
  motif <- convert_motifs(motif, BPPARAM = BPPARAM)

  in_type <- motif["type"]

  if (in_type == type) return(motif)

  if (missing(pseudocount)) pseudocount <- motif["pseudocount"]
  bkg <- motif["bkg"]

  if (length(motif["nsites"]) == 0) nsites <- 100 else nsites <- motif["nsites"]

  alph <- rownames(motif@motif)

  # PCM in:
  if (in_type == "PCM") {
    if (type == "PPM") {
      motif@motif <- apply(motif["motif"], 2, pcm_to_ppmC,
                           pseudocount = pseudocount)
      motif["type"] <- "PPM"
    } else if (type == "PWM") {
      motif@motif <- apply(motif["motif"], 2, pcm_to_ppmC,
                           pseudocount = pseudocount)
      motif@motif <- apply(motif["motif"], 2, ppm_to_pwmC,
                           bkg = bkg,
                           pseudocount = pseudocount,
                           nsites = nsites)
      motif["type"] <- "PWM"
    } else if (type == "ICM") {
      motif@motif <- apply(motif["motif"], 2, pcm_to_ppmC,
                           pseudocount = pseudocount)
      if (nsize_correction) {
        motif@motif <- apply(motif["motif"], 2, ppm_to_icm,
                             bkg = bkg, 
                             nsites = nsites,
                             schneider_correction = nsize_correction,
                             relative_entropy = relative_entropy)
      } else {
        motif@motif <- apply(motif["motif"], 2, ppm_to_icmC,
                             bkg = bkg, 
                             relative_entropy = relative_entropy)
      }
      motif["type"] <- "ICM"
    }
  }

  # PPM in:
  if (in_type == "PPM") {
    if (type == "PCM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pcmC,
                           nsites = nsites)
      motif["type"] <- "PCM"
    } else if (type == "PWM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pwmC,
                           bkg = bkg,
                           pseudocount = pseudocount,
                           nsites = nsites)
      motif["type"] <- "PWM"
    } else if (type == "ICM") {
      if (nsize_correction) {
        motif@motif <- apply(motif["motif"], 2, ppm_to_icm,
                             bkg = bkg, 
                             nsites = nsites,
                             schneider_correction = nsize_correction,
                             relative_entropy = relative_entropy)
      } else {
        motif@motif <- apply(motif["motif"], 2, ppm_to_icmC,
                             bkg = bkg, 
                             relative_entropy = relative_entropy)
      }
      motif["type"] <- "ICM"
    }
  }

  # PWM in:
  if (in_type == "PWM") {
    if (type == "PCM") {
      motif@motif <- apply(motif["motif"], 2, pwm_to_ppmC,
                           bkg = bkg)
      motif@motif <- apply(motif["motif"], 2, ppm_to_pcmC,
                           nsites = nsites)
      motif["type"] <- "PCM"
    } else if (type == "PPM") {
      motif@motif <- apply(motif["motif"], 2, pwm_to_ppmC,
                           bkg = bkg)
      motif["type"] <- "PPM"
    } else if (type == "ICM") {
      motif@motif <- apply(motif["motif"], 2, pwm_to_ppmC,
                           bkg = bkg)
      if (nsize_correction) {
        motif@motif <- apply(motif["motif"], 2, ppm_to_icm,
                             bkg = bkg, 
                             nsites = nsites,
                             schneider_correction = nsize_correction,
                             relative_entropy = relative_entropy)
      } else {
        motif@motif <- apply(motif["motif"], 2, ppm_to_icmC,
                             bkg = bkg, 
                             relative_entropy = relative_entropy)
      }
      motif["type"] <- "ICM"
    }
  }

  # ICM in:
  if (in_type == "ICM") {
    motif@motif <- apply(motif["motif"], 2, icm_to_ppmC)
    if (type == "PCM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pcmC,
                           nsites = nsites)
      motif["type"] <- "PCM"
    } else if (type == "PPM") {
      motif["type"] <- "PPM"
    } else if (type == "PWM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pwmC,
                           bkg = bkg,
                           pseudocount = pseudocount,
                           nsites = nsites)
      motif["type"] <- "PWM"
    }
  }

  rownames(motif@motif) <- alph

  motif <- .internal_convert(motif, CLASS_IN, BPPARAM = BPPARAM)
  motif

}
