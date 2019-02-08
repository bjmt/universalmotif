#' Convert \linkS4class{universalmotif} type.
#'
#' Switch between position count matrix (PCM), position probability matrix
#' (PPM), position weight matrix (PWM), and information count matrix (ICM)
#' types.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#' @param pseudocount `numeric(1)` Correction to be applied to prevent `-Inf`
#'   from apearing in PWM matrices. If missing, the pseudocount stored in the
#'   \linkS4class{universalmotif} 'pseudocount' slot will be
#'   used, which defaults to 0.8; the suggested value from
#'   \insertCite{pseudo;textual}{universalmotif}.
#' @param nsize_correction `logical(1)` If true, the ICM
#'   at each position will be corrected to account
#'   for small sample sizes. Only used if `relative_entropy = FALSE`.
#' @param relative_entropy `logical(1)` If true, the ICM will be
#'   calculated as relative entropy. See details.
#'
#' @return See [convert_motifs()] for possible output motif objects.
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
#'    counts at that position (`letter_count / position_total`). 
#'    As a result, each position will sum to 1. Letters with counts of 0 will
#'    thus have a probability of 0, which can be undesirable when searching for
#'    motifs in a set of sequences. To avoid this a pseudocount can be added 
#'    (`(letter_count + pseudocount) / (position_total + pseudocount)`).
#'
#'    Position weight matrix (PWM; \insertCite{pwm;textual}{universalmotif}),
#'    also known as position-specific weight
#'    matrix (PSWM), position-specific scoring matrix (PSSM), or
#'    log-odds matrix. At each position, each letter is represented by it's
#'    log-likelihood (`log2(letter_probability / background_probility)`),
#'    which is normalized using the background letter frequencies. A PWM matrix
#'    is constructed from a PPM; if any position has 0-probability letters to
#'    which pseudocounts were not added, then the final log-likelihodd of these
#'    letters will be `-Inf`.
#'
#'    Information content matrix (ICM; \insertCite{icm;textual}{universalmotif}).
#'    An ICM is a PPM where each letter probability is multiplied by the total
#'    information content at that position. The information content of each
#'    position is determined as: `totalIC - Hi`, where the total information
#'    totalIC
#'    
#'    `totalIC <- log2(alphabet_length)`, and the Shannon entropy 
#'    \insertCite{shannon}{universalmotif} for a specific 
#'    position (Hi)
#'
#'    `Hi <- -sum(sapply(alphabet_frequencies, function(x) x * log(2))`.
#'
#'    As a result, the total sum or height of each position is representative of
#'    it's sequence conservation, measured in the unit 'bits', which is a unit
#'    of energy (\insertCite{bits;textual}{universalmotif}; see
#'    \url{https://fr-s-schneider.ncifcrf.gov/logorecommendations.html}
#'    for more information). However not all programs will calculate
#'    information content the same; some will 'correct' the total information
#'    content at each position using a correction factor as described by
#'    \insertCite{correction;textual}{universalmotif}. This correction can
#'    applied by setting `nsize_correction = TRUE`, however it will only
#'    be applied if the 'nsites' slot is not empty. This is done using
#'    `TFBSTools:::schneider_correction`
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
#'    `letter_freq * log2(letter_freq / bkg_freq)`. When calculating 
#'    information content using Shannon entropy, the maximum content for
#'    each position will always be `log2(alphabet_length)`; this does
#'    not hold for information content calculated as relative entropy.
#'    Please note that conversion from ICM assumes the information content
#'    was _not_ calculated as relative entropy.
#'
#' @examples
#' jaspar.pcm <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                       package = "universalmotif"))
#'
#' ## pseudocounts default to 0.8
#' jaspar.pwm <- convert_type(jaspar.pcm, type = "PPM")
#' 
#' ## setting pseudocounts to 0 will prevent any correction from being
#' ## applied to PPM/PWM matrices
#' jaspar.pwm <- convert_type(jaspar.pcm, type = "PWM", pseudocount = 0)
#'
#' @references
#'    \insertRef{kl}{universalmotif}
#'
#'    \insertRef{pseudo}{universalmotif}
#'
#'    \insertRef{correction}{universalmotif}
#'
#'    \insertRef{icm}{universalmotif}
#'
#'    \insertRef{bits}{universalmotif}
#'
#'    \insertRef{shannon}{universalmotif}
#'
#'    \insertRef{pwm}{universalmotif}
#'
#'    \insertRef{tfbstools}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [convert_motifs()]
#' @export
convert_type <- function(motifs, type, pseudocount, nsize_correction = FALSE, 
                         relative_entropy = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!type %in% c("PCM", "PPM", "PWM", "ICM")) {
    type_check <- paste0(" * Incorrect 'type': expected `PCM`, `PPM`, `PWM` or `ICM`; got `",
                         type, "`")
    all_checks <- c(all_checks, type_check)
  }
  char_check <- check_fun_params(list(type = args$type), 1, FALSE, "character")
  num_check <- check_fun_params(list(pseudocount = args$pseudocount),
                                1, TRUE, "numeric")
  logi_check <- check_fun_params(list(nsize_correction = args$nsize_correction,
                                      relative_entropy = args$relative_entropy),
                                 c(1, 1), c(FALSE, FALSE), "logical")
  all_checks <- c(all_checks, char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

   if (!type %in% c("PCM", "PPM", "PWM", "ICM")) {
    stop("unrecognized 'type'")
  }

  if (is.list(motifs)) was_list <- TRUE else was_list <- FALSE

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, "character")
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  if (missing(pseudocount)) pseudocount <- NULL

  motifs <- lapply(motifs, function(x) convert_type_single(x, type, pseudocount,
                                                             nsize_correction,
                                                             relative_entropy)) 

  motifs <- .internal_convert(motifs, unique(CLASS_IN))
  if (length(motifs) == 1 && !was_list) motifs <- motifs[[1]]
  motifs

}

convert_type_single <- function(motif, type, pseudocount, nsize_correction = FALSE, 
                         relative_entropy = FALSE) {

  if (nrow(motif["motif"]) != length(motif["bkg"]))
    stop("convert_type_sing: length of bkg must match alphabet length [INTERNAL]")
    
  in_type <- motif["type"]

  if (in_type == type) return(motif)

  if (is.null(pseudocount)) pseudocount <- motif["pseudocount"]
  bkg <- motif["bkg"]

  if (length(motif["nsites"]) == 0) nsites <- 100 else nsites <- motif["nsites"]

  motif <- switch(in_type,
                  "PCM" = convert_from_pcm(motif, type, pseudocount, bkg, nsites,
                                           nsize_correction, relative_entropy),
                  "PPM" = convert_from_ppm(motif, type, pseudocount, bkg, nsites,
                                           nsize_correction, relative_entropy),
                  "PWM" = convert_from_pwm(motif, type, pseudocount, bkg, nsites,
                                           nsize_correction, relative_entropy),
                  "ICM" = convert_from_icm(motif, type, pseudocount, bkg, nsites))

  motif

}

convert_from_pcm <- function(motif, type, pseudocount, bkg, nsites, 
                             nsize_correction = FALSE, relative_entropy = FALSE) {

  alph <- rownames(motif@motif)

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

  rownames(motif@motif) <- alph

  motif

}

convert_from_ppm <- function(motif, type, pseudocount, bkg, nsites,
                             nsize_correction = FALSE, relative_entropy = FALSE) {

  alph <- rownames(motif@motif)

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

  rownames(motif@motif) <- alph

  motif

}

convert_from_pwm <- function(motif, type, pseudocount, bkg, nsites,
                             nsize_correction = FALSE, relative_entropy = FALSE) {

  alph <- rownames(motif@motif)

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

  rownames(motif@motif) <- alph

  motif

}

convert_from_icm <- function(motif, type, pseudocount, bkg, nsites) {

  alph <- rownames(motif@motif)

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
  
  rownames(motif@motif) <- alph

  motif

}
