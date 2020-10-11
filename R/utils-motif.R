#' Motif-related utility functions.
#'
#' @param allow.nonfinite `logical(1)` If `FALSE`, then apply a pseudocount if
#'    non-finite values are found in the PWM. Note that if the motif has a
#'    pseudocount greater than zero and the motif is not currently of type PWM,
#'    then this parameter has no effect as the pseudocount will be
#'    applied automatically when the motif is converted to a PWM internally. This
#'    value is set to `FALSE` by default in order to stay consistent with
#'    pre-version 1.8.0 behaviour.
#' @param allow.zero `logical(1)` If `FALSE`, apply a pseudocount if zero values
#'    are found in the background frequencies.
#' @param alphabet `character(1)` One of `c('DNA', 'RNA')`.
#' @param bkg `numeric` Should be the same length as the alphabet length.
#' @param bkg1 `numeric` Vector of background probabilities for the first column.
#'    Only relevant if `method = "ALLR"`.
#' @param bkg2 `numeric` Vector of background probabilities for the second column.
#'    Only relevant if `method = "ALLR"`.
#' @param delete `logical(1)` Clear gap information from motif. If `FALSE`, then
#'    it can be reactivated  simply with `add_gap(motif)`.
#' @param gaploc `numeric` Motif gap locations. The gap occurs immediately after
#'    every position value. If missing, uses `round(ncol(motif) / 2)`.
#' @param letter `character(1)` Any DNA, RNA, or AA IUPAC letter. Ambiguity letters
#'    are accepted.
#' @param match `character` Sequence string to calculate score from.
#' @param method `character(1)` Column comparison metric. See [compare_motifs()]
#'    for details.
#' @param mingap `numeric` Minimum gap size. Must have one value for every location.
#'    If missing, set to 1.
#' @param maxgap `numeric` Maximum gap size. Must have one value for every location.
#'    If missing, set to 5.
#' @param motif Motif object to calculate scores from, or add/remove gap, or round.
#' @param motifs `list` A list of \linkS4class{universalmotif} motifs.
#' @param na.rm `logical` Remove columns where all values are `NA`.
#' @param nsites `numeric(1)` Number of sites motif originated from.
#' @param nsites1 `numeric(1)` Number of sites for the first column. Only relevant
#'    if `method = "ALLR"`.
#' @param nsites2 `numeric(1)` Number of sites for the second column. Only relevant
#'    if `method = "ALLR"`.
#' @param pct.tolerance `numeric(1)` or `character(1)` The minimum tolerated
#'    proportion each letter must represent per position in order not to be
#'    rounded off, either as a numeric value from 0 to 1 or a percentage written as
#'    a string from "0%" to "100%".
#' @param position `numeric` A numeric vector representing the frequency or
#'    probability for each alphabet letter at a specific position.
#' @param pseudocount `numeric(1)` Used to prevent zeroes in motif matrix.
#' @param pval `character(1)` String-formatted p-value.
#' @param relative_entropy `logical(1)` Calculate information content as
#'    relative entropy or Kullback-Leibler divergence.
#' @param schneider_correction `logical(1)` Apply sample size correction.
#' @param score `numeric(1)` Logodds motif score.
#' @param smooth `logical(1)` Apply pseudocount correction.
#' @param threshold `numeric(1)` Any number of numeric values between 0 and 1
#'    representing score percentage.
#' @param threshold.type `character` For `"total"`, a threshold of zero
#'    represents the minimum possible score. This means the range of scores that
#'    can be extracted is from the minimum to the maximum possible scores. For
#'    `"fromzero"`, a threshold of zero is a score of zero. This means the range
#'    of scores is from zero to the maximum. The `"total"` threshold type can
#'    only be used if no non-finite values are present in the PWM.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM' 'ICM')`.
#' @param use.freq `numeric(1)` Use regular motif or the respective `multifreq`
#'    representation.
#' @param x `numeric` First column for comparison.
#' @param y `numeric` Second column for comparison.
#'
#' @return
#'    For [consensus_to_ppm()] and [consensus_to_ppmAA()]: a numeric
#'    vector of length 4 and 20, respectively.
#'
#'    For [get_consensus()] and [get_consensusAA()]: a character vector
#'    of length 1.
#'
#'    For [get_matches()]: a `character` vector of motif matches.
#'
#'    For [motif_score()]: a named `numeric` vector of motif scores.
#'
#'    For [log_string_pval()]: a `numeric` vector of length 1.
#'
#'    For [position_icscore()]: a `numeric` vector of length 1.
#'
#'    For [ppm_to_icm()], [icm_to_ppm()], [pcm_to_ppm()],
#'    [ppm_to_pcm()], [ppm_to_pwm()], and [pwm_to_ppm()]: a `numeric`
#'    vector with length equal to input `numeric` vector.
#'
#'    For [prob_match()]: a `numeric` vector of probabilities.
#'
#'    For [round_motif()]: the input motif, rounded.
#'
#'    For [score_match()]: a `numeric` vector with the match motif score.
#'
#'    For [summarise_motifs()]: a `data.frame` with columns representing
#'    the [universalmotif-class] slots.
#'
#' @examples
#' data(examplemotif)
#' examplemotif0 <- examplemotif
#' examplemotif0["pseudocount"] <- 0
#' 
#' #######################################################################
#' ## add_gap
#' ## Add gap information to a motif.
#' m <- create_motif()
#' # Add a gap size 5-8 between positions 4 and 5:
#' m <- add_gap(m, gaploc = 4, mingap = 5, maxgap = 8)
#'
#' #######################################################################
#' ## compare_columns
#' ## Compare two numeric vectors using the metrics from compare_motifs()
#' compare_columns(c(0.5, 0.1, 0.1, 0.2), c(0.7, 0.1, 0.1, 0.1), "PCC")
#'
#' #######################################################################
#' ## consensus_to_ppm
#' ## Do the opposite of get_consensus. Note that loss of information is
#' ## inevitable. Generates a sequence matrix.
#' sapply(c("A", "G", "T", "B"), consensus_to_ppm)
#'
#' #######################################################################
#' ## consensus_to_ppmAA
#' ## Do the opposite of get_consensusAA and generate a motif matrix.
#' sapply(c("V", "A", "L"), consensus_to_ppmAA)
#'
#' #######################################################################
#' ## get_consensus
#' ## Get a consensus string from a DNA/RNA motif.
#' m <- create_motif()["motif"]
#' apply(m, 2, get_consensus)
#'
#' #######################################################################
#' ## get_consensusAA
#' ## Get a consensus string from an amino acid motif. Unless each position
#' ## is clearly dominated by a single amino acid, the resulting string will
#' ## likely be useless.
#' m <- create_motif(alphabet = "AA")["motif"]
#' apply(m, 2, get_consensusAA, type = "PPM")
#'
#' #######################################################################
#' ## get_match
#' ## Get all possible motif matches above input score
#' get_matches(examplemotif, 0)
#' get_matches(examplemotif0, 0, allow.nonfinite = TRUE)
#'
#' #######################################################################
#' ## get_scores
#' ## Get all possible scores for a motif
#' length(get_scores(examplemotif))
#' get_scores(examplemotif)
#' get_scores(examplemotif0, allow.nonfinite = TRUE)
#'
#' #######################################################################
#' ## icm_to_ppm
#' ## Do the opposite of ppm_to_icm.
#' m <- create_motif(type = "ICM")["motif"]
#' apply(m, 2, icm_to_ppm)
#'
#' #######################################################################
#' ## motif_score
#' ## Calculate motif score from different thresholds
#' m <- normalize(examplemotif)
#' motif_score(m, c(0, 0.8, 1))
#' motif_score(examplemotif0, c(0, 0.8, 1), allow.nonfinite = TRUE,
#'    threshold.type = "fromzero")
#'
#' #######################################################################
#' ## log_string_pval
#' ## Get the log of a string-formatted p-value
#' log_string_pval("1e-400")
#'
#' #######################################################################
#' ## pcm_to_ppm
#' ## Go from a count type motif to a probability type motif.
#' m <- create_motif(type = "PCM", nsites = 50)["motif"]
#' apply(m, 2, pcm_to_ppm, pseudocount = 1)
#'
#' #######################################################################
#' ## position_icscore
#' ## Similar to ppm_to_icm, except this calculates the position sum.
#' m <- create_motif()["motif"]
#' apply(m, 2, position_icscore, type = "PPM", bkg = rep(0.25, 4))
#'
#' #######################################################################
#' ## ppm_to_icm
#' ## Convert one column from a probability type motif to an information
#' ## content type motif.
#' m <- create_motif(nsites = 100, pseudocount = 0.8)["motif"]
#' apply(m, 2, ppm_to_icm, nsites = 100, bkg = rep(0.25, 4))
#'
#' #######################################################################
#' ## ppm_to_pcm
#' ## Do the opposite of pcm_to_ppm.
#' m <- create_motif()["motif"]
#' apply(m, 2, ppm_to_pcm, nsites = 50)
#'
#' #######################################################################
#' ## ppm_to_pwm
#' ## Go from a probability type motif to a weight type motif.
#' m <- create_motif()["motif"]
#' apply(m, 2, ppm_to_pwm, nsites = 100, bkg = rep(0.25, 4))
#'
#' #######################################################################
#' ## prob_match, prob_match_bkg
#' ## Calculate probability of a particular match based on background
#' ## frequencies
#' prob_match(examplemotif, "TATATAT")
#' ## Since this motif has a uniform background, the probability of
#' ## finding any motif hit within the sequence is equal
#' prob_match(examplemotif, "TATATAG")
#' m <- examplemotif
#' m["bkg"] <- c(0.3, 0.2, 0.2, 0.3)
#' prob_match(m, "TATATAT")
#' ## The prob_match_bkg alternative allows you to simply pass along the
#' ## background frequencies
#' prob_match_bkg(c(A=0.3, C=0.2, G=0.2, T=0.3), c("TATATAT", "TATATAG"))
#'
#' #######################################################################
#' ## pwm_to_ppm
#' ## Do the opposite of ppm_to_pwm.
#' m <- create_motif(type = "PWM")["motif"]
#' apply(m, 2, pwm_to_ppm, bkg = rep(0.25, 4))
#'
#' #######################################################################
#' ## Note that not all type conversions can be done directly; for those
#' ## type conversions which are unavailable, universalmotif just chains
#' ## together others (i.e. from PCM -> ICM => pcm_to_ppm -> ppm_to_icm)
#'
#' #######################################################################
#' ## round_motif
#' ## Round down letter scores to 0
#' m <- create_motif()
#' ## Remove letters from positions which are less than 5% of the total
#' ## position:
#' round_motif(m, pct.tolerance = 0.05)
#'
#' #######################################################################
#' ## score_match
#' ## Calculate score of a particular match
#' score_match(examplemotif, "TATATAT")
#' score_match(examplemotif, "TATATAG")
#' score_match(examplemotif0, "TATATAT", allow.nonfinite = TRUE)
#' score_match(examplemotif0, "TATATAG", allow.nonfinite = TRUE)
#'
#' #######################################################################
#' ## summarise_motifs
#' ## Create a data.frame of information based on a list of motifs.
#' m1 <- create_motif()
#' m2 <- create_motif()
#' m3 <- create_motif()
#' summarise_motifs(list(m1, m2, m3))
#'
#' #######################################################################
#' ## ungap
#' ## Unset motif's gap status. Does not delete actual gap data unless
#' ## delete = TRUE.
#' m <- create_motif()
#' m <- add_gap(m, 3, 2, 4)
#' m <- ungap(m)
#' # Restore gap data:
#' m <- add_gap(m)
#'
#' @seealso [create_motif()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name utils-motif
NULL

# TODO: rewrite examples

#' @rdname utils-motif
#' @export
add_gap <- function(motif, gaploc = ncol(motif) %/% 2, mingap = 1, maxgap = 5) {

  args <- as.list(environment())
  num_check <- check_fun_params(
    list(gaploc = args$gaploc, mingap = args$mingap, maxgap = args$maxgap),
    0, FALSE, TYPE_NUM
  )
  if (length(num_check)) stop(all_checks_collapse(num_check))

  maxlen <- max(c(length(gaploc), length(mingap), length(maxgap)))

  motif@gapinfo@gaploc <- rep(as.integer(gaploc), length.out = maxlen)
  motif@gapinfo@mingap <- rep(as.integer(mingap), length.out = maxlen)
  motif@gapinfo@maxgap <- rep(as.integer(maxgap), length.out = maxlen)

  if (any(c(gaploc, mingap, maxgap) < 0))
    stop("'gaploc', 'mingap', 'maxgap' must be positive numbers")

  if (any(mingap > maxgap))
    stop("'mingap' cannot be greater than the corresponding 'maxgap'")

  if (any(maxgap == 0))
    stop("'maxgap' must be greater than 0")

  motif@gapinfo@isgapped <- TRUE

  validObject_universalmotif(motif)

  motif

}

#' @rdname utils-motif
#' @export
compare_columns <- function(x, y, method,
                            bkg1 = rep(1 / length(x), length(x)),
                            bkg2 = rep(1 / length(y), length(y)),
                            nsites1 = 100, nsites2 = 100) {

  method <- match.arg(method, COMPARE_METRICS)
  if (length(x) != length(y))
    stop("length(x) does not match length(y) [", length(x), ", ", length(y), "]")
  compare_columns_cpp(x, y, bkg1, bkg2, nsites1, nsites2, method[1])

}

#' @rdname utils-motif
#' @export
consensus_to_ppm <- function(letter) {
  if (!letter %in% DNA_ALPHABET)
    stop(letter, " is not a DNA IUPAC symbol")
  consensus_to_ppmC(letter)
}

#' @rdname utils-motif
#' @export
consensus_to_ppmAA <- function(letter) {
  if (!letter %in% AA_ALPHABET)
    stop(letter, " is not an AA IUPAC symbol")
  consensus_to_ppmAAC(letter)
}

#' @rdname utils-motif
#' @export
get_consensus <- function(position, alphabet = "DNA", type = "PPM",
                          pseudocount = 1) {

  if (!type %in% c("PCM", "PPM", "PWM", "ICM"))
    stop("type must be one of ICM, PCM, PPM, PWM")
  if (!alphabet %in% c("DNA", "RNA"))
    stop("alphabet must be one of DNA, RNA")

  get_consensusC(position, alphabet, type, pseudocount)

}

#' @rdname utils-motif
#' @export
get_consensusAA <- function(position, type = "PPM", pseudocount = 0) {

  if (!type %in% c("PCM", "PPM", "PWM", "ICM"))
    stop("type must be one of ICM, PCM, PPM, PWM")

  get_consensusAAC(position, type, pseudocount)

}

#' @rdname utils-motif
#' @export
get_matches <- function(motif, score, allow.nonfinite = FALSE) {

  motif <- convert_motifs(motif)
  if (motif@type != "PWM")
    motif <- convert_type_internal(motif, "PWM")
  if (is.list(motif) && length(motif) == 1)
    motif <- motif[[i]]
  else if (is.list(motif))
    stop("a single motif must be input")
  if (any(is.infinite(motif@motif)) && !allow.nonfinite) {
    message(wmsg("Note: found -Inf values in PWM motif, normalizing. ",
      "Set `allow.nonfinite = TRUE` to prevent this behaviour."))
    motif <- normalize(motif)
  }

  if (any(is.infinite(motif@motif))) {
    score.range <- c(-Inf, motif_score(motif, threshold = 1,
        allow.nonfinite = TRUE, threshold.type = "fromzero"))
  } else {
    score.range <- motif_score(motif, allow.nonfinite = allow.nonfinite)
  }

  if (score > score.range[2])
    stop(wmsg("input score is greater than max possible score ",
              round(score.range[2], 3)))
  if (score < score.range[1])
    stop(wmsg("input score is less than min possible score ",
              round(score.range[1], 3)))

  score <- as.integer(score * 1000)

  alph <- rownames(motif@motif)

  score.mat <- matrix(suppressWarnings(as.integer(motif@motif * 1000)),
    nrow = nrow(motif@motif))
  if (anyNA(score.mat)) {
    score.mat[is.na(score.mat)] <- as.integer(min_max_ints()$min / ncol(score.mat))
  }

  alph.sort <- apply(score.mat, 2, order, decreasing = TRUE)
  for (i in seq_len(ncol(score.mat))) {
    score.mat[, i] <- score.mat[alph.sort[, i], i]
  }

  col.sort <- order(apply(score.mat, 2, max), decreasing = TRUE)
  score.mat <- score.mat[, col.sort]

  all.paths <- branch_and_bound_cpp_exposed(score.mat, score)

  all.paths <- all.paths[, order(col.sort), drop = FALSE]

  all.paths <- paths_alph_unsort(all.paths, alph.sort - 1)

  paths_to_alph(all.paths, alph)

}

#' @rdname utils-motif
#' @export
get_scores <- function(motif, allow.nonfinite = FALSE) {

  motif <- convert_motifs(motif)
  if (is.list(motif) && length(motif) > 1)
    stop("please only input a single motif")
  else if (is.list(motif)) motif <- motif[[1]]

  motif <- convert_type_internal(motif, "PWM")
  if (any(is.infinite(motif@motif)) && !allow.nonfinite) {
    message(wmsg("Note: found -Inf values in motif PWM, adding a pseudocount. ",
      "Set `allow.nonfinite = TRUE` to prevent this behaviour."))
    motif <- normalize(motif)
  }

  m <- motif@motif
  m[is.infinite(m)] <- NA
  m <- matrix(as.integer(m * 1000), nrow = nrow(m))

  ans <- expand_scores(m)
  ans[ans <= min_max_ints()$min] <- -Inf
  ans[is.finite(ans)] <- ans[is.finite(ans)] / 1000
  sort(ans, decreasing = TRUE)

}

#' @rdname utils-motif
#' @export
icm_to_ppm <- function(position) {
  icm_to_ppmC(position)
}

#' @rdname utils-motif
#' @export
motif_score <- function(motif, threshold = c(0, 1), use.freq = 1,
  allow.nonfinite = FALSE, threshold.type = c("total", "fromzero")) {

  threshold.type <- match.arg(threshold.type, c("total", "fromzero"))

  if (any(threshold > 1) || any(threshold < 0))
    stop("For 'threshold', please only use values between 0 and 1")

  motif <- convert_motifs(motif)

  if (is.list(motif) && length(motif) > 1)
    stop("Please only input a single motif")
  else if (is.list(motif))
    motif <- motif[[1]]

  if (!is(motif, "universalmotif"))
    stop("Unknown motif object")

  if (use.freq == 1) {

    motif <- convert_type_internal(motif, "PWM")

    if (any(is.infinite(motif@motif)) && !allow.nonfinite) {
      message(wmsg("Note: found -Inf values in motif PWM, adding a pseudocount. ",
        "Set `allow.nonfinite = TRUE` to prevent this behaviour."))
      motif <- normalize(motif)
    } else if (any(is.infinite(motif@motif)) && threshold.type == "total") {
      stop(wmsg("Score calculates cannot be performed for `threshold.type = \"total\"`",
        "if non-finite values are present in the PWM and `allow.nonfinite = TRUE`"))
    }

    mat <- matrix(suppressWarnings(as.integer(motif@motif * 1000)), nrow = nrow(motif@motif))

  } else {

    if (!as.character(use.freq[1]) %in% names(motif@multifreq))
      stop("missing appropriate multifreq slot [", use.freq, "]")

    mat <- motif@multifreq[[as.character(use.freq[1])]]

    mat <- MATRIX_ppm_to_pwm(mat, bkg = motif@bkg[rownames(mat)],
                             pseudocount = motif@pseudocount,
                             nsites = motif@nsites)

    mat <- matrix(suppressMessages(as.integer(mat * 1000)), nrow = nrow(mat))

  }

  if (threshold.type == "total") {
    s.max <- sum(apply(mat, 2, max))
    s.min <- sum(apply(mat, 2, min))
    s.total <- abs(s.max) + abs(s.min)
    out <- s.total * threshold - abs(s.min)
  } else {
    s.max <- sum(apply(mat, 2, max, na.rm = TRUE))
    out <- s.max * threshold
  }

  s.max <- sum(apply(mat, 2, max, na.rm = TRUE))

  names(out) <- paste0(as.character(threshold * 100), "%")

  out / 1000

}

#' @rdname utils-motif
#' @export
log_string_pval <- function(pval) {
  pval_str2double(as.character(pval))
}

#' @rdname utils-motif
#' @export
pcm_to_ppm <- function(position, pseudocount = 0) {
  pcm_to_ppmC(position, pseudocount)
}

#' @rdname utils-motif
#' @export
position_icscore <- function(position, bkg = numeric(), type = "PPM",
                             pseudocount = 1, nsites = 100,
                             relative_entropy = FALSE,
                             schneider_correction = FALSE) {

  if (is.null(bkg) || missing(bkg) || length(bkg) == 0) {
    bkg <- rep(1 / length(position), length(position))
  }
  if (!type %in% c("ICM", "PPM", "PWM", "PCM"))
    stop("type must be one of ICM, PCM, PPM, PWM")

  if (relative_entropy && schneider_correction)
    stop("relative_entropy and schneider_correction cannot both be TRUE")

  if (!schneider_correction)
    position_icscoreC(position, bkg, type, pseudocount, nsites, relative_entropy)
  else
    sum(ppm_to_icm(position, bkg, TRUE, nsites))

}

#' @rdname utils-motif
#' @export
ppm_to_icm <- function(position, bkg = numeric(), schneider_correction = FALSE,
                       nsites = 100, relative_entropy = FALSE) {
  # NOTE: Basic IC computation assumes uniform bkg frequencies!
  #       For different bkg frequencies: Relative entropy or Kullback-Leibler
  #       (KL) divergence
  if (is.null(bkg) || missing(bkg) || length(bkg) == 0) {
    bkg <- rep(1 / length(position), length(position))
  }

  if (relative_entropy && schneider_correction)
    stop("relative_entropy and schneider_correction cannot both be TRUE")

  if (relative_entropy) {

    ppm_to_icmC(position, bkg, TRUE)

  } else {

    if (schneider_correction) {
      if (length(position) != 4)
        stop("schneider correction is only available for motifs where nrow(motif) == 4")
      p <- as.integer(ppm_to_pcm(position, nsites = nsites))
      names(bkg) <- DNA_BASES
      if (requireNamespace("TFBSTools", quietly = TRUE)) {
        total_ic <- TFBSTools::toICM(matrix(p, nrow = 4, dimnames = list(DNA_BASES)),
                                     pseudocounts = 0, schneider = TRUE,
                                     bg = bkg)
        total_ic <- sum(as.vector(total_ic))
      } else {
        stop("The 'TFBSTools' package is required for 'schneider_correction'")
      }

    } else {

      height_after <- -sum(position * log2(position), na.rm = TRUE)
      total_ic <- log2(length(position)) - height_after

    }

    ic <- position * total_ic
    ic

  }

}

#' @rdname utils-motif
#' @export
ppm_to_pcm <- function(position, nsites = 100) {
  ppm_to_pcmC(position, nsites)
}

#' @rdname utils-motif
#' @export
ppm_to_pwm <- function(position, bkg = numeric(), pseudocount = 1, nsites = 100,
                       smooth = TRUE) {
  if (missing(bkg) || length(bkg) == 0)
    bkg <- rep(1 / length(position), length(position))
  if (length(nsites) == 0) nsites <- 100
  if (smooth && pseudocount != 0) {
    position <- ppm_to_pcm(position, nsites = nsites)
    position <- pcm_to_ppm(position, pseudocount = pseudocount)
  }
  for (i in seq_along(position)) {
    position[i] <- log2(position[i] / bkg[i])
  }
  position
}

#' @rdname utils-motif
#' @export
prob_match <- function(motif, match, allow.zero = TRUE) {

  if (missing(motif) || missing(match))
    stop("motif and/or match are missing")

  if (!is.character(match) || length(match) < 1)
    stop("match must be a non-empty character vector")

  motif <- convert_motifs(motif)

  if (is.list(motif) && length(motif) != 1)
    stop("Please only input a single motif")
  else if (is.list(motif))
    motif <- motif[[1]]

  if (!is(motif, "universalmotif"))
    stop("Unknown motif class")

  if (motif@type != "PPM")
    motif <- convert_type_internal(motif, "PPM")

  bkg <- motif@bkg[seq_len(nrow(motif))]

  if (!allow.zero && any(bkg == 0)) {
    message(wmsg("Note: found zero values in motif background frequencies, ",
        "applying pseudocount. Set `allow.nonfinite = TRUE` to prevent this behaviour."))
    pseudo <- motif@pseudocount
    if (pseudo == 0) {
      message("Note: motif has a pseudocount of zero, using 1")
      pseudo <- 1
    }
    nsites <- motif@nsites
    if (length(nsites) == 0) {
      message("Note: motif has no nsites info, using 100")
      nsites <- 100
    }
    bkg <- bkg + (pseudo / nsites) / length(bkg)
  }

  match <- lapply(match, safeExplode)
  if (any(vapply(match, length, integer(1)) != ncol(motif)))
    stop(wmsg("motif length [", ncol(motif), "] and match length [",
              paste0(vapply(match, length, integer(1)), collapse = ","),
              "] are not equal"))

  if (!all(unique(unlist(match)) %in% rownames(motif@motif)))
    stop("Found letters in match not found in motif")

  prob <- rep(1, length(match))

  # mat <- matrix(motif@motif,
  #               nrow = nrow(motif@motif),
  #               dimnames = dimnames(motif@motif))

  # for (i in seq_len(ncol(motif))) {
  #   for (j in seq_along(prob)) {
  #     # prob[j] <- prob[j] * mat[match[[j]][i], i]
  #     prob[j] <- prob[j] * bkg[match[[j]][i]]
  #   }
  # }

  for (i in seq_along(prob)) {
    for (j in seq_len(ncol(motif))) {
      prob[i] <- prob[i] * bkg[match[[i]][j]]
    }
  }

  prob

}

#' @rdname utils-motif
#' @export
prob_match_bkg <- function(bkg, match) {

  if (!is.character(match) || length(match) < 1)
    stop("match must be a non-empty character vector")

  if (!is.numeric(bkg) || length(bkg) < 1)
    stop("bkg must be a non-empty numeric vector")

  if (is.null(names(bkg)))
    stop("bkg must be a named vector")

  match <- lapply(match, safeExplode)
  lets <- sort_unique_cpp(unname(unlist(match)))

  if (any(!lets %in% names(bkg)))
    stop("found letters in match not in bkg vector names")

  prob <- rep(1, length(match))

  for (i in seq_along(prob)) {
    for (j in seq_along(match[[i]])) {
      prob[i] <- prob[i] * bkg[match[[i]][j]]
    }
  }

  prob

}

#' @rdname utils-motif
#' @export
pwm_to_ppm <- function(position, bkg = numeric()) {
  if (missing(bkg) || length(bkg) == 0)
    bkg <- rep(1 / length(position), length(position))
  position <- vapply(position, function(x) 2 ^ x, numeric(1))
  if (sum(position) > 0.99 && sum(position) < 1.01) return(position)
  for (i in seq_along(position)) position[i] <- position[i] * bkg[i]
  if (sum(position) > 0.99 && sum(position) < 1.01) return(position)
  message("Note: position does not add up to 1; normalizing..")
  pos_missing <- sum(position)
  position <- position / pos_missing
  position
}

#' @rdname utils-motif
#' @export
round_motif <- function(motif, pct.tolerance = 0.05) {
  if (!is(motif, "universalmotif"))
    stop("'motif' must be a 'universalmotif' object")
  if (is.character(pct.tolerance))
    pct.tolerance <- as.numeric(gsub("%", "", pct.tolerance, fixed = TRUE)) / 100
  validObject_universalmotif(motif)
  type <- motif@type
  motif <- convert_type_single(motif, "PPM", 0)
  motif@motif <- round_motif_cpp(motif@motif, pct.tolerance)
  convert_type(motif, type)
}

#' @rdname utils-motif
#' @export
score_match <- function(motif, match, allow.nonfinite = FALSE) {

  if (missing(motif) || missing(match))
    stop("motif and/or match are missing")

  if (!is.character(match) || length(match) < 1)
    stop("match must be a non-empty character vector")

  motif <- convert_motifs(motif)

  if (is.list(motif) && length(motif) != 1)
    stop("Please only input a single motif")
  else if (is.list(motif))
    motif <- motif[[1]]

  if (!is(motif, "universalmotif"))
    stop("Unknown motif class")

  if (motif@type != "PWM")
    motif <- convert_type_internal(motif, "PWM")

  if (!allow.nonfinite && any(is.infinite(motif@motif))) {
    message(wmsg("Note: found -Inf values in motif PWM, applying pseudocount. ",
      "Set `allow.nonfinite = TRUE` to prevent this behaviour."))
    motif <- normalize(motif)
  }

  match <- lapply(match, safeExplode)
  if (any(vapply(match, length, integer(1)) != ncol(motif)))
    stop(wmsg("motif length [", ncol(motif), "] and match length [",
              paste0(vapply(match, length, integer(1)), collapse = ","),
              "] are not equal"))

  if (!all(unique(unlist(match)) %in% rownames(motif@motif)))
    stop("Found letters in match not found in motif")

  score <- rep(0, length(match))

  mat <- matrix(motif@motif,
                nrow = nrow(motif@motif),
                dimnames = dimnames(motif@motif))

  for (i in seq_len(ncol(mat))) {
    for (j in seq_along(score)) {
    score[j] <- score[j] + mat[match[[j]][i], i]
    }
  }

  score[is.finite(score)] <- as.integer(score[is.finite(score)] * 1000) / 1000
  score

}

#' @rdname utils-motif
#' @export
summarise_motifs <- function(motifs, na.rm = TRUE) {

  # ~0.05 seconds for entire MotifDb library
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  classcheck <- vapply(motifs, function(x) !is(x, "universalmotif"), logical(1))
  if (any(classcheck)) stop("all motifs must be 'universalmotif'")

  # Very strange bug where it fails if there's only a single motif
  len1 <- if (length(motifs) == 1) TRUE else FALSE
  if (len1) motifs <- c(motifs, motifs)

  out <- summarise_motifs_cpp(motifs)
  out <- out[, c("name", "altname", "family", "organism", "consensus", "alphabet",
                 "strand", "icscore", "nsites", "bkgsites", "pval", "qval", "eval")]
  if (na.rm) out <- Filter(function(x) !all(is.na(x)), out)

  if (len1) out <- out[1, ]

  out

}

#' @rdname utils-motif
#' @export
ungap <- function(motif, delete = FALSE) {

  if (delete) {
    motif@gapinfo <- new("universalmotif_gapped")
  } else {
    motif@gapinfo@isgapped <- FALSE
  }

  motif

}
