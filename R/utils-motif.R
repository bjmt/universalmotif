#' Motif-related utility functions.
#'
#' @param alphabet `character(1)` One of `c('DNA', 'RNA')`.
#' @param bkg `Numeric` Should be the same length as the alphabet length.
#' @param BP `logical(1)` Use the \pkg{BiocParallel} package. See
#'    [BiocParallel::register()] to change the default backend.
#' @param db.motifs `list` Database motifs.
#' @param letter `character(1)` Any DNA, RNA, or AA IUPAC letter. Ambiguity letters
#'    are accepted.
#' @param match `character(1)` Sequence string to calculate score from.
#' @param method `character(1)` One of `c('PCC', 'MPCC', 'EUCL', 'MEUCL',
#'    'SW', 'MSW', 'KL', 'MKL')`. See [compare_motifs()].
#' @param min.mean.ic `numeric(1)` See [compare_motifs()].
#' @param min.overlap `numeric(1)` Minimum required motif overlap. See
#'    [compare_motifs()].
#' @param motif Motif object to calculate scores from.
#' @param motifs `list` A list of \linkS4class{universalmotif} motifs.
#' @param na.rm `logical` Remove columns where all values are `NA`.
#' @param normalise.scores `logical(1)` See [compare_motifs()].
#' @param nsites `numeric(1)` Number of sites motif originated from.
#' @param position `numeric` A numeric vector representing the frequency or
#'    probability for each alphabet letter at a specific position.
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param pseudocount `numeric(1)` Used to prevent zeroes in motif matrix.
#' @param rand.tries `numeric(1)` Number of random motifs to create for
#'    P-value computation.
#' @param relative_entropy `logical(1)` Calculate information content as
#'    relative entropy or Kullback-Leibler divergence.
#' @param schneider_correction `logical(1)` Apply sample size correction.
#' @param score `numeric(1)` Logodds motif score.
#' @param shuffle.db `logical(1)` Shuffle `db.motifs` rather than
#'    generate random motifs with [create_motif()].
#' @param shuffle.k `numeric(1)` See [shuffle_motifs()].
#' @param shuffle.leftovers `character(1)` See [shuffle_motifs()].
#' @param shuffle.method `character(1)` See [shuffle_motifs()].
#' @param smooth `logical(1)` Apply pseudocount correction.
#' @param threshold `numeric(1)` Any number of numeric values between 0 and 1
#'    representing score percentage.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM' 'ICM')`.
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
#'    For [make_DBscores()]: a `data.frame` with score distributions for the
#'    input database.
#'
#'    For [motif_score()]: a named `numeric` vector of motif scores.
#'
#'    For [position_icscore()]: a `numeric` vector of length 1.
#'
#'    For [ppm_to_icm()], [icm_to_ppm()], [pcm_to_ppm()],
#'    [ppm_to_pcm()], [ppm_to_pwm()], and [pwm_to_ppm()]: a `numeric`
#'    vector with length equal to input `numeric` vector.
#'
#'    For [score_match()]: a `numeric` vector with the match motif score.
#'
#'    For [summarise_motifs()]: a `data.frame` with columns representing
#'    the [universalmotif-class] slots.
#'
#' @examples
#' #######################################################################
#' ## Setting up some variables
#' data(examplemotif)
#' m <- normalize(examplemotif)
#' motif <- create_motif(nsites = 100, pseudocount = 0.8)["motif"]
#' motif.icm <- apply(motif, 2, ppm_to_icm, nsites = 100,
#'                    bkg = c(0.25, 0.25, 0.25, 0.25))
#' motif.ppm <- apply(motif.icm, 2, icm_to_ppm)
#' motif.consensus <- apply(motif.ppm, 2, get_consensus)
#' motif.aa <- create_motif(alphabet = "AA")["motif"]
#' motif.aa.consensus <- apply(motif.aa, 2, get_consensusAA, type = "PPM")
#' #######################################################################
#' 
#' #######################################################################
#' ## consensus_to_ppm
#' ## Do the opposite of get_consensus. Note that loss of information is
#' ## inevitable.
#' motif.ppm4 <- sapply(motif.consensus, consensus_to_ppm)
#'
#' #######################################################################
#' ## consensus_to_ppmAA
#' ## Do the opposite of get_consensusAA.
#' motif.aa2 <- sapply(motif.aa.consensus, consensus_to_ppmAA)
#'
#' #######################################################################
#' ## get_consensus
#' ## Get a consensus string from a DNA/RNA motif.
#' motif.consensus <- apply(motif.ppm, 2, get_consensus)
#'
#' #######################################################################
#' ## get_consensusAA
#' ## Get a consensus string from an amino acid motif. Unless each position
#' ## is clearly dominated by a single amino acid, the resulting string will
#' ## likely be useless.
#' motif.aa <- create_motif(alphabet = "AA")["motif"]
#' motif.aa.consensus <- apply(motif.aa, 2, get_consensusAA, type = "PPM")
#'
#' #######################################################################
#' ## get_match
#' ## Get all possible motif matches above input score
#' get_matches(m, 10)
#'
#' #######################################################################
#' ## icm_to_ppm
#' ## Do the opposite of ppm_to_icm.
#' motif.ppm <- apply(motif.icm, 2, icm_to_ppm)
#'
#' #######################################################################
#' ## make_DBscores
#' ## Generate P-value database for use with compare_motifs. Note that these
#' ## must be created individually for all combinations of methods and
#' ## normalisation.
#' \dontrun{
#' library(MotifDb)
#' motifs <- convert_motifs(MotifDb[1:100])
#' make_DBscores(motifs, method = "PCC")
#' }
#'
#' #######################################################################
#' ## motif_score
#' ## Calculate motif score from different thresholds
#' data(examplemotif)
#' m <- normalize(examplemotif)
#' motif_score(m, c(0, 0.8, 1))
#'
#' #######################################################################
#' ## pcm_to_ppm
#' ## Go from a count type motif to a probability type motif.
#' motif.pcm <- create_motif(type = "PCM", nsites = 50)["motif"]
#' motif.ppm2 <- apply(motif.pcm, 2, pcm_to_ppm, pseudocount = 1)
#'
#' #######################################################################
#' ## position_icscore
#' ## Similar to ppm_to_icm, except this calculates a sum for the position.
#' ic.scores <- apply(motif.ppm, 2, position_icscore, type = "PPM",
#'                    bkg = c(0.25, 0.25, 0.25, 0.25))
#'
#' #######################################################################
#' ## ppm_to_icm
#' ## Convert one column from a probability type motif to an information
#' ## content type motif.
#' motif <- create_motif(nsites = 100, pseudocount = 0.8)["motif"]
#' motif.icm <- apply(motif, 2, ppm_to_icm, nsites = 100,
#'                    bkg = c(0.25, 0.25, 0.25, 0.25))
#'
#' #######################################################################
#' ## ppm_to_pcm
#' ## Do the opposite of pcm_to_ppm.
#' motif.pcm2 <- apply(motif.ppm2, 2, ppm_to_pcm, nsites = 50)
#'
#' #######################################################################
#' ## ppm_to_pwm
#' ## Go from a probability type motif to a weight type motif.
#' motif.pwm <- apply(motif.ppm, 2, ppm_to_pwm, nsites = 100,
#'                    bkg = c(0.25, 0.25, 0.25, 0.25))
#'
#' #######################################################################
#' ## pwm_to_ppm
#' ## Do the opposite of ppm_to_pwm.
#' motif.ppm3 <- apply(motif.pwm, 2, pwm_to_ppm,
#'                     bkg = c(0.25, 0.25, 0.25, 0.25))
#'
#' #######################################################################
#' ## Note that not all type conversions can be done directly; for those
#' ## type conversions which are unavailable, universalmotif just chains
#' ## together others (i.e. from PCM -> ICM => pcm_to_ppm -> ppm_to_icm)
#'
#' #######################################################################
#' ## score_match
#' ## Calculate score of a particular match
#' score_match(m, "TATATAT")
#' score_match(m, "TATATAG")
#'
#' #######################################################################
#' ## summarise_motifs
#' ## Create a data.frame of information based on a list of motifs.
#' m1 <- create_motif()
#' m2 <- create_motif()
#' m3 <- create_motif()
#' summarise_motifs(list(m1, m2, m3))
#'
#' @seealso [create_motif()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name utils-motif
NULL

#' @rdname utils-motif
#' @export
motif_score <- function(motif, threshold = c(0, 1)) {

  if (any(threshold < 0) || any(threshold > 1))
    stop("For 'threshold', please only use values between 0 and 1")

  motif <- convert_motifs(motif)

  if (is.list(motif) && length(motif) > 1)
    stop("Please only input a single motif")
  else if (is.list(motif))
    motif <- motif[[1]]

  if (!is(motif, "universalmotif"))
    stop("Unknown motif object")

  motif <- convert_type_internal(motif, "PWM")

  if (any(is.infinite(motif@motif))) {
    warning("Found -Inf values in motif PWM, adding a pseudocount of 1",
            immediate. = TRUE)
    motif <- normalize(motif)
  }

  mat <- matrix(as.integer(motif@motif * 1000), nrow = nrow(motif@motif))

  s.max <- sum(apply(mat, 2, max))
  s.min <- sum(apply(mat, 2, min))

  s.total <- abs(s.max) + abs(s.min)

  out <- s.total * threshold - abs(s.min)
  names(out) <- paste0(as.character(threshold * 100), "%")

  out / 1000

}

#' @rdname utils-motif
#' @export
score_match <- function(motif, match) {

  if (missing(motif) || missing(match))
    stop("motif and/or match are missing")

  if (!is.character(match) || length(match) != 1)
    stop("match must be a single string")

  motif <- convert_motifs(motif)

  if (is.list(motif) && length(motif) != 1)
    stop("Please only input a single motif")
  else if (is.list(motif))
    motif <- motif[[1]]

  if (!is(motif, "universalmotif"))
    stop("Unknown motif class")

  if (motif@type != "PWM")
    motif <- convert_type_internal(motif, "PWM")

  if (any(is.infinite(motif@motif))) {
    warning("Found -Inf values in motif PWM, adding a pseudocount of 1",
            immediate. = TRUE)
    motif <- normalize(motif)
  }

  match <- safeExplode(match)
  if (length(match) != ncol(motif))
    stop(wmsg("motif length [", ncol(motif), "] and match length [",
              length(match), "] are not equal"))

  if (!all(match %in% rownames(motif@motif)))
    stop("Found letters in match not found in motif")

  score <- 0

  mat <- matrix(as.integer(motif@motif * 1000),
                nrow = nrow(motif@motif),
                dimnames = dimnames(motif@motif))

  for (i in seq_len(ncol(mat))) {
    score <- score + mat[match[i], i]
  }

  score / 1000

}

#' @rdname utils-motif
#' @export
get_matches <- function(motif, score) {

  motif <- convert_motifs(motif)
  if (motif@type != "PWM")
    motif <- convert_type_internal(motif, "PWM")
  if (any(is.infinite(motif@motif))) {
    warning(wmsg("found -Inf values in PWM motif, normalizing"),
            immediate. = TRUE)
    motif <- normalize(motif)
  }

  score.range <- motif_score(motif)
  if (score > score.range[2])
    stop(wmsg("input score is greater than max possible score ",
              round(score.range[2], 3)))
  if (score < score.range[1])
    stop(wmsg("input score is less than min possible score ",
              round(score.range[1], 3)))

  score <- as.integer(score * 1000)

  alph <- rownames(motif@motif)

  score.mat <- matrix(as.integer(motif@motif * 1000), nrow = nrow(motif@motif))

  alph.sort <- apply(score.mat, 2, order, decreasing = TRUE)
  for (i in seq_len(ncol(score.mat))) {
    score.mat[, i] <- score.mat[alph.sort[, i], i]
  }

  col.sort <- order(apply(score.mat, 2, max), decreasing = TRUE)
  score.mat <- score.mat[, col.sort]

  all.paths <- branch_and_bound_kmers(score.mat, score)

  all.paths <- all.paths[, order(col.sort), drop = FALSE]

  all.paths <- paths_alph_unsort(all.paths, alph.sort)

  paths_to_alph(all.paths, alph)

}

#' @rdname utils-motif
#' @export
ppm_to_icm <- function(position, bkg, schneider_correction = FALSE, nsites,
                       relative_entropy = FALSE) {
  # NOTE: Basic IC computation assumes uniform bkg frequencies!
  #       For different bkg frequencies: Relative entropy or Kullback-Leibler
  #       (KL) divergence
  if (is.null(bkg) || missing(bkg)) {
    bkg <- rep(1 / length(position), length(position))
  if (length(nsites) == 0) nsites <- 100
  }
  if (relative_entropy) {
    for (i in seq_along(position)) {
      position[i] <- position[i] * log2(position[i] / bkg[i])
      if (is.na(position[i]) || position[i] < 0) position[i] <- 0
    }
    position
  } else {
    height_after <- -sum(vapply(position, function(x) {
                                            y <- x * log2(x)
                                            ifelse(is.na(y), 0, y)
                                          }, numeric(1)))
    total_ic <- log2(length(position)) - height_after
    if (schneider_correction && !missing(nsites)) {
      correction <- ppm_to_pcm(position, nsites = nsites)
      if (requireNamespace("TFBSTools", quietly = TRUE)) {
        correction <- TFBSTools:::schneider_correction(matrix(correction), bkg)
      } else {
        stop("The 'TFBSTools' package is required for 'schneider_correction'")
      }
      total_ic <- total_ic + correction
    }
    ic <- position * total_ic
    ic
  }
}

#' @rdname utils-motif
#' @export
icm_to_ppm <- function(position) {
  total_ic <- sum(position)
  ppm <- position / total_ic
  ppm
}

#' @rdname utils-motif
#' @export
pcm_to_ppm <- function(position, pseudocount = 0.8) {
  possum <- sum(position)
  num_letters <- length(position)
  if (pseudocount != 0) {
    pos <- vapply(position, function(x)
                  (x + (pseudocount / num_letters)) / (possum + pseudocount),
                  double(1))
  } else {
    pos <- vapply(position, function(x) x / possum, double(1))
  }
  return(pos)
}

#' @rdname utils-motif
#' @export
ppm_to_pcm <- function(position, nsites = 100) {
  if (length(nsites) == 0 || missing(nsites)) nsites <- 100
  pos <- vapply(position, function(x) round(x * nsites), numeric(1))
  if (sum(pos) != nsites) {
    fix <- nsites - sum(pos)
    pos[which(range(pos)[2] == pos)[1]] <- pos[which(range(pos)[2] == pos)[1]] + fix
  }
  return(pos)
}

#' @rdname utils-motif
#' @export
ppm_to_pwm <- function(position, bkg, pseudocount = 0.8, nsites = 100,
                       smooth = TRUE) {
  if (missing(bkg)) bkg <- rep(1 / length(position), length(position))
  if (length(nsites) == 0) nsites <- 100
  if (smooth && pseudocount != 0) {
    position <- ppm_to_pcm(position, nsites = nsites)
    position <- pcm_to_ppm(position, pseudocount = pseudocount)
  }
  for (i in seq_along(position)) {
    position[i] <- log2(position[i] / bkg[i])
  }
  return(position)
}

#' @rdname utils-motif
#' @export
pwm_to_ppm <- function(position, bkg) {
  if (missing(bkg)) bkg <- rep(1 / length(position), length(position))
  position <- vapply(position, function(x) 2 ^ x, numeric(1))
  if (sum(position) > 0.99 && sum(position) < 1.01) return(position)
  for (i in seq_along(position)) position[i] <- position[i] * bkg[i]
  if (sum(position) > 0.99 && sum(position) < 1.01) return(position)
  warning("position does not add up to 1; normalizing..")
  pos_missing <- sum(position)
  position <- position / pos_missing
  position
}

#' @rdname utils-motif
#' @export
position_icscore <- function(position, bkg, type, pseudocount = 0.8, nsites = 100,
                             relative_entropy = FALSE) {

  motif <- position
  if (missing(bkg)) bkg <- rep(1 / length(position), length(position))

  if (length(nsites) == 0) nsites <- 100
  if (length(motif) != length(bkg)) {
    bkg <- rep(1 / length(motif), length(motif))
  }

  if (type == "PCM") {
      motif <- pcm_to_ppm(motif, pseudocount)
  }
  if (type == "PWM") motif <- pwm_to_ppm(motif, bkg = bkg)
  if (type == "ICM") return(sum(motif))

  if (type == "PPM") {
    motif <- ppm_to_pcm(position = motif, nsites = nsites)
    motif <- pcm_to_ppm(position = motif, pseudocount = pseudocount)
  }

  # ic <- vector(length = length(motif))
  #
  # for (i in seq_along(motif)) {
    # if (motif[i] == 0) ic[i] <- 0 else {
      # ic[i] <- motif[i] * log2(motif[i] / bkg[i])
    # }
  # }
#
  # sum(ic)

  if (relative_entropy) {
    for (i in seq_along(motif)) {
      motif[i] <- motif[i] * log2(motif[i] / bkg[i])
      if (is.na(motif[i]) || motif[i] < 0) motif[i] <- 0
    }
    total_ic <- sum(motif)
  } else {
    bkg <- rep(1 / length(position), length(position))
    height_after <- -sum(vapply(motif, function(x) {
                                         y <- x * log2(x)
                                         ifelse(is.na(y), 0, y)
                                       }, numeric(1)))
    total_ic <- log2(length(motif)) - height_after
  }
  total_ic

}

#' @rdname utils-motif
#' @export
get_consensus <- function(position, alphabet = "DNA", type = "PPM",
                          pseudocount = 0.8) {

  pos <- position

  if (type == "PCM") {
    pos <- pcm_to_ppm(pos, pseudocount)
    type <- "PPM"
  }

  if (type == "PWM") {
    pos <- pwm_to_ppm(pos)
    type <- "PPM"
  }

  if (type == "ICM") {
    warning("get_consensus cannot handle ICM type")
    return(character(0))
  }

  if (type == "PPM") {

    if (alphabet == "DNA") names(pos) <- DNA_BASES
    if (alphabet == "RNA") names(pos) <- DNA_BASES

    # single letter consensus:

    if (pos[1] > 0.5 && pos[1] > sort(pos)[3] * 2) return("A")
    if (pos[2] > 0.5 && pos[2] > sort(pos)[3] * 2) return("C")
    if (pos[3] > 0.5 && pos[3] > sort(pos)[3] * 2) return("G")
    if (pos[4] > 0.5 && pos[4] > sort(pos)[3] * 2) {
      ifelse(alphabet == "DNA", return("T"), return("U"))
    }

    # two letter consensus:

    if (pos[1] > 0.5) {
      if (names(sort(pos)[3]) == "C" && sort(pos)[3] > 0.25) return("M")
      if (names(sort(pos)[3]) == "G" && sort(pos)[3] > 0.25) return("R")
      if (names(sort(pos)[3]) %in% c("T", "U") &&
          sort(pos)[3] > 0.25) return("W")
    }

    if (pos[2] > 0.5) {
      if (names(sort(pos)[3]) == "A" && sort(pos)[3] > 0.25) return("M")
      if (names(sort(pos)[3]) == "G" && sort(pos)[3] > 0.25) return("S")
      if (names(sort(pos)[3]) %in% c("T", "U") &&
          sort(pos)[3] > 0.25) return("Y")
    }

    if (pos[3] > 0.5) {
      if (names(sort(pos)[3]) == "A" && sort(pos)[3] > 0.25) return("R")
      if (names(sort(pos)[3]) == "C" && sort(pos)[3] > 0.25) return("S")
      if (names(sort(pos)[3]) %in% c("T", "U") &&
          sort(pos)[3] > 0.25) return("K")
    }

    if (pos[4] > 0.5) {
      if (names(sort(pos)[3]) == "A" && sort(pos)[3] > 0.25) return("W")
      if (names(sort(pos)[3]) == "C" && sort(pos)[3] > 0.25) return("Y")
      if (names(sort(pos)[3]) == "G" && sort(pos)[3] > 0.25) return("K")
    }

    if ((pos[1] + pos[2]) > 0.75) return("M")
    if ((pos[1] + pos[3]) > 0.75) return("R")
    if ((pos[1] + pos[4]) > 0.75) return("W")

    if ((pos[2] + pos[3]) > 0.75) return("S")
    if ((pos[2] + pos[4]) > 0.75) return("Y")

    if ((pos[3] + pos[4]) > 0.75) return("K")

    # three letter consensus:

    if (all(pos[c(1, 2, 4)] > 0.25)) return("H")
    if (all(pos[c(2, 3, 4)] > 0.25)) return("B")
    if (all(pos[c(1, 2, 3)] > 0.25)) return("V")
    if (all(pos[c(1, 3, 4)] > 0.25)) return("D")

    # no consensus:

    return("N")

  }

}

#' @rdname utils-motif
#' @export
consensus_to_ppm <- function(letter) {
  if (letter == "A") return(c(0.997, 0.001, 0.001, 0.001))
  if (letter == "C") return(c(0.001, 0.997, 0.001, 0.001))
  if (letter == "G") return(c(0.001, 0.001, 0.997, 0.001))
  if (letter %in% c("T", "U")) return(c(0.001, 0.001, 0.001, 0.997))
  if (letter == "R") return(c(0.499, 0.001, 0.499, 0.001))
  if (letter == "Y") return(c(0.001, 0.499, 0.001, 0.499))
  if (letter == "M") return(c(0.499, 0.499, 0.001, 0.001))
  if (letter == "K") return(c(0.001, 0.001, 0.499, 0.499))
  if (letter == "S") return(c(0.001, 0.499, 0.499, 0.001))
  if (letter == "W") return(c(0.499, 0.001, 0.001, 0.499))
  if (letter == "H") return(c(0.333, 0.333, 0.001, 0.333))
  if (letter == "B") return(c(0.001, 0.333, 0.333, 0.333))
  if (letter == "V") return(c(0.333, 0.333, 0.333, 0.001))
  if (letter == "D") return(c(0.333, 0.001, 0.333, 0.333))
  if (letter %in% c("N", "+", "-", ".")) return(c(0.25, 0.25, 0.25, 0.25))
  stop(letter, " is not an IUPAC symbol")
}

#' @rdname utils-motif
#' @export
consensus_to_ppmAA <- function(letter) {
  if (letter %in% c("X", ".", "-", "+")) return(rep(0.05, 20))
  if (letter == "B") return(c(rep(0.001, 2), 0.491, rep(0.001, 8), 0.491,
                              rep(0.001, 8)))
  if (letter == "Z") return(c(rep(0.001, 3), 0.491, rep(0.001, 9), 0.491,
                              rep(0.001, 6)))
  if (letter == "J") return(c(rep(0.001, 7), 0.491, 0.001, 0.491,
                              rep(0.001, 10)))
  i <- which(AA_STANDARD == letter)
  c(rep(0.001, i - 1), 0.981, rep(0.001, 20 - i))
}

#' @rdname utils-motif
#' @export
get_consensusAA <- function(position, type, pseudocount) {
  motif <- position
  if (type == "PCM") {
    motif <- pcm_to_ppm(motif, pseudocount)
    type <- "PPM"
  }
  if (type == "PWM") {
    motif <- pwm_to_ppm(motif)
    type <- "PPM"
  }
  if (type == "ICM") {
    warning("get_consensusAA cannot handle ICM type")
    return(character(0))
  }
  if (motif[3] >= 0.4 && motif[12] >= 0.4) return("B")
  if (motif[4] >= 0.4 && motif[14] >= 0.4) return("Z")
  if (motif[8] >= 0.4 && motif[10] >= 0.4) return("J")
  if (all(motif == 0.05) || all(motif < 0.1)) return("X")
  .aa <- order(motif, decreasing = TRUE)
  if (motif[.aa[1]] == motif[.aa[2]]) return("X")
  if (motif[.aa[1]] > 0.2) return(AA_STANDARD[.aa[1]])
  "X"
}

#' @rdname utils-motif
#' @export
summarise_motifs <- function(motifs, na.rm = TRUE) {
  # ~0.05 seconds for entire MotifDb library
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  classcheck <- vapply(motifs, function(x) !is(x, "universalmotif"), logical(1))
  if (any(classcheck)) stop("all motifs must be 'universalmotif'")
  out <- summarise_motifs_cpp(motifs)
  out <- out[, c("name", "altname", "family", "organism", "consensus", "alphabet",
                 "strand", "icscore", "nsites", "bkgsites", "pval", "qval", "eval")]
  if (na.rm) out <- Filter(function(x) !all(is.na(x)), out)
  out
}
