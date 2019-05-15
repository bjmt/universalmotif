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
#' @param method `character` Any of `c('PCC', 'MPCC', 'EUCL', 'MEUCL',
#'    'SW', 'MSW', 'KL', 'MKL')`. See [compare_motifs()]. If multiple methods
#'    are provided, each result will be in its own list.
#' @param min.mean.ic `numeric(1)` See [compare_motifs()].
#' @param min.position.ic `numeric(1)` Minimum information content required between
#'    individual alignment positions for it to be counted in the final alignment
#'    score. It is recommended to use this together with `normalise.scores = TRUE`,
#'    as this will help punish scores resulting from only a fraction of an
#'    alignment.
#' @param min.overlap `numeric(1)` Minimum required motif overlap. See
#'    [compare_motifs()].
#' @param motif Motif object to calculate scores from.
#' @param motifs `list` A list of \linkS4class{universalmotif} motifs.
#' @param na.rm `logical` Remove columns where all values are `NA`.
#' @param normalise.scores `logical(1)` See [compare_motifs()].
#' @param nsites `numeric(1)` Number of sites motif originated from.
#' @param nthreads `numeric(1)` Run [compare_motifs()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
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
#' @param shuffle.method `character(1)` See [shuffle_motifs()].
#' @param smooth `logical(1)` Apply pseudocount correction.
#' @param threshold `numeric(1)` Any number of numeric values between 0 and 1
#'    representing score percentage.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM' 'ICM')`.
#' @param widths `numeric` Motif widths to use in P-value database calculation.
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
#'    input database, or a `list` with a `data.frame` for each method and
#'    an additional `list` entry logging function parameters.
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

  all.paths <- branch_and_bound_cpp_exposed(score.mat, score)

  all.paths <- all.paths[, order(col.sort), drop = FALSE]

  all.paths <- paths_alph_unsort(all.paths, alph.sort - 1)

  paths_to_alph(all.paths, alph)

}

#' @rdname utils-motif
#' @export
icm_to_ppm <- function(position) {
  icm_to_ppmC(position)
}

#' @rdname utils-motif
#' @export
make_DBscores <- function(db.motifs, method, shuffle.db = TRUE,
                          shuffle.k = 3, shuffle.method = "linear",
                          rand.tries = 100, widths = 5:30,
                          min.position.ic = 0,
                          normalise.scores = TRUE, min.overlap = 6,
                          min.mean.ic = 0.25, progress = TRUE, BP = FALSE,
                          nthreads = 1) {

  if (length(method) > 1) {
    out <- vector("list", length(method) + 1)
    names(out) <- c(method, "args")
    for (m in method) {
      out[[m]] <- make_DBscores(db.motifs, m, shuffle.db, shuffle.k, shuffle.method,
                                rand.tries, widths, min.position.ic,
                                normalise.scores, min.overlap, min.mean.ic,
                                progress, BP, nthreads)
    }
    out$args <- args[-1]
    return(out)
  }

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(method = args$method,
                                      shuffle.method = args$shuffle.method),
                                 c(0, 1), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(shuffle.k = args$shuffle.k,
                                     rand.tries = args$rand.tries,
                                     min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     nthreads = args$nthreads,
                                     min.position.ic = args$min.position.ic),
                                numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(shuffle.db = args$shuffle.db,
                                      progress = args$progress, BP = args$BP,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  db.motifs <- convert_motifs(db.motifs)
  db.motifs <- convert_type(db.motifs, "PPM")

  db.ncols <- vapply(db.motifs, function(x) ncol(x@motif), numeric(1))

  rand.mots <- vector("list", length(widths))
  names(rand.mots) <- as.character(widths)
  rand.ncols <- rep(FALSE, length(widths))
  names(rand.ncols) <- as.character(widths)
  if (shuffle.db) {

    for (i in widths) {
      if (!any(db.ncols == i)) next;
      tmp <- list()
      while (length(tmp) < rand.tries) {
        tmp <- c(tmp, shuffle_motifs(db.motifs[db.ncols == i], k = shuffle.k,
                                     method = shuffle.method))
      }
      tmp <- tmp[seq_len(rand.tries)]
      rand.mots[[as.character(i)]] <- tmp
      rand.ncols[as.character(i)] <- TRUE
    }

  } else {
    for (i in widths) {
      rand.mots[[as.character(i)]] <- lapply(rand.tries, function(x) create_motif(i))
    }
  }

  totry <- expand.grid(list(target = as.numeric(widths[rand.ncols]),
                            subject = sort(unique(db.ncols))))[, 2:1]
  totry$mean <- rep(NA_real_, nrow(totry))
  totry$sd <- rep(NA_real_, nrow(totry))

  res <- vector("list", length(unique(db.ncols)))
  allcomps <- vector("list", length(unique(db.ncols)))
  names(res) <- as.character(sort(unique(db.ncols)))
  names(allcomps) <- as.character(sort(unique(db.ncols)))

  rand.mots <- unlist(rand.mots, recursive = FALSE)
  randcols <- vapply(rand.mots, function(x) ncol(x@motif), integer(1))

  if (progress) print_pb(0)
  counter <- 1
  total <- length(unique(db.ncols))

  for (i in sort(unique(db.ncols))) {

    tmp <- db.motifs[db.ncols == i]

    tmpall <- c(tmp, rand.mots)
    mtmp <- lapply(tmpall, function(x) x@motif)
    btmp <- lapply(tmpall, function(x) x@bkg)
    comps <- get_comp_indices(seq_along(tmp), length(tmpall))
    comps <- comps[!comps[, 2] %in% seq_along(tmp), ]
    allcomps[[as.character(i)]] <- rep(randcols, length(tmp))
    res[[as.character(i)]] <- compare_motifs_cpp(mtmp, comps[, 1] - 1, comps[, 2] - 1,
                                                 method, min.overlap, FALSE,
                                                 btmp, 1, FALSE, min.mean.ic,
                                                 normalise.scores, nthreads,
                                                 min.position.ic)

    if (progress) update_pb(counter, total)
    counter <- counter + 1

  }

  for (i in seq_len(nrow(totry))) {
    tmp <- res[[as.character(totry[i, 1])]]
    totry$mean[i] <- mean(tmp[allcomps[[as.character(totry[i, 1])]] == totry[i, 2]])
    totry$sd[i] <- sd(tmp[allcomps[[as.character(totry[i, 1])]] == totry[i, 2]])
  }

  totry$method <- rep(method, nrow(totry))
  totry$normalised <- rep(normalise.scores, nrow(totry))
  totry

}

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
    warning("Found -Inf values in motif PWM, adding a pseudocount",
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
pcm_to_ppm <- function(position, pseudocount = 0) {
  pcm_to_ppmC(position, pseudocount)
}

#' @rdname utils-motif
#' @export
position_icscore <- function(position, bkg = 0, type = "PPM", pseudocount = 1,
                             nsites = 100, relative_entropy = FALSE) {

  if (!type %in% c("ICM", "PPM", "PWM", "PCM"))
    stop("type must be one of ICM, PCM, PPM, PWM")

  position_icscoreC(position, bkg, type, pseudocount, nsites, relative_entropy)

}

#' @rdname utils-motif
#' @export
ppm_to_icm <- function(position, bkg, schneider_correction = FALSE, nsites = 100,
                       relative_entropy = FALSE) {
  # NOTE: Basic IC computation assumes uniform bkg frequencies!
  #       For different bkg frequencies: Relative entropy or Kullback-Leibler
  #       (KL) divergence
  if (is.null(bkg) || missing(bkg)) {
    bkg <- rep(1 / length(position), length(position))
  }
  if (relative_entropy) {
    ppm_to_icmC(position, bkg, TRUE)
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
ppm_to_pcm <- function(position, nsites = 100) {
  ppm_to_pcmC(position, nsites)
}

#' @rdname utils-motif
#' @export
ppm_to_pwm <- function(position, bkg, pseudocount = 1, nsites = 100,
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
