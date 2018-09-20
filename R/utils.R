#' Utility functions.
#'
#' Various small functions used for motif creation.
#'
#' @param position `numeric` A numeric vector representing the frequency or
#'    probability for each alphabet letter at a specific position.
#' @param bkg `Numeric` Should be the same length as the alphabet length.
#' @param nsites `numeric(1)` Number of sites motif originated from.
#' @param schneider_correction `logical(1)` Apply sample size correction.
#' @param relative_entropy `logical(1)` Calculate information content as
#'    relative entropy or Kullback-Leibler divergence.
#' @param pseudocount `numeric(1)` Used to prevent zeroes in motif matrix.
#' @param smooth `logical(1)` Apply pseudocount correction.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM' 'ICM')`.
#' @param alphabet `character(1)` One of `c('DNA', 'RNA')`.
#' @param letter `character(1)` Any DNA, RNA, or AA IUPAC letter. Ambiguity letters
#'    are accepted.
#' @param db.motifs `list` Database motifs.
#' @param method `character(1)` One of `c('Pearson', 'Euclidean', 'KL')`.
#' @param shuffle.db `logical(1)` Shuffle `db.motifs` rather than
#'    generate random motifs with [create_motif()].
#' @param shuffle.k `numeric(1)`
#' @param shuffle.method `character(1)`
#' @param shuffle.leftovers `character(1)`
#' @param rand.tries `numeric(1)`
#' @param normalise.scores `logical(1)`
#' @param min.overlap `numeric(1)` Minimum required motif overlap.
#' @param min.mean.ic `numeric(1)`
#' @param progress `logical(1)` Show progress.
#' @param BP `logical(1)` Use BiocParallel.
#' @param motifs `list` A list of \linkS4class{universalmotif} motifs.
#' @param na.rm `logical` Remove columns where all values are \code{NA}.
#'
#' @return
#'    For [ppm_to_icm()], [icm_to_ppm()], [pcm_to_ppm()],
#'    [ppm_to_pcm()], [ppm_to_pwm()], and [pwm_to_ppm()]: a `numeric`
#'    vector with length equal to input `numeric` vector.
#'
#'    For [consensus_to_ppm()] and [consensus_to_ppmAA()]: a numeric
#'    vector of length 4 and 20, respectively.
#'
#'    For [position_icscore()]: a `numeric` vector of length 1.
#'
#'    For [get_consensus()] and [get_consensusAA()]: a character vector
#'    of length 1.
#'
#'    For [make_DBscores()]: a `data.frame` with score distributions for the
#'    input database.
#'
#'    For [summarise_motifs()]: a `data.frame` with columns representing
#'    the [universalmotif-class] slots.
#'
#' @seealso [create_motif()], [compare_motifs()]
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name utilities
NULL

DNA_DI <- c("AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT")

# DNA_TRI <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG",
             # "AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC",
             # "CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA",
             # "GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT",
             # "GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG",
             # "TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT")

#' @rdname utilities
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

#' @rdname utilities
#' @export
icm_to_ppm <- function(position) {
  total_ic <- sum(position)
  ppm <- position / total_ic
  ppm
}

#' @rdname utilities
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

#' @rdname utilities
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

#' @rdname utilities
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

#' @rdname utilities
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

#' @rdname utilities
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

#' @rdname utilities
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

#' @rdname utilities
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

#' @rdname utilities
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

#' @rdname utilities
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

#' @rdname utilities
#' @export
summarise_motifs <- function(motifs, na.rm = TRUE) {
  if (!is.list(motifs) || !is(motifs, "MotifList")) motifs <- list(motifs)
  classcheck <- vapply(motifs, function(x) !is(x, "universalmotif"), logical(1))
  if (any(classcheck)) stop("all motifs must be 'universalmotif'")
  out <- do.call(rbind, lapply(motifs, as.data.frame))
  out <- out[, c("name", "altname", "family", "organism", "consensus", "alphabet",
                 "strand", "icscore", "nsites", "bkgsites", "pval", "qval", "eval")]
  if (na.rm) out <- Filter(function(x) !all(is.na(x)), out)
  out
}

.internal_convert <- function(motifs, class = NULL) {

  if (is.null(class)) {
    CLASS_PKG <- attributes(class(motifs))$package
    CLASS <- class(motifs)
    CLASS_IN <- paste(CLASS_PKG, CLASS, sep = "-")
    return(CLASS_IN)
  } else {
    if (length(class) == 1 && class[1] != "universalmotif-universalmotif") {
      tryCatch(motifs <- convert_motifs(motifs, class),
               error = function(e) message("motifs converted to class 'universalmotif'"))
    } else if (length(class) > 1) message("motifs converted to class 'universalmotif'")
    return(motifs)
  }

}

# for a motif of length 4, the transition matrix is something like this:

#       bkg pos1 pos2 pos3 pos4
#  bkg    0    1    0    0    0
# pos1    0    0    1    0    0
# pos2    0    0    0    1    0
# pos3    0    0    0    0    1
# pos4    1    0    0    0    0

lapply_ <- function(X, FUN, ..., BP = FALSE, PB = FALSE) {

  FUN <- match.fun(FUN)

  if (!BP) {

    if (!PB) {

      out <- lapply(X, FUN, ...)

    } else {

      out <- vector("list", length(X))
      max <- length(X)
      print_pb(0)
      if (is.list(X)) {
        for (i in seq_along(X)) {
          out[[i]] <- do.call(FUN, list(X[[i]], ...))
          update_pb(i, max)
        }
      } else {
        for (i in seq_along(X)) {
          out[[i]] <- do.call(FUN, list(X[i], ...))
          update_pb(i, max)
        }
      }

    }

  } else {

    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      out <- BiocParallel::bplapply(X, FUN, ...)
    } else {
      stop("'BiocParallel' is not installed")
    }
    # BPPARAM <- BiocParallel::bpparam()
    # if (PB) BPPARAM$progressbar <- TRUE
    # out <- BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)

  }

  out

}

mapply_ <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE,
                    USE.NAMES = TRUE, BP = FALSE, PB = FALSE) {

  FUN <- match.fun(FUN)

  if (!BP) {

    if (!PB) {

      out <- mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY,
                    USE.NAMES = USE.NAMES)

    } else {

      # not sure how to implement USE.NAMES here, get error sometimes
      dots <- list(...)
      dots.len <- vapply(dots, length, numeric(1))
      dots.len.max <- max(dots.len)
      dots <- lapply(dots, rep, length.out = dots.len.max)
      out <- vector("list", dots.len.max)

      print_pb(0)
      for (i in seq_len(dots.len.max)) {
        dots.i <- mapply(function(dots, i) {
                           if (is.list(dots)) dots[[i]]
                           else dots[i]
                    }, dots, i, SIMPLIFY = FALSE)
        out[[i]] <- do.call(FUN, c(dots.i, MoreArgs))
        update_pb(i, dots.len.max)
      }

      if (SIMPLIFY && length(dots))
        out <- simplify2array(out, higher = (SIMPLIFY == "array"))

    }

  } else {

    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      BPPARAM <- BiocParallel::bpparam()
      if (PB) BPPARAM$progressbar <- TRUE
      out <- BiocParallel::bpmapply(FUN, ..., MoreArgs = MoreArgs,
                                    SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES,
                                    BPPARAM = BPPARAM)
    } else {
      stop("'BiocParallel' is not installed")
    }

  }

  out

}
