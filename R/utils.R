#' Utility functions.
#'
#' Various small functions used for motif creation. Note that all of the
#' functions described here (with the following exceptions: [make_DBscores()],
#' [summarise_motifs()]) are here only for demonstration purposes; internally
#' the \pkg{universalmotif} package uses faster C++ code for type conversion
#' and consensus manipulation.
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
#' @param method `character(1)` One of `c('PCC', 'MPCC', 'EUCL', 'MEUCL',
#'    'SW', 'MSW', 'KL', 'MKL')`. See [compare_motifs()].
#' @param shuffle.db `logical(1)` Shuffle `db.motifs` rather than
#'    generate random motifs with [create_motif()].
#' @param shuffle.k `numeric(1)` See [shuffle_motifs()].
#' @param shuffle.method `character(1)` See [shuffle_motifs()].
#' @param shuffle.leftovers `character(1)` See [shuffle_motifs()].
#' @param rand.tries `numeric(1)` Number of random motifs to create for
#'    P-value computation.
#' @param normalise.scores `logical(1)` See [compare_motifs()].
#' @param min.overlap `numeric(1)` Minimum required motif overlap. See
#'    [compare_motifs()].
#' @param min.mean.ic `numeric(1)` See [compare_motifs()].
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param BP `logical(1)` Use the \pkg{BiocParallel} package. See
#'    [BiocParallel::register()] to change the default backend.
#' @param motifs `list` A list of \linkS4class{universalmotif} motifs.
#' @param na.rm `logical` Remove columns where all values are `NA`.
#' @param lets `character` Letters to generate k-lets from.
#' @param k `integer(1)` K-let size.
#' @param motif Motif object to calculate scores from.
#' @param threshold `numeric(1)` Any number of numeric values between 0 and 1
#'    representing score percentage.
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
#' @examples
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
#' ## ppm_to_icm
#' ## Convert one column from a probability type motif to an information
#' ## content type motif.
#' motif <- create_motif(nsites = 100, pseudocount = 0.8)["motif"]
#' motif.icm <- apply(motif, 2, ppm_to_icm, nsites = 100,
#'                    bkg = c(0.25, 0.25, 0.25, 0.25))
#'
#' ## icm_to_ppm
#' ## Do the opposite of ppm_to_icm.
#' motif.ppm <- apply(motif.icm, 2, icm_to_ppm)
#'
#' ## pcm_to_ppm
#' ## Go from a count type motif to a probability type motif.
#' motif.pcm <- create_motif(type = "PCM", nsites = 50)["motif"]
#' motif.ppm2 <- apply(motif.pcm, 2, pcm_to_ppm, pseudocount = 1)
#'
#' ## ppm_to_pcm
#' ## Do the opposite of pcm_to_ppm.
#' motif.pcm2 <- apply(motif.ppm2, 2, ppm_to_pcm, nsites = 50)
#'
#' ## ppm_to_pwm
#' ## Go from a probability type motif to a weight type motif.
#' motif.pwm <- apply(motif.ppm, 2, ppm_to_pwm, nsites = 100,
#'                    bkg = c(0.25, 0.25, 0.25, 0.25))
#'
#' ## pwm_to_ppm
#' ## Do the opposite of ppm_to_pwm.
#' motif.ppm3 <- apply(motif.pwm, 2, pwm_to_ppm,
#'                     bkg = c(0.25, 0.25, 0.25, 0.25))
#'
#' ## Note that not all type conversions can be done directly; for those
#' ## type conversions which are unavailable, universalmotif just chains
#' ## together others (i.e. from PCM -> ICM => pcm_to_ppm -> ppm_to_icm)
#'
#' ## position_icscore
#' ## Similar to ppm_to_icm, except this calculates a sum for the position.
#' ic.scores <- apply(motif.ppm, 2, position_icscore, type = "PPM",
#'                    bkg = c(0.25, 0.25, 0.25, 0.25))
#'
#' ## get_consensus
#' ## Get a consensus string from a DNA/RNA motif.
#' motif.consensus <- apply(motif.ppm, 2, get_consensus)
#'
#' ## consensus_to_ppm
#' ## Do the opposite of get_consensus. Note that loss of information is
#' ## inevitable.
#' motif.ppm4 <- sapply(motif.consensus, consensus_to_ppm)
#'
#' ## get_consensusAA
#' ## Get a consensus string from an amino acid motif. Unless each position
#' ## is clearly dominated by a single amino acid, the resulting string will
#' ## likely be useless.
#' motif.aa <- create_motif(alphabet = "AA")["motif"]
#' motif.aa.consensus <- apply(motif.aa, 2, get_consensusAA, type = "PPM")
#'
#' ## consensus_to_ppmAA
#' ## Do the opposite of get_consensusAA.
#' motif.aa2 <- sapply(motif.aa.consensus, consensus_to_ppmAA)
#'
#' ## summarise_motifs
#' ## Create a data.frame of information based on a list of motifs.
#' m1 <- create_motif()
#' m2 <- create_motif()
#' m3 <- create_motif()
#' summarise_motifs(list(m1, m2, m3))
#'
#' ## get_klets
#' ## Generate all possible k-lets for a set of letters
#' get_klets(c("A", "C", "G", "T"), 3)
#'
#' ## motif_score
#' ## Calculate motif score from different thresholds
#' m <- create_motif()
#' motif_score(m, c(0, 0.8, 1))
#'
#' @seealso [create_motif()], [compare_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name utilities
NULL

# INTERNAL CONSTANTS -----------------------------------------------------------

DNA_DI <- c("AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT")

DNA_TRI <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG",
             "AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC",
             "CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA",
             "GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT",
             "GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG",
             "TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT")

# TYPES:
#
#        NULL   0
#      symbol   1
#     closure   3
# environment   4
#    promises   5
#    logical   10
#    integer   13
#       real   14
#    complex   15
#     string   16
#        ...   17
#        ANY   18
#    generic   19
# expression   20
#  byte code   21
#    ext ptr   22
#        raw   24
#         S4   25

TYPE_LOGI <- 10L
TYPE_NUM  <- 14L
TYPE_CHAR <- 16L
TYPE_S4   <- 25L

UNIVERSALMOTIF_SLOTS <- c(

  "name",
  "altname",
  "family",
  "organism",
  "motif",
  "alphabet",
  "type",
  "icscore",
  "nsites",
  "pseudocount",
  "bkg",
  "bkgsites",
  "consensus",
  "strand",
  "pval",
  "qval",
  "eval",
  "multifreq",
  "extrainfo"

)

COMPARE_METRICS <- c("PCC", "MPCC", "EUCL", "MEUCL", "SW", "MSW", "KL", "MKL")

# EXPORTED UTILITIES -----------------------------------------------------------

#' @rdname utilities
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
    warning("Found -Inf values in motif PWM, adding a pseudocount of 1")
    motif <- normalize(motif)
  }

  s.max <- sum(apply(motif@motif, 2, max))
  s.min <- sum(apply(motif@motif, 2, min))

  s.total <- abs(s.max) + abs(s.min)

  out <- s.total * threshold - abs(s.min)
  names(out) <- paste0(as.character(threshold * 100), "%")

  out

}

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

# INTERNAL UTILITIES -----------------------------------------------------------

.internal_convert <- function(motifs, class = NULL) {

  if (is.null(class)) {

    CLASS <- class(motifs)
    CLASS_PKG <- attributes(CLASS)$package
    CLASS_IN <- collapse_cpp(c(CLASS_PKG, "-", CLASS))

    CLASS_IN

  } else {

    if (length(class) == 1 && class[1] != "universalmotif-universalmotif") {

      tryCatch(motifs <- convert_motifs(motifs, class),
               error = function(e) message("motifs converted to class 'universalmotif'"))

    } else if (length(class) > 1)
      message("motifs converted to class 'universalmotif'")

    motifs

  }

}

# for a motif of length 4, the transition matrix is something like this:
#       bkg pos1 pos2 pos3 pos4
#  bkg    0    1    0    0    0
# pos1    0    0    1    0    0
# pos2    0    0    0    1    0
# pos3    0    0    0    0    1
# pos4    1    0    0    0    0

# Inspired from S4Vectors::wmsg
wmsg2 <- function(..., exdent = 0, indent = 0)
  paste0(strwrap(paste0(..., collapse = ""), exdent = exdent, indent = indent),
         collapse = "\n")

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
