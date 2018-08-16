#' Utility functions.
#'
#' Various small functions used for motif creation.
#'
#' @param position \code{numeric} A numeric vector representing the frequency or
#'    probability for each alphabet letter at a specific position.
#' @param bkg \code{Numeric} Should be the same length as the alphabet length.
#' @param nsites \code{numeric(1)} Number of sites motif originated from.
#' @param schneider_correction \code{logical(1)} Apply sample size correction.
#' @param relative_entropy \code{logical(1)} Calculate information content as
#'    relative entropy or Kullback-Leibler divergence.
#' @param pseudocount \code{numeric(1)} Used to prevent zeroes in motif matrix.
#' @param smooth \code{logical(1)} Apply pseudocount correction.
#' @param type \code{character(1)} One of \code{c('PCM', 'PPM', 'PWM' 'ICM')}.
#' @param alphabet \code{character(1)} One of \code{c('DNA', 'RNA')}.
#' @param letter \code{character(1)} Any DNA, RNA, or AA IUPAC letter. Ambiguity letters
#'    are accepted.
#' @param db.motifs \code{list} Database motifs.
#' @param method \code{character(1)} One of \code{c('Pearson', 'Euclidean', 'KL')}.
#' @param shuffle.db \code{logical(1)} Shuffle \code{db.motifs} rather than
#'    generate random motifs with \code{\link{create_motif}}
#' @param shuffle.k \code{numeric(1)}
#' @param shuffle.method \code{character(1)}
#' @param shuffle.leftovers \code{character(1)}
#' @param rand.tries \code{numeric(1)}
#' @param min.overlap \code{numeric(1)} Minimum required motif overlap.
#' @param min.mean.ic \code{numeric(1)}
#' @param progress_bar \code{logical(1)} Show progress.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return 
#'    For \code{ppm_to_icm}, \code{icm_to_ppm}, \code{pcm_to_ppm},
#'    \code{ppm_to_pcm}, \code{ppm_to_pwm}, and \code{pwm_to_ppm}: a numeric
#'    vector with length equal to input numeric vector.
#'
#'    For \code{consensus_to_ppm} and \code{consensus_to_ppmAA}: a numeric
#'    vector of length 4 and 20, respectively.
#'
#'    For \code{position_icscore}: a numeric vector of length 1.
#'
#'    For \code{get_consensus} and \code{get_consensusAA}: a character vector
#'    of length 1.
#'
#'    For \code{make_DBscores}: a data.frame with score distributions for the
#'    input database.
#'
#' @seealso \code{\link{create_motif}}
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
ppm_to_icm <- function(position, bkg,
                       schneider_correction = FALSE, nsites,
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
      correction <- TFBSTools:::schneider_correction(matrix(correction), bkg)
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
position_icscore <- function(position, bkg, type,
                             pseudocount = 0.8, nsites = 100,
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

.internal_convert <- function(motifs, class = NULL, BPPARAM = SerialParam()) {

  if (is.null(class)) {
    CLASS_PKG <- attributes(class(motifs))$package
    CLASS <- class(motifs)
    CLASS_IN <- paste(CLASS_PKG, CLASS, sep = "-")
    return(CLASS_IN)
  } else {
    if (length(class) == 1 && class != "universalmotif-universalmotif") {
      tryCatch(motifs <- convert_motifs(motifs, class, BPPARAM = BPPARAM),
               error = function(e) message("motifs converted to class 'universalmotif'"))
    }
    return(motifs)
  }
  
}

add_multi <- function(bkg, sequences, k) {

  seq.width <- unique(width(sequences))
  if (seq.width < k - 1) {
    warning("motif is not long enough for k = ", k)
    return(matrix())
  }

  emissions <- matrix(nrow = 4^k, ncol = seq.width - k + 1)
  
  multi_rows <- matrix(nrow = k, ncol = 4^k)
  for (i in seq_len(k)) {
    j <- rep(DNA_BASES, each = 4^(k - i + 1) / 4)
    if (length(j) != 4^k) j <- rep(j, 4^k / length(j))
    multi_rows[i, ] <- j
  }
  multi_rows <- apply(multi_rows, 2, paste, collapse = "")

  rownames(emissions) <- multi_rows
  colnames(emissions) <- seq_len(ncol(emissions))

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])
  seqs.split <- t(seqs.split)

  for (i in seq_len(seq.width - k + 1)) {
    current.seqs <- seqs.split[, i:(i + k - 1)]
    current.seqs <- apply(current.seqs, 1, paste, collapse = "")
    current.seqs <- DNAStringSet(current.seqs)
    emissions.i <- colSums(oligonucleotideFrequency(current.seqs, k, 1))
    emissions.i <- emissions.i / sum(emissions.i)
    emissions[, i] <- emissions.i
  }

  emissions

}

add_multi_ANY <- function(sequences, k, alph) {

  seq.width <- unique(width(sequences))
  if (seq.width < k - 1) {
    warning("motif is not long enough for k = ", k)
    return(matrix())
  }

  alph.len <- length(alph)
  emissions <- matrix(rep(0, alph.len^k * (seq.width - k + 1)),
                      nrow = alph.len^k, ncol = seq.width - k + 1)

  alph.comb <- as.matrix(expand.grid(rep(list(alph), k)))
  alph.comb <- apply(alph.comb, 1, paste, collapse = "")
  alph.comb <- sort(alph.comb)

  rownames(emissions) <- alph.comb
  colnames(emissions) <- seq_len(ncol(emissions))

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])

  seq.list <- lapply(seq_len(ncol(seqs.split)),
                     function(x) single_to_k(seqs.split[, x], k))

  seq.list.i <- character(length(seq.list[[1]]))
  for (i in seq_along(seq.list[[1]])) {
    for (j in seq_along(seq.list)) {
      seq.list.i[j] <- seq.list[[j]][i]
    }
    seq.list.t <- as.matrix(table(seq.list.i))
    for (j in seq_len(nrow(seq.list.t))) {
      emissions[rownames(emissions) == rownames(seq.list.t)[j], i] <- seq.list.t[j, ]
    }
    emissions[, i] <- emissions[, i] / sum(emissions[, i])
  }

  emissions

}

# for a motif of length 4, the transition matrix is something like this:

#       bkg pos1 pos2 pos3 pos4 
#  bkg    0    1    0    0    0
# pos1    0    0    1    0    0
# pos2    0    0    0    1    0
# pos3    0    0    0    0    1
# pos4    1    0    0    0    0

check_input_params <- function(char, charlen, num, numlen, logi, logilen) {

  fails <- c()

  if (!missing(char)) {

    if (missing(charlen)) charlen <- rep(1, length(char))

    char.check1 <- vapply(char, is.character, logical(1))
    char.check2 <- mapply(function(x, y) {
                            if (y == 0) TRUE else length(x) == y
                     }, char, charlen)

    if (any(!char.check1)) {
      j <- length(fails) + 1
      for (i in which(!char.check1)) {
        fails[j] <- paste0(" * Incorrect type for: `", names(char)[i],
                          "`; expected `character`, got `",
                          mode(char[[i]]), "`\n")
        j <- j + 1
      }
    }

    if (any(!char.check2)) {
      j <- length(fails) + 1
      for (i in which(!char.check2)) {
        fails[j] <- paste0(" * Incorrect vector length for: `", names(char)[i],
                          "`; expected length ", charlen[i], ", got ",
                          length(char[[i]]), "\n")
        j <- j + 1
      }
    }

  }

  if (!missing(num)) {

    if (missing(numlen)) numlen <- rep(1, length(num))

    num.check1 <- vapply(num, is.numeric, logical(1))
    num.check2 <- mapply(function(x, y) {
                            if (y == 0) TRUE else length(x) == y
                     }, num, numlen)

    if (any(!num.check1)) {
      j <- length(fails) + 1
      for (i in which(!num.check1)) {
        fails[j] <- paste0(" * Incorrect type for: `", names(num)[i],
                          "`; expected `numeric`, got `",
                          mode(num[[i]]), "`\n")
        j <- j + 1
      }
    }

    if (any(!num.check2)) {
      j <- length(fails) + 1
      for (i in which(!num.check2)) {
        fails[j] <- paste0(" * Incorrect vector length for: `", names(num)[i],
                          "`; expected length ", numlen[i], ", got ",
                          length(num[[i]]), "\n")
        j <- j + 1
      }
    }

  }

  if (!missing(logi)) {

    if (missing(logilen)) logilen <- rep(1, length(logi))

    logi.check1 <- vapply(logi, is.logical, logical(1))
    logi.check2 <- mapply(function(x, y) {
                            if (y == 0) TRUE else length(x) == y
                     }, logi, logilen)

    if (any(!logi.check1)) {
      j <- length(fails) + 1
      for (i in which(!logi.check1)) {
        fails[j] <- paste0(" * Incorrect type for: `", names(logi)[i],
                          "`; expected `logical`, got `",
                          mode(logi[[i]]), "`\n")
        j <- j + 1
      }
    }

    if (any(!logi.check2)) {
      j <- length(fails) + 1
      for (i in which(!logi.check2)) {
        fails[j] <- paste0(" * Incorrect vector length for: `", names(logi)[i],
                          "`; expected length ", logilen[i], ", got ",
                          length(logi[[i]]), "\n")
        j <- j + 1
      }
    }

  }

  if (length(fails) > 0) {
    fails <- c("\n", fails)
    stop(fails, call. = FALSE) 
  } else return()

}
