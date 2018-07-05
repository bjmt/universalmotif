#' Utility functions.
#'
#' @param position Numeric. A numeric vector representing the frequency or
#'    probability for each alphabet letter at a specific position.
#' @param bkg Numeric. Should be the same length as the alphabet length.
#' @param nsites Numeric. Number of sites motif originated from.
#' @param schneider_correction Logical. Apply sample size correction.
#' @param relative_entropy Logical. Calculate information content as
#'    relative entropy or Kullback-Leibler divergence.
#' @param pseudocount Numeric. Used to prevent zeroes in motif matrix.
#' @param smooth Logical. Apply pseudocount correction.
#' @param type Character. 'PCM', 'PPM', 'PWM' or 'ICM'.
#' @param alphabet Character. 'DNA' or 'RNA'.
#' @param letter Character. Any DNA, RNA, or AA IUPAC letter. Ambiguity letters
#'    are accepted.
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
#
# DNA_TETRA <- matrix(c(rep(DNA_TRI, 4), rep("A", 64), rep("C", 64),
                      # rep("G", 64), rep("T", 64)), ncol = 2)
# DNA_TETRA <- apply(DNA_TETRA, 1, paste, collapse = "")
# DNA_TETRA <- sort(DNA_TETRA)
#
# DNA_PENTA <- matrix(c(rep(DNA_TETRA, 4), rep("A", 256), rep("C", 256),
                      # rep("G", 256), rep("T", 256)), ncol = 2)
# DNA_PENTA <- apply(DNA_PENTA, 1, paste, collapse = "")
# DNA_PENTA <- sort(DNA_PENTA)
#
# DNA_HEXA <- matrix(c(rep(DNA_PENTA, 4,), rep("A", 1024), rep("C", 1024),
                     # rep("G", 1024), rep("T", 1024)), ncol = 2)
# DNA_HEXA <- apply(DNA_HEXA, 1, paste, collapse = "")
# DNA_HEXA <- sort(DNA_HEXA)

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

.internal_convert <- function(motifs, class = NULL, BPPARAM = bpparam()) {

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

.my_create_first <- function(bkg, sequences) {

  seq.width <- unique(width(sequences))
  if (seq.width < 2) return(matrix())

  emissions <- matrix(nrow = 16, ncol = seq.width - 1)
  rownames(emissions) <- DNA_DI
  colnames(emissions) <- seq_len(seq.width)[-seq.width]

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])
  seqs.split <- t(seqs.split)

  for (i in seq_len(seq.width - 1)) {
    current.seqs <- seqs.split[, i:(i + 1)]
    current.seqs <- apply(current.seqs, 1, paste, collapse = "")
    current.seqs <- DNAStringSet(current.seqs)
    emissions.i <- colSums(dinucleotideFrequency(current.seqs))
    emissions.i <- emissions.i / sum(emissions.i) #+ 0.00001
    emissions[, i] <- emissions.i
  }

  emissions

}

.my_create_second <- function(bkg, sequences) {

  seq.width <- unique(width(sequences))
  if (seq.width < 3) return(matrix())

  emissions <- matrix(nrow = 64, ncol = seq.width - 2)
  rownames(emissions) <- DNA_TRI
  colnames(emissions) <- seq_len(seq.width)[-c(seq.width - 1, seq.width)]

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])
  seqs.split <- t(seqs.split)

  for (i in seq_len(seq.width - 2)) {
    current.seqs <- seqs.split[, i:(i + 2)]
    current.seqs <- apply(current.seqs, 1, paste, collapse = "")
    current.seqs <- DNAStringSet(current.seqs)
    emissions.i <- colSums(trinucleotideFrequency(current.seqs))
    emissions.i <- emissions.i / sum(emissions.i)
    emissions[, i] <- emissions.i
  }

  emissions

}

.my_create_third <- function(bkg, sequences) {

  seq.width <- unique(width(sequences))
  if (seq.width < 4) return(matrix())

  emissions <- matrix(nrow = 256, ncol = seq.width - 3)
  rownames(emissions) <- DNA_TETRA
  colnames(emissions) <- seq_len(seq.width)[-c(seq.width - 2, seq.width - 1,
                                               seq.width)]

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])
  seqs.split <- t(seqs.split)

  for (i in seq_len(seq.width - 3)) {
    current.seqs <- seqs.split[, i:(i + 3)]
    current.seqs <- apply(current.seqs, 1, paste, collapse = "")
    current.seqs <- DNAStringSet(current.seqs)
    emissions.i <- colSums(oligonucleotideFrequency(current.seqs, 4, 1))
    emissions.i <- emissions.i / sum(emissions.i)
    emissions[, i] <- emissions.i
  }

  emissions

}

.my_create_fourth <- function(bkg, sequences) {

  seq.width <- unique(width(sequences))
  if (seq.width < 5) return(matrix())

  emissions <- matrix(nrow = 1024, ncol = seq.width - 4)
  rownames(emissions) <- DNA_PENTA
  colnames(emissions) <- seq_len(seq.width)[-c(seq.width - 3, seq.width - 2,
                                               seq.width - 1, seq.width)]

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])
  seqs.split <- t(seqs.split)

  for (i in seq_len(seq.width - 4)) {
    current.seqs <- seqs.split[, i:(i + 4)]
    current.seqs <- apply(current.seqs, 1, paste, collapse = "")
    current.seqs <- DNAStringSet(current.seqs)
    emissions.i <- colSums(oligonucleotideFrequency(current.seqs, 5, 1))
    emissions.i <- emissions.i / sum(emissions.i)
    emissions[, i] <- emissions.i
  }

  emissions

}

.my_create_fifth <- function(bkg, sequences) {

  seq.width <- unique(width(sequences))
  if (seq.width < 5) return(matrix())

  emissions <- matrix(nrow = 4096, ncol = seq.width - 5)
  rownames(emissions) <- DNA_HEXA
  colnames(emissions) <- seq_len(seq.width)[-c(seq.width - 4, seq.width - 3,
                                               seq.width - 2, seq.width - 1,
                                               seq.width)]

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])
  seqs.split <- t(seqs.split)

  for (i in seq_len(seq.width - 5)) {
    current.seqs <- seqs.split[, i:(i + 5)]
    current.seqs <- apply(current.seqs, 1, paste, collapse = "")
    current.seqs <- DNAStringSet(current.seqs)
    emissions.i <- colSums(oligonucleotideFrequency(current.seqs, 6, 1))
    emissions.i <- emissions.i / sum(emissions.i)
    emissions[, i] <- emissions.i
  }

  emissions

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

# for a motif of length 4, the transition matrix is:

#      pos0 pos1 pos2 pos3 pos4 
# pos0    0    1    0    0    0
# pos1    0    0    1    0    0
# pos2    0    0    0    1    0
# pos3    0    0    0    0    1
# pos4    1    0    0    0    0
