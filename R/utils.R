ppm_to_icm <- function(position, bkg = c(0.25, 0.25, 0.25, 0.25),
                       IC_floor = FALSE, IC_ceiling = FALSE) {
  # NOTE: different background frequencies can result in IC higher than 2!
  # not sure how to deal with this..
  if (is.null(bkg) || missing(bkg)) {
    bkg <- rep(1 / length(position), length(position))
  }
  for (i in seq_along(position)) {
    position[i] <- position[i] * log2(position[i] / bkg[i])
    if (IC_floor && position[i] < 0) position[i] <- 0
    if (IC_ceiling && position[i] > 2) position[i] <- 2
  }
  return(position)
}

pcm_to_ppm <- function(position, possum, pseudoweight = 0.8) {
  if (missing(possum)) possum <- sum(position)
  pos <- vapply(position, function(x)
                (x + (pseudoweight / 4)) / (possum + pseudoweight),
                double(1))
  return(pos)
}

ppm_to_pcm <- function(position, nsites = 100) {
  if (length(nsites) == 0 || missing(nsites)) nsites <- 100
  pos <- vapply(position, function(x) round(x * nsites), numeric(1))
  if (sum(pos) != nsites) {
    fix <- nsites - sum(pos)
    pos[which(range(pos)[2] == pos)[1]] <- pos[which(range(pos)[2] == pos)[1]] + fix
  }
  return(pos)
}

ppm_to_pwm <- function(position, background = c(0.25, 0.25, 0.25, 0.25),
                       pseudoweight = 0.8, nsites = 100, smooth = TRUE) {
  if (length(nsites) == 0) nsites <- 100
  if (smooth) {
    position <- ppm_to_pcm(position, nsites = nsites)
    position <- pcm_to_ppm(position, pseudoweight = pseudoweight)
  }
  for (i in seq_along(position)) {
    position[i] <- log2(position[i] / background[i])
  }
  return(position)
}

pwm_to_ppm <- function(position, background = c(0.25, 0.25, 0.25, 0.25)) {
  for (i in seq_along(position)) {
    position[i] <- 2 ^ position[i]
    position[i] <- position[i] * background[i]
  }
  return(position)
}

position_icscore <- function(motif, bkg = c(0.25, 0.25, 0.25, 0.25), type,
                             pseudoweight = 0.8, nsites = 100) {

  if (length(nsites) == 0) nsites <- 100

  if (type == "PCM") {
      possum <- sum(motif)
      motif <- pcm_to_ppm(motif, possum, pseudoweight)
  }
  if (type == "PWM") motif <- pwm_to_ppm(motif, background = bkg)
  if (type == "ICM") return(sum(motif))

  if (type == "PPM") {
    motif <- ppm_to_pcm(position = motif, nsites = nsites)
    motif <- pcm_to_ppm(position = motif, pseudoweight = pseudoweight)
  }

  ic1 <- motif[1] * log2(motif[1] / bkg[1])
  ic2 <- motif[2] * log2(motif[2] / bkg[2])
  ic3 <- motif[3] * log2(motif[3] / bkg[3])
  ic4 <- motif[4] * log2(motif[4] / bkg[4])

  sum(ic1, ic2, ic3, ic4)

}

get_consensus <- function(mot_matrix, alphabet = "DNA", type = "PPM",
                          pseudoweight = 0.8) {

  pos <- mot_matrix

  if (type == "PCM") {
    possum <- sum(pos)
    pos <- pcm_to_ppm(pos, possum, pseudoweight)
    type <- "PPM"
  }

  if (type == "PWM") {
    pos <- pwm_to_ppm(pos)
    type <- "PPM"
  }

  if (type == "ICM") {
    stop("get_consensus cannot handle ICM type")
  }

  if (type == "PPM") {

    if (alphabet == "DNA") names(pos) <- c("A", "C", "G", "T")
    if (alphabet == "RNA") names(pos) <- c("A", "C", "G", "U")

    # single letter consensus:

    if (pos[1] > 0.5 && pos[1] > sort(pos)[3] * 2) return("A")
    if (pos[2] > 0.5 && pos[2] > sort(pos)[3] * 2) return("C")
    if (pos[3] > 0.5 && pos[3] > sort(pos)[3] * 2) return("G")
    if (pos[4] > 0.5 && pos[4] > sort(pos)[3] * 2) ifelse(alphabet == "DNA",
                                                          return("T"),
                                                          return("U"))

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

consensus_to_ppmAA <- function(letter) {
  if (letter %in% c("X", ".")) return(rep(0.05, 20))
  if (letter == "B") return(c(rep(0.001, 2), 0.491, rep(0.001, 8), 0.491,
                              rep(0.001, 8)))
  if (letter == "Z") return(c(rep(0.001, 3), 0.491, rep(0.001, 9), 0.491,
                              rep(0.001, 6)))
  AAs <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N",
           "P", "Q", "R", "S", "T", "V", "W", "Y")
  i <- which(AAs == letter)
  c(rep(0.001, i - 1), 0.981, rep(0.001, 20 - i))
}

get_consensusAA <- function(motif) {
  if (motif[3] >= 0.4 && motif[12] >= 0.4) return("B")
  if (motif[4] >= 0.4 && motif[14] >= 0.4) return("Z")
  if (all(motif == 0.05)) return("X")
  AAs <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P",
           "Q", "R", "S", "T", "V", "W", "Y")
  AAs[order(motif, decreasing = TRUE)[1]]
}

.internal_convert <- function(motifs, class = NULL) {

  if (is.null(class)) {
    if (class(motifs) == "list") {
      CLASS_PKG <- unique(vapply(motifs, function(x) attributes(class(x))$package,
                                 character(1)))
      CLASS <- unique(vapply(motifs, function(x) class(x), character(1)))
      CLASS_IN <- paste(CLASS_PKG, CLASS, sep = "-")
    } else {
      CLASS_PKG <- attributes(class(motifs))$package
      CLASS <- class(motifs)
      CLASS_IN <- paste(CLASS_PKG, CLASS, sep = "-")
    }
    return(CLASS_IN)
  } else {
    if (length(class) == 1 && class != "universalmotif-universalmotif") {
      tryCatch(motifs <- convert_motifs(motifs, class),
               error = function(e) message("motifs converted to class 'universalmotif'"))
    }
    return(motifs)
  }
  
}
