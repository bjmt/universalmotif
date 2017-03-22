######################################################################
## Benjamin Tremblay
##
## Assorted shared internal functions.
##
######################################################################

check_filter_args <- function(x) {
  for (i in which(is.na(x))) stop("'", names(x[i]),
                                  "' must be NULL or numeric", call. = FALSE)
  for (i in which(!vapply(x, is.null, logical(1)))) {
    if (!is.numeric(x[[i]]) || length(x[[i]]) != 1) {
      stop("'", names(x[i]), "' must be NULL or numeric of length 1",
           call. = FALSE)
    }
  }
}

check_logi_args <- function(x) {
  for (i in which(is.na(x))) stop("'", names(x[i]), "' must be logical",
                                  call. = FALSE)
  for (i in which(!vapply(x, is.logical, logical(1)))) {
    stop("'", names(x[i]), "' must be logical", call. = FALSE)
  }
  for (i in which(!vapply(x, function(x) length(x) == 1, logical(1)))) {
    stop("'", names(x[i]), "' must be a logical of length 1", call. = FALSE)
  }
}

check_out_class <- function(x) {
  if (!x %in% c("matrix-1", "matrix-2")) {
    stop("\"", x, "\" is not a supported 'out_class'", call. = FALSE)
  }
}

position_icscore <- function(motif, bkg = c(0.25, 0.25, 0.25, 0.25), type,
                             pseudoweight = 0.8, nsites = 100) {

  if (type == "PCM") {
      possum <- sum(motif)
      motif <- pcm_to_ppm(motif, possum, pseudoweight)
  }
  if (type == "PWM") motif <- pwm_to_ppm(motif, background = bkg)
  if (type == "ICM") return(sum(motif))

  if (type == "PPM") {
    motif <- ppm_to_pcm(motif, nsites = nsites)
    motif <- pcm_to_ppm(motif, pseudoweight = pseudoweight)
  }

  ic1 <- motif[1] * log2(motif[1] / bkg[1])
  ic2 <- motif[2] * log2(motif[2] / bkg[2])
  ic3 <- motif[3] * log2(motif[3] / bkg[3])
  ic4 <- motif[4] * log2(motif[4] / bkg[4])

  ic <- sum(ic1, ic2, ic3, ic4)

  return(ic)

}

ppm_to_icm <- function(position, bkg = c(0.25, 0.25, 0.25, 0.25),
                       IC_floor = FALSE) {
  # NOTE: different background frequencies can result in IC higher than 2!
  # not sure how to deal with this..
  if (is.null(bkg)) {
    bkg <- rep(1 / length(position), length(position))
  }
  for (i in seq_along(position)) {
    position[i] <- position[i] * log2(position[i] / bkg[i])
    if (IC_floor && position[i] < 0) position[i] <- 0
  }
  return(position)
}

# is this possible??
# icm_to_pwm <- function(position, bkg = c(0.25, 0.25, 0.25, 0.25)) {
#   for (i in seq_along(position)) {
#     position[i] <- position[i] /
#   }
# }

pcm_to_ppm <- function(position, possum, pseudoweight = 0.8) {
  if (missing(possum)) possum <- sum(position)
  pos <- vapply(position, function(x)
                (x + (pseudoweight / 4)) / (possum + pseudoweight),
                double(1))
  return(pos)
}

ppm_to_pcm <- function(position, nsites = 100) {
  if (length(nsites) == 0) nsites <- 100
  pos <- vapply(position, function(x) round(x * nsites), numeric(1))
  if (sum(pos) != nsites) {
    fix <- nsites - sum(pos)
    pos[which(range(pos)[2] == pos)] <- pos[which(range(pos)[2] == pos)] + fix
  }
  return(pos)
}

ppm_to_pwm <- function(position, background = c(0.25, 0.25, 0.25, 0.25),
                       pseudoweight = 0.8, nsites = 100) {
  if (length(nsites) == 0) nsites <- 100
  position <- ppm_to_pcm(position, nsites = nsites)
  position <- pcm_to_ppm(position, pseudoweight = pseudoweight)
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
    warning("get_consensus cannot handle ICM type for now", call. = FALSE)
    return(NULL)
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
  if (letter == "A") return(c(1, 0, 0, 0))
  if (letter == "C") return(c(0, 1, 0, 0))
  if (letter == "G") return(c(0, 0, 1, 0))
  if (letter == "T" || letter == "U") return(c(0, 0, 0, 1))
  if (letter == "R") return(c(0.5, 0, 0.5, 0))
  if (letter == "Y") return(c(0, 0.5, 0, 0.5))
  if (letter == "M") return(c(0.5, 0.5, 0, 0))
  if (letter == "K") return(c(0, 0, 0.5, 0.5))
  if (letter == "S") return(c(0, 0.5, 0.5, 0))
  if (letter == "W") return(c(0.5, 0, 0, 0.5))
  if (letter == "H") return(c(0.333, 0.333, 0.001, 0.333))
  if (letter == "B") return(c(0.001, 0.333, 0.333, 0.333))
  if (letter == "V") return(c(0.333, 0.333, 0.333, 0.001))
  if (letter == "D") return(c(0.333, 0.001, 0.333, 0.333))
  if (letter == "N") return(c(0.25, 0.25, 0.25, 0.25))
  if (letter == "+") return(c(0.25, 0.25, 0.25, 0.25))
  if (letter == "-") return(c(0.25, 0.25, 0.25, 0.25))
  if (letter == ".") return(c(0.25, 0.25, 0.25, 0.25))
  stop("not an IUPAC symbol")
}

withinlistvapply <- function(X, FUN, FUN.VALUE, ...) {
  if (!is.list(X)) stop("only works across lists")
  thelength <- unique(vapply(X, function(x) length(x), integer(1)))
  if (length(thelength) != 1) stop("list entries are not of the same length")
  tmplist <- as.list(seq_len(thelength)) 
  for (i in seq_len(thelength)) {
    for (j in seq_along(X)) {
      tmplist[[i]][j] <- X[[j]][i]
    }
  }
  finalvector <- vapply(X = tmplist, FUN = FUN, FUN.VALUE = FUN.VALUE, ...)
  return(finalvector)
}

string_to_pfm <- function(thestring, alphabet) {
  alphtable <- table(strsplit(thestring, split = "")[[1]])
  if (alphabet == "DNA") {
    pospfm <- c("A" = alphtable["A"], "C" = alphtable["C"],
                "G" = alphtable["G"], "T" = alphtable["T"])
  }
  if (alphabet == "RNA") {
    pospfm <- c("A" = alphtable["A"], "C" = alphtable["C"],
                "G" = alphtable["G"], "U" = alphtable["U"])
  }
  pospfm[is.na(pospfm)] <- 0
  return(pospfm)
}
