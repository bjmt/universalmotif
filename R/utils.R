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
                             pseudoweight = 0.8) {

  if (type == "PCM") {
      possum <- sum(motif)
      motif <- pcm_to_ppm(motif, possum, pseudoweight)
      type <- "PPM"
  }
  if (type == "PPM") {
    if (motif[1] == 0) ic1 <- 0 else ic1 <- motif[1] * log2(motif[1] / bkg[1])
    if (motif[2] == 0) ic2 <- 0 else ic2 <- motif[2] * log2(motif[2] / bkg[2])
    if (motif[3] == 0) ic3 <- 0 else ic3 <- motif[3] * log2(motif[3] / bkg[3])
    if (motif[4] == 0) ic4 <- 0 else ic4 <- motif[4] * log2(motif[4] / bkg[4])
  }

  ic <- sum(ic1, ic2, ic3, ic4)

  return(ic)

}
#
# TFBSTools:
    # toPWM: can generate a prob matrix or log-prob matrix
    # toICM: generates IC scores for each position from 0 to 2

pcm_to_ppm <- function(position, possum, pseudoweight) {
  pos <- vapply(position, function(x)
                (x + pseudoweight / 4) / (possum + pseudoweight),
                double(1))
}

get_consensus <- function(mot_matrix, alphabet = "DNA", type = "PPM",
                          pseudoweight = 0.8) {

  pos <- mot_matrix

  if (type == "PCM") {
    possum <- sum(pos)
    pos <- pcm_to_ppm(pos, possum, pseudoweight)
    type <- "PPM"
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
