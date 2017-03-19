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

# motif_icscore <- function(motif, bkg) {
#
#   # make use of TFBSTools::toICM?
#
# }
#
# TFBSTools:
    # toPWM: can generate a prob matrix or log-prob matrix
    # toICM: generates IC scores for each position from 0 to 2

get_consensus <- function(mot_matrix, alphabet = "DNA", type = "PPM") {

  pos <- mot_matrix

  if (type == "PCM") {
    possum <- sum(pos)
    pos <- vapply(pos, function(x) x / possum, double(1))
  }

  if (alphabet == "DNA") names(pos) <- c("A", "C", "G", "T")
  if (alphabet == "RNA") names(pos) <- c("A", "C", "G", "U")

  # single letter consensus:

  if (pos[1] > 0.5 && pos[1] > sort(pos)[3] * 2) return("A")
  if (pos[2] > 0.5 && pos[2] > sort(pos)[3] * 2) return("C")
  if (pos[3] > 0.5 && pos[3] > sort(pos)[3] * 2) return("G")
  if (pos[4] > 0.5 && pos[4] > sort(pos)[3] * 2) ifelse(alphabet = "DNA",
                                                        return("T"),
                                                        return("U"))

  # two letter consensus:

  if (pos[1] > 0.5) {
    if (names(sort(pos)[3]) == "C") return("M")
    if (names(sort(pos)[3]) == "G") return("R")
    if (names(sort(pos)[3]) %in% c("T", "U")) return("W")
  }

  if (pos[2] > 0.5) {
    if (names(sort(pos)[3]) == "A") return("M")
    if (names(sort(pos)[3]) == "G") return("S")
    if (names(sort(pos)[3]) %in% c("T", "U")) return("Y")
  }

  if (pos[3] > 0.5) {
    if (names(sort(pos)[3]) == "A") return("R")
    if (names(sort(pos)[3]) == "C") return("S")
    if (names(sort(pos)[3]) %in% c("T", "U")) return("K")
  }

  if (pos[4] > 0.5) {
    if (names(sort(pos)[3]) == "A") return("W")
    if (names(sort(pos)[3]) == "C") return("Y")
    if (names(sort(pos)[3]) == "G") return("K")
  }

  # three letter consensus:

  if (all(pos[c(1, 2, 4)] > 0.25)) return("H")
  if (all(pos[c(2, 3, 4)] > 0.25)) return("B")
  if (all(pos[c(1, 2, 3)] > 0.25)) return("V")
  if (all(pos[c(1, 3, 4)] > 0.25)) return("D")

  # no consensus: 

  return("N")

}
