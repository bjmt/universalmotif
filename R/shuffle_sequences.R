#' Shuffle input sequences.
#'
#' @param sequences XStringSet objects.
#' @param inter Logical. Shuffle letters between sequences? If \code{FALSE},
#'              then the alphabet frequencies for each sequence are maintained.
#' @param keepDI Logical. Keep dinucleotide frequencies. Only supported for
#'               DNA and RNA.
#'
#' @return XStringSet object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
shuffle_sequences <- function(sequences, inter = FALSE, keepDI = FALSE) {

  alph <- sequences@elementType
  if (alph == "DNAString") {
    alphabet <- "DNA"
    alph.split <- DNA_BASES
  } else if (alph == "RNAString") {
    alphabet <- "RNA"
    alph.split <- RNA_BASES
  } else if (alph == "AAString") {
    alphabet <- "AA"
  } else {
    alphabet <- "B"
  }

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
  seqs <- as.character(sequences)
  widths <- width(sequences)
  seqs <- lapply(seqs, function(x) strsplit(x, "")[[1]])

  if (!keepDI) {

    if (!inter) {
      seqs <- mapply(function(x, y) sample(x, y), seqs, widths,
                    SIMPLIFY = FALSE)
      seqs <- lapply(seqs, function(x) paste(x, collapse = ""))
    } else if (inter) {
      seqs.all <- unlist(seqs)
      seqs.all <- sample(seqs.all, length(seqs.all))
      seqs.i <- mapply(function(x, y) as.factor(rep(x, y)),
                      seq_len(length(sequences)), widths, SIMPLIFY = FALSE)
      seqs.i <- unlist(seqs.i)
      seqs <- split(seqs.all, seqs.i)
      seqs <- lapply(seqs, function(x) paste(x, collapse = ""))
    } else stop("'inter' must be TRUE or FALSE")

  } else if (keepDI) {

    if (!inter) {
      monofreqs <- table(unlist(seqs))
      total.letters <- sum(monofreqs)
      monofreqs <- vapply(monofreqs, function(x) x / total.letters,
                          numeric(1))
      ditrans <- oligonucleotideTransitions(sequences, as.prob = TRUE)
      seqs.out <- vector("list", length(sequences))
      for (i in seq_len(length(sequences))) {
        seqs.out[[i]] <- rep(NA, widths[i])
        seqs.out[[i]][1] <- sample(alph.split, 1, prob = monofreqs)
        for (j in 2:widths[i]) {
          previous.nuc <- seqs.out[[i]][j - 1]
          curr.prob <- ditrans[previous.nuc, ]
          seqs.out[[i]][j] <- sample(alph.split, 1, prob = curr.prob)
        }
      }
      seqs <- lapply(seqs.out, function(x) paste(x, collapse = ""))
    } else if (inter) {
      seqs.out <- vector("list", length(sequences))
      for (i in seq_len(length(sequences))) {
        monofreqs <- table(seqs[[i]])
        monofreqs <- monofreqs / widths[i]
        ditrans <- oligonucleotideTransitions(sequences[i], as.prob = TRUE)
        seq.out <- rep(NA, widths[i])
        seq.out[1] <- sample(alph.split, 1, prob = monofreqs)
        for (j in 2:widths[i]) {
          previous.nuc <- seq.out[j - 1]
          curr.prob <- ditrans[previous.nuc, ]
          seq.out[j] <- sample(alph.split, 1, prob = curr.prob)
        }
        seqs.out[[i]] <- seq.out
      }
      seqs <- lapply(seqs.out, function(x) paste(x, collapse = ""))
    } else stop("'inter' must be TRUE or FALSE")

  } 

  seqs <- unlist(seqs)

  if (alphabet == "DNA") {
    seqs <- DNAStringSet(seqs)
  } else if (alphabet == "RNA") {
    seqs <- RNAStringSet(seqs)
  } else if (alphabet == "AA") {
    seqs <- AAStringSet(seqs)
  } else seqs <- BStringSet(seqs)

  names(seqs) <- seq.names

  seqs

}
