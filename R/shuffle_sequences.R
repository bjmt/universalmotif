#' Shuffle input sequences.
#'
#' @param sequences XStringSet objects.
#' @param inter Logical. Shuffle letters between sequences? If \code{FALSE},
#'              then the alphabet frequencies for each sequence are maintained.
#' @param keepDI Logical. Keep dinucleotide frequencies. Only supported for
#'               DNA and RNA.
#' @param keepTRI Logical. Keep trinucleotide frequencies. Only supported for
#'                DNA and RNA.
#'
#' @return XStringSet object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
shuffle_sequences <- function(sequences, inter = FALSE, keepDI = FALSE,
                              keepTRI = FALSE) {

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

  if (keepDI && keepTRI) stop("only one of 'keepDI' and 'keepTRI' can be TRUE")

  if (!keepDI && !keepTRI) {

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
      ditrans <- oligonucleotideTransitions(sequences, as.prob = TRUE)
      difreq <- colSums(dinucleotideFrequency(sequences, as.prob = TRUE))
      seqs.out <- vector("list", length(sequences))
      for (i in seq_len(length(sequences))) {
        first.di <- sample(names(difreq), 1, prob = difreq)
        first.di <- strsplit(first.di, "")[[1]]
        seqs.out[[i]] <- rep(NA, widths[i])
        seqs.out[[i]][1] <- first.di[1]
        seqs.out[[i]][2] <- first.di[2]
        for (j in 3:widths[i]) {
          previous.nuc <- seqs.out[[i]][j - 1]
          curr.prob <- ditrans[previous.nuc, ]
          curr.prob[is.na(curr.prob)] <- 0.01
          seqs.out[[i]][j] <- sample(alph.split, 1, prob = curr.prob)
        }
      }
      seqs <- lapply(seqs.out, function(x) paste(x, collapse = ""))
    } else if (inter) {
      seqs.out <- vector("list", length(sequences))
      for (i in seq_len(length(sequences))) {
        difreq <- dinucleotideFrequency(sequences[i], as.prob = TRUE)[1, ]
        ditrans <- oligonucleotideTransitions(sequences[i], as.prob = TRUE)
        seq.out <- rep(NA, widths[i])
        first.di <- sample(names(difreq), 1, prob = difreq)
        first.di <- strsplit(first.di, "")[[1]]
        seq.out[1] <- first.di[1]
        seq.out[2] <- first.di[2]
        for (j in 3:widths[i]) {
          previous.nuc <- seq.out[j - 1]
          curr.prob <- ditrans[previous.nuc, ]
          curr.prob[is.na(curr.prob)] <- 0.01
          seq.out[j] <- sample(alph.split, 1, prob = curr.prob)
        }
        seqs.out[[i]] <- seq.out
      }
      seqs <- lapply(seqs.out, function(x) paste(x, collapse = ""))
    } else stop("'inter' must be TRUE or FALSE")

  } else if (keepTRI) {
    
    if (!inter) {
      ditrans <- oligonucleotideTransitions(sequences, as.prob = TRUE)
      tritrans <- oligonucleotideTransitions(sequences, 2, 1, as.prob = TRUE)
      trifreq <- colSums(trinucleotideFrequency(sequences, as.prob = TRUE))
      seqs.out <- vector("list", length(sequences))
      for (i in seq_len(length(sequences))) {
        seqs.out[[i]] <- rep(NA, widths[i])
        first.tri <- sample(names(trifreq), 1, prob = trifreq)
        first.tri <- strsplit(first.tri, "")[[1]]
        seqs.out[[i]][1] <- first.tri[1]
        seqs.out[[i]][2] <- first.tri[2]
        seqs.out[[i]][3] <- first.tri[3]
        for (j in 4:widths[i]) {
          previous.nuc1 <- seqs.out[[i]][j - 1]
          previous.nuc2 <- seqs.out[[i]][j - 2]
          previous.nuc <- paste0(previous.nuc2, previous.nuc2)
          curr.prob <- tritrans[previous.nuc, ]
          curr.prob[is.na(curr.prob)] <- 0
          seqs.out[[i]][j] <- sample(alph.split, 1, prob = curr.prob)
        }
      }
      seqs <- lapply(seqs.out, function(x) paste(x, collapse = ""))
    } else if (inter) {
      seqs.out <- vector("list", length(sequences))
      for (i in seq_len(length(sequences))) {
        ditrans <- oligonucleotideTransitions(sequences[i], as.prob = TRUE)
        trifreq <- trinucleotideFrequency(sequences[i], as.prob = TRUE)[1, ]
        tritrans <- oligonucleotideTransitions(sequences[i], 2, 1, as.prob = TRUE)
        seq.out <- rep(NA, widths[i])
        first.tri <- sample(names(trifreq), 1, prob = trifreq)
        first.tri <- strsplit(first.tri, "")[[1]]
        seq.out[1] <- first.tri[1]
        seq.out[2] <- first.tri[2]
        seq.out[3] <- first.tri[3]
        for (j in 4:widths[i]) {
          previous.nuc1 <- seq.out[j - 1]
          previous.nuc2 <- seq.out[j - 2]
          previous.nuc <- paste0(previous.nuc2, previous.nuc1)
          curr.prob <- tritrans[previous.nuc, ]
          curr.prob[is.na(curr.prob)] <- 0.01
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
