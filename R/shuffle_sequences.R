#' Shuffle input sequences.
#'
#' @param sequences XStringSet objects.
#' @param k Numeric. k-let size.
#' @param method Character. One of 'markov', 'linear', and 'random'.
#' @param leftovers Character. For \code{method = 'random'}.
#' @param BPPARAM See \code{\link[BiocParallel]{SerialParam}}.
#'
#' @return XStringSet object.
#'
#' @details
#'    If one of \code{keepDI} and \code{keepTRI} is \code{TRUE}, then
#'    the Markov model \insertCite{markovmodel}{universalmotif} is used to
#'    generate sequences which will maintain (on average) the di/trinucleotide
#'    frequencies. Please note that this method is not a 'true' shuffling, and
#'    for short sequences (e.g. <100bp) this can result in slightly more
#'    dissimilar sequences versus true shuffling. See
#'    \insertCite{markovmodel;textual}{universalmotif} and 
#'    \insertCite{markovmodel2;textual}{universalmotif} for a discussion on the
#'    topic.
#'
#' @references
#'    \insertRef{markovmodel}{universalmotif}
#'
#'    \insertRef{markovmodel2}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
shuffle_sequences <- function(sequences, k = 1, method = "linear",
                               leftovers = "asis", BPPARAM = SerialParam()) {

  alph <- class(sequences)

  seq.names <- names(sequences)

  if (k == 1) {
    sequences <- as.character(sequences)
    sequences <- bplapply(sequences, shuffle_k1, BPPARAM = BPPARAM)
  } else if (method == "markov") {
    sequences <- bplapply(sequences, shuffle_markov, k = k, BPPARAM = BPPARAM)
  } else if (method == "random") {
    sequences <- as.character(sequences)
    sequences <- bplapply(sequences, shuffle_random, k = k, leftover = leftovers,
                          BPPARAM = BPPARAM)
  } else if (method == "linear") {
    sequences <- as.character(sequences)
    sequences <- bplapply(sequences, shuffle_linear, k = k, BPPARAM = BPPARAM)
  } else stop("incorrect 'k' and 'method' combo")

  sequences <- unlist(sequences)

  if (alph == "DNAStringSet") {
    sequences <- DNAStringSet(sequences)
  } else if (alph == "RNAStringSet") {
    sequences <- RNAStringSet(sequences)
  } else if (alph == "AAStringSet") {
    sequences <- AAStringSet(sequences)
  } else sequences <- BStringSet(sequences)

  if (!is.null(seq.names)) names(sequences) <- seq.names

  sequences

}

# this creates truly random k-lets; unfortunately this has the side effect of
# leaving behind 'leftover'-lets smaller than k
shuffle_random <- function(sequence, k, leftover.strat) {

  seq.len <- nchar(sequence)

  seqs1 <- strsplit(sequence, "")[[1]]
  seqs2 <- rep(TRUE, seq.len)

  seqs.k <- lapply(seq_len(k),
                   function(k) {
                     k <- k - 1
                     if (k > 0) {
                       seqs.k <- seqs1[-c(1:k)]
                       seqs.k <- c(seqs.k, rep(NA, k))
                       matrix(seqs.k)
                     } else matrix(seqs1)
                   })

  seqs.k <- do.call(cbind, seqs.k)
  seqs.k <- seqs.k[-c((nrow(seqs.k) - k + 2):nrow(seqs.k)), ]

  new.seq <- matrix(nrow = k, ncol = round((nrow(seqs.k) + 0.1) / 2))
  for (i in 1:nrow(seqs.k)) {
    if (all(is.na(seqs.k))) break
    repeat {
      j <- sample(1:nrow(seqs.k), 1)  # this repeat loop can be removed
      let <- seqs.k[j, ]
      if (!any(is.na(let))) break
    }

    new.seq[, i] <- let

    del1 <- j - k + 1
    del2 <- j + k - 1
    if (del1 < 0) del1 <- 0
    if (del2 > nrow(seqs.k)) del2 <- nrow(seqs.k)

    seqs.k[del1:del2, ] <- NA
    seqs2[j:(j + k - 1)] <- FALSE

  }

  new.seq <- as.character(new.seq)
  new.seq <- new.seq[!is.na(new.seq)]

  leftover <- seqs1[seqs2]
  
  if (length(leftover) > 0 ){
    if (leftover.strat == "last") {
      new.seq <- c(new.seq, leftover)
    } else if (leftover.strat == "asis") {
      leftover <- seqs1
      leftover[!seqs2] <- NA
      toadd <- length(leftover) - length(new.seq)
      toadd.left <- round((toadd + 0.1) / 2)
      toadd.right <- toadd - toadd.left
      toadd.left <- rep(NA, toadd.left)
      toadd.right <- rep(NA, toadd.right)
      new.seq <- c(toadd.left, new.seq, toadd.right)
      new.seq <- matrix(c(leftover, new.seq), nrow = 2, byrow = TRUE)
      new.seq <- as.character(new.seq)
      new.seq <- new.seq[!is.na(new.seq)]
    } else if (leftover.strat == "first") {
      new.seq <- c(leftover, new.seq)
    } else if (leftover.strat == "split") {
      if (length(leftover) == 1) {
        new.seq <- c(leftover, new.seq)
      } else {
        left.len <- length(leftover)
        left.left <- round((left.len + 0.1) / 2)
        left.right <- left.len - left.left
        new.seq <- c(leftover[1:left.left], new.seq,
                     leftover[left.right:length(leftover)])
      }
    } else if (leftover.strat == "discard") {
      # do nothing
    } else stop("unknown 'leftovers' arg")
  }


  paste(new.seq, collapse = "")

}

# this function won't leave any 'leftover' letters behind, but the k-lets are
# predictable; the sequence is split up linearly every 'k'-letters
shuffle_linear <- function(sequence, k) {

  seq.len <- nchar(sequence)
  
  seq1 <- strsplit(sequence, "")[[1]]

  seq.mod <- seq.len %% k
  if (seq.mod == 1) {
    left.keep <- 1
    right.keep <- 0
    seq2 <- seq1[-left.keep]
    left.keep <- seq1[1]
    right.keep <- character(0)
  } else if (seq.mod == 0) {
    left.keep <- character(0)
    right.keep <- character(0)
    seq2 <- seq1
  } else {
    left.keep <- round((seq.mod + 0.1) / 2)
    right.keep <- seq.mod - left.keep
    seq2 <- seq1[-c(1:left.keep, right.keep:seq.len)]
    left.keep <- seq1[1:left.keep]
    right.keep <- seq1[right.keep:length(seq1)]
  }

  seq.len2 <- length(seq2)

  tosplit <- seq.len2 / k
  tosplit <- tosplit - 1

  tosplit.i <- vapply(0:tosplit, function(x) 1 + k * x, numeric(1))

  seq.split <- lapply(tosplit.i, function(x) matrix(seq2[x:(x + k -1)]))
  seq.split <- do.call(cbind, seq.split)

  new.i <- sample(seq_len(ncol(seq.split)), ncol(seq.split))
  seq.split <- seq.split[, new.i]

  seq.split <- as.character(seq.split)

  new.seq <- c(left.keep, seq.split, right.keep)
  new.seq <- paste(new.seq, collapse = "")

  new.seq

}

shuffle_markov <- function(sequence, k) {
  
  sequence <- DNAStringSet(sequence)
  freqs <- colSums(oligonucleotideFrequency(sequence, width = k,
                                            as.prob = TRUE))
  trans <- oligonucleotideTransitions(sequence, k - 1, 1, as.prob = TRUE)
  seqout <- rep(NA, width(sequence))
  first.k <- sample(names(freqs), 1, prob = freqs)
  first.k <- strsplit(first.k, "")[[1]]
  seqout[1:k] <- first.k
  for (i in (k + 1):width(sequence)) {
    previous.k <- seqout[(i - k + 1):(i - 1)]
    previous.k <- paste(previous.k, collapse = "")
    curr.prob <- trans[previous.k, ]
    curr.prob[is.na(curr.prob)] <- 0.0001
    seqout[i] <- sample(DNA_BASES, 1, prob = curr.prob)
  }
  seqout <- paste(seqout, collapse = "")
  seqout

}

shuffle_k1 <- function(sequence) {
  
  sequence <- as.character(sequence)
  sequence <- strsplit(sequence, "")[[1]]
  sequence <- sample(sequence, length(sequence))
  sequence <- paste(sequence, collapse = "")
  sequence

}
