#' Shuffle input sequences.
#'
#' Given a set of input sequences, shuffle the letters within those
#' sequences with any k-let size.
#'
#' @param sequences [Biostrings::XStringSet-class] For `method = 'markov'`,
#'    [Biostrings::DNAStringSet-class] and [Biostrings::RNAStringSet-class] only.
#' @param k `numeric(1)` K-let size.
#' @param method `character(1)` One of `c('markov', 'linear', 'random')`.
#'    Only relevant is `k > 1`. See details.
#' @param leftovers `character(1)` For `method = 'random'`. One of
#'    `c('asis', 'first', 'split', 'discard')`. See details.
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [shuffle_sequences()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for large jobs (such as
#'    shuffling billions of letters). Furthermore, the behaviour of `progress = TRUE`
#'    is changed if `BP = TRUE`; the default \pkg{BiocParallel} progress bar will
#'    be shown (which unfortunately is much less informative).
#'
#' @return [Biostrings::XStringSet-class] The input sequences will be returned with 
#'    identical names and lengths.
#'
#' @details
#'    If `method = 'markov'`, then the Markov model is used to
#'    generate sequences which will maintain (on average) the k-let
#'    frequencies. Please note that this method is not a 'true' shuffling, and
#'    for short sequences (e.g. <100bp) this can result in slightly more
#'    dissimilar sequences versus true shuffling. See
#'    \insertCite{markovmodel;textual}{universalmotif} and 
#'    \insertCite{markovmodel2;textual}{universalmotif} for a discussion on the
#'    topic.
#'
#'    If `method = 'linear'`, then the input sequences are split linearly
#'    every `k` letters; for example, for `k = 3` 'ACAGATAGACCC' becomes
#'    'ACA GAT AGA CCC'; afterwhich these `3`-lets are shuffled randomly. If
#'    `method = 'random'`, then `k`-lets are picked from the sequence
#'    completely randomly. This however can leave 'leftover' letters, where
#'    lone letter islands smaller than `k` are left. There are a few options
#'    provided to deal with these: `leftovers = 'asis'` will leave these
#'    letter islands in place; `leftovers = 'first'` will place these
#'    letters at the beginning of the sequence; `leftovers = 'split'`
#'    will place half of the leftovers at the beginning and end of the 
#'    sequence; `leftovers = 'discard'` simply gets rid of the leftovers.
#'
#'    Do note however, that the `method` parameter is only relevant for `k > 1`.
#'    For this, a simple `sample` call is performed.
#'
#' @references
#'    \insertRef{markovmodel2}{universalmotif}
#'
#'    \insertRef{markovmodel}{universalmotif}
#'
#' @examples
#' sequences <- create_sequences()
#' sequences.shuffled <- shuffle_sequences(sequences, k = 2)
#'
#' @seealso [create_sequences()], [scan_sequences()], [enrich_motifs()],
#'    [shuffle_motifs()]
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
shuffle_sequences <- function(sequences, k = 1, method = "linear",
                               leftovers = "asis", progress = FALSE,
                               BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(method = args$method,
                                      leftovers = args$leftovers),
                                 numeric(), logical(), "character")
  num_check <- check_fun_params(list(k = args$k), 1, FALSE, "numeric")
  s4_check <- check_fun_params(list(sequences = args$sequences),
                               numeric(), logical(), "S4")
  logi_check <- check_fun_params(list(progress = args$progress, BP = args$BP),
                                 numeric(), logical(), "logical")
  all_checks <- c(char_check, num_check, s4_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  alph <- class(sequences)

  seq.names <- names(sequences)

  if (k == 1) {
    sequences <- as.character(sequences)
    sequences <- lapply_(sequences, shuffle_k1, PB = progress, BP = BP)
  } else if (method == "markov") {
    sequences <- lapply_(sequences, shuffle_markov, k = k, PB = progress, BP = BP)
  } else if (method == "random") {
    sequences <- as.character(sequences)
    sequences <- lapply_(sequences, shuffle_random, k = k, leftover = leftovers,
                         PB = progress, BP = BP)
  } else if (method == "linear") {
    sequences <- as.character(sequences)
    sequences <- lapply_(sequences, shuffle_linear, k = k, PB = progress,
                         BP = BP)
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
shuffle_random <- function(sequence, k, leftover.strat, mode = 1) {

  if (mode == 1) {
    seq.len <- nchar(sequence)
    seqs1 <- strsplit(sequence, "")[[1]]
  } else if (mode == 2) {
    seq.len <- length(sequence)
    seqs1 <- sequence
  }

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

  if (mode == 1) new.seq <- collapse_cpp(new.seq)

  new.seq 

}

# this function won't leave any 'leftover' letters behind, but the k-lets are
# predictable; the sequence is split up linearly every 'k'-letters
shuffle_linear <- function(sequence, k, mode = 1) {

  if (mode == 1) {
    seq.len <- nchar(sequence)
    seq1 <- strsplit(sequence, "")[[1]]
  } else if (mode == 2){
    seq.len <- length(sequence)
    seq1 <- sequence
  }

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
    right.keep <- seq.len - right.keep + 1
    seq2 <- seq1[-c(1:left.keep, right.keep:seq.len)]
    left.keep <- seq1[1:left.keep]
    right.keep <- seq1[right.keep:length(seq1)]
  }

  seq.len2 <- length(seq2)

  tosplit <- seq.len2 / k
  tosplit <- tosplit - 1

  tosplit.i <- vapply(0:tosplit, function(x) 1 + k * x, numeric(1))

  seq.split <- lapply(tosplit.i, function(x) matrix(seq2[x:(x + k - 1)]))
  seq.split <- do.call(cbind, seq.split)

  new.i <- sample(seq_len(ncol(seq.split)), ncol(seq.split))
  seq.split <- seq.split[, new.i]

  seq.split <- as.character(seq.split)

  new.seq <- c(left.keep, seq.split, right.keep)
  if (mode == 1) {
    new.seq <- collapse_cpp(new.seq)
  }

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
  seqout <- collapse_cpp(seqout)
  seqout

}

shuffle_k1 <- function(sequence) {
  
  sequence <- as.character(sequence)
  sequence <- strsplit(sequence, "")[[1]]
  sequence <- sample(sequence, length(sequence))
  sequence <- collapse_cpp(sequence)
  sequence

}
