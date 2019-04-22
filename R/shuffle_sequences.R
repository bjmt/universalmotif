#' Shuffle input sequences.
#'
#' Given a set of input sequences, shuffle the letters within those
#' sequences with any k-let size.
#'
#' @param sequences \code{\link{XStringSet}} Set of sequences to shuffle. Works
#'    with any set of characters.
#' @param k `numeric(1)` K-let size.
#' @param method `character(1)` One of `c('euler', 'markov', 'linear', 'random')`.
#'    Only relevant is `k > 1`. See details. The `'random'` method is deprecated
#'    and will be removed in the next minor version.
#' @param leftovers `character(1)` For `method = 'random'`. One of
#'    `c('asis', 'first', 'split', 'discard')`.
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [shuffle_sequences()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for large jobs (such as
#'    shuffling billions of letters). Furthermore, the behaviour of `progress = TRUE`
#'    is changed if `BP = TRUE`; the default \pkg{BiocParallel} progress bar will
#'    be shown (which unfortunately is much less informative).
#'
#' @return \code{\link{XStringSet}} The input sequences will be returned with 
#'    identical names and lengths.
#'
#' @details
#'    If `method = 'markov'`, then the Markov model is used to
#'    generate sequences which will maintain (on average) the k-let
#'    frequencies. Please note that this method is not a 'true' shuffling, and
#'    for short sequences (e.g. <100bp) this can result in slightly more
#'    dissimilar sequences versus true shuffling. See
#'    \insertCite{markovmodel;textual}{universalmotif} for a discussion on the
#'    topic.
#'
#'    If `method = 'euler'`, then the sequence shuffling method proposed by
#'    \insertCite{markovmodel2;textual}{universalmotif} is used. As opposed
#'    to the 'markov' method, this one preserves exact k-let frequencies. This
#'    is done by creating a k-let edge graph, then following a
#'    random Eulerian walk through the graph. Not all walks will use up all
#'    available letters however, so the cycle-popping algorithm proposed by
#'    \insertCite{eulerAlgo;textual}{universalmotif} is used to find a
#'    random Eulerian path. A side effect of using this method is that the
#'    starting and ending sequence letters will remain unshuffled.
#'
#'    If `method = 'linear'`, then the input sequences are split linearly
#'    every `k` letters; for example, for `k = 3` 'ACAGATAGACCC' becomes
#'    'ACA GAT AGA CCC'; after which these `3`-lets are shuffled randomly.
#'
#'    Do note however, that the `method` parameter is only relevant for `k > 1`.
#'    For `k = 1`, a simple `sample` call is performed.
#'
#' @references
#'    \insertRef{markovmodel2}{universalmotif}
#'
#'    \insertRef{markovmodel}{universalmotif}
#'
#'    \insertRef{eulerAlgo}{universalmotif}
#'
#' @examples
#' sequences <- create_sequences()
#' sequences.shuffled <- shuffle_sequences(sequences, k = 2)
#'
#' @seealso [create_sequences()], [scan_sequences()], [enrich_motifs()],
#'    [shuffle_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
shuffle_sequences <- function(sequences, k = 1, method = "euler",
                               leftovers = "asis", progress = FALSE,
                               BP = FALSE) {

  # Idea: Moving-window markov shuffling. Get k-let frequencies in windows,
  #       and generate new letters based on local probabilities.

  # Idea: For k>1, do like k=1 and call utf8ToInt(sequence) then work with
  #       ints. This could make looking at the previous k-let in the
  #       markov and euler methods cheaper, instead of having to paste
  #       the previous strings. Maybe something like:
  #       as.integer(factor(utf8ToInt(sequence)))
  #       Then follow the method used in scan_sequences() to index k-lets.
  #
  #       utf8ToInt() is about 10X faster than strsplit() and safeExplode().

  # TODO: I don't think the k>1 methods are safe to use with sequences
  #       containing non-standard letters (e.g. 'N' in DNA sequences).

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!method %in% c("markov", "linear", "random", "euler")) {
    method_check <- paste0(" * Incorrect 'shuffle.method': expected `euler`, `markov`, ",
                           "`linear` or `random`; got `",
                           method, "`")
    method_check <- wmsg2(method_check)
    all_checks <- c(all_checks, method_check)
  }
  if (!leftovers %in% c("asis", "first", "split", "discard")) {
    leftovers_check <- paste0(" * Incorrect 'shuffle.leftovers': expected ",
                              "`asis`, `first`, `split` or `discard`; got `",
                              leftovers, "`")
    leftovers_check <- wmsg2(leftovers_check)
    all_checks <- c(all_checks, leftovers_check)
  }
  char_check <- check_fun_params(list(method = args$method,
                                      leftovers = args$leftovers),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(k = args$k), 1, FALSE, TYPE_NUM)
  s4_check <- check_fun_params(list(sequences = args$sequences),
                               numeric(), logical(), TYPE_S4)
  logi_check <- check_fun_params(list(progress = args$progress, BP = args$BP),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(all_checks, char_check, num_check, s4_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  alph <- seqtype(sequences)

  seq.names <- names(sequences)

  if (k == 1) {
    sequences <- lapply_(as.character(sequences), shuffle_k1,
                         PB = progress, BP = BP)
  } else {
    switch(method,
      "euler" = {
      sequences <- lapply_(as.character(sequences), shuffle_euler, k = k,
                           PB = progress, BP = BP)
      },
      "markov" = {
        if (seqtype(sequences) %in% c("DNA", "RNA"))
          sequences <- lapply_(sequences, shuffle_markov, k = k,
                               PB = progress, BP = BP)
        else
          sequences <- lapply_(as.character(sequences), shuffle_markov_any,
                               k = k, PB = progress, BP = BP)
      },
      "random" = {
        warning("The 'random' method is deprecated, please use 'linear' or 'markov'",
                immediate. = TRUE)
        sequences <- lapply_(as.character(sequences), shuffle_random, k = k,
                             leftover = leftovers, PB = progress, BP = BP)
      },
      "linear" = {
        sequences <- lapply_(as.character(sequences), shuffle_linear, k = k,
                             PB = progress, BP = BP)
      },
      stop("incorrect 'k' and 'method' combination")
    )
  }

  sequences <- unlist(sequences)

  sequences <- switch(alph, "DNA" = DNAStringSet(sequences),
                      "RNA" = RNAStringSet(sequences),
                      "AA" = AAStringSet(sequences),
                      BStringSet(sequences))

  if (!is.null(seq.names)) names(sequences) <- seq.names

  sequences

}

shuffle_k1 <- function(sequence) {
  intToUtf8(sample(utf8ToInt(sequence)))
}

# this creates truly random k-lets; unfortunately this has the side effect of
# leaving behind 'leftover'-lets smaller than k
shuffle_random <- function(sequence, k, leftover.strat = "asis", mode = 1) {
  # - benchmark timings: >200 times slower than shuffle_k1/shuffle_linear
  # - runtime decreases with increasing k!
  # -  _very_  memory inefficient!!! lots of gc

  if (mode == 1) {
    seq.len <- nchar(sequence)
    seqs1 <- safeExplode(sequence)
  } else if (mode == 2) {
    seq.len <- length(sequence)
    seqs1 <- sequence
  }

  seqs.k <- lapply(seq_len(k),
                   function(k) {
                     k <- k - 1
                     if (k > 0) {
                       seqs.k <- seqs1[-seq_len(k)]
                       seqs.k <- c(seqs.k, character(k))
                       matrix(seqs.k)
                     } else matrix(seqs1)
                   })
  seqs.k <- do.call(cbind, seqs.k)
  seqs.k.n <- nrow(seqs.k)
  seqs.k <- seqs.k[-c((seqs.k.n - k + 2):seqs.k.n), ]

  seqs.k.n <- nrow(seqs.k)
  seqs.k.n.len <- seq_len(seqs.k.n)
  new.seq <- matrix(character(round((seqs.k.n + 0.1) / 2) * k), nrow = k)

  seqs.k.new.i <- sample(seqs.k.n.len)

  # major time + mem sink:
  new.seq <- shuffle_random_loop(seqs.k.n - 1, k, seqs.k.new.i - 1, new.seq,
                                 seqs.k)

  seqs2 <- !seqs1 %in% new.seq

  new.seq <- as.character(new.seq)

  leftover <- seqs1[seqs2]

  if (length(leftover) > 0)
    switch(leftover.strat,
      "last" = {
        new.seq <- c(new.seq, leftover)
      },
      "asis" = {
        new.seq <- new.seq[new.seq != ""]
        leftover <- seqs1
        leftover[!seqs2] <- ""
        toadd <- length(leftover) - length(new.seq)
        toadd.left <- round((toadd + 0.1) / 2)
        toadd.right <- toadd - toadd.left
        toadd.left <- character(toadd.left)
        toadd.right <- character(toadd.right)
        new.seq <- c(toadd.left, new.seq, toadd.right)
        new.seq <- matrix(c(leftover, new.seq), nrow = 2, byrow = TRUE)
        new.seq <- as.character(new.seq)
      },
      "first" = {
        new.seq <- c(leftover, new.seq)
      },
      "split" = {
        if (length(leftover) == 1) {
          new.seq <- c(leftover, new.seq)
        } else {
          left.len <- length(leftover)
          left.left <- round((left.len + 0.1) / 2)
          left.right <- left.len - left.left
          new.seq <- c(leftover[seq_len(left.left)], new.seq,
                       leftover[left.right:length(leftover)])
        }
      },
      "discard" = {
        # do nothing
      },
      stop("unknown 'leftovers' arg")
    )

  if (mode == 1) new.seq <- collapse_cpp(new.seq)
  else new.seq <- new.seq[new.seq != ""]

  new.seq

}

# this function won't leave any 'leftover' letters behind, but the k-lets are
# predictable; the sequence is split up linearly every 'k'-letters
shuffle_linear <- function(sequence, k, mode = 1) {
  # - benchmark timings: about 1.25 times slower than shuffle_k1
  # - runtime is independent of k!

  # somehow the safeExplode() + collapse_cpp() combo is faster than
  # utf8ToInt() + intToUtf8()!

  if (mode == 1) {
    seq.len <- nchar(sequence)
    seq1 <- safeExplode(sequence)
  } else if (mode == 2) {
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
    seq2 <- seq1[-c(seq_len(left.keep), right.keep:seq.len)]
    left.keep <- seq1[seq_len(left.keep)]
    right.keep <- seq1[right.keep:length(seq1)]
  }

  seq.split <- matrix(seq2, nrow = k)

  seq.split.ncol <- ncol(seq.split)
  new.i <- sample.int(seq.split.ncol, seq.split.ncol)
  seq.split <- seq.split[, new.i]

  seq.split <- as.character(seq.split)

  # TODO: try and insert left.keep and right.keep somewhere inside seq.split?
  new.seq <- c(left.keep, seq.split, right.keep)
  if (mode == 1) new.seq <- collapse_cpp(new.seq)

  new.seq

}

shuffle_markov <- function(sequence, k) {
  # - benchmark timings: >10 times slower than shuffle_k1/shuffle_linear
  # - runtime is indenpendent of k!
  # - appears to be more memory efficient than even shuffle_k1...

  sequence <- DNAStringSet(sequence)
  seq.width <- width(sequence)

  # the Biostrings functions here sometimes cause C stack related crashes...
  freqs <- oligonucleotideFrequency(sequence, width = k, as.prob = TRUE)
  freqs <- colMeans(freqs)

  # using both oligonucleotideFrequency + oligonucleotideTransitions could be
  # avoided here, see letter_freqs; regardless, the extra cost is quite minor
  trans <- oligonucleotideTransitions(sequence, k - 1, 1, as.prob = TRUE)
  trans <- t(trans)
  trans[is.na(trans)] <- 1 / ncol(trans) / 1000  # if set to zero, sometimes
                                                 # get all probs = 0
  seqout <- character(seq.width)
  first.k <- sample(names(freqs), 1, prob = freqs)
  first.k <- safeExplode(first.k)
  seqout[seq_len(k)] <- first.k

  trans.cols <- colnames(trans)
  names(trans.cols) <- trans.cols
  # main time + mem sink:
  seqout <- shuffle_markov_loop(k, seq.width, k, seqout, DNA_BASES, trans,
                                trans.cols)

  seqout

}

shuffle_markov_any <- function(sequence, k) {
  # - almost as fast as shuffle_markov, but much less memory efficient
  # - runtime increases with increasing k

  seq1 <- safeExplode(sequence)
  lets.uniq <- sort(unique(seq1))
  seq.width <- length(seq1)

  let.info <- letter_freqs(seq1, k)
  freqs <- let.info$frequencies$freqs
  trans <- t(let.info$transitions)
  trans[trans == 0] <- 1 / nrow(trans) / 1000

  seqout <- character(seq.width)
  first.k <- sample(let.info$frequencies$lets, 1, prob = freqs)
  first.k <- safeExplode(first.k)
  seqout[seq_len(k)] <- first.k

  trans.cols <- colnames(trans)
  names(trans.cols) <- trans.cols
  seqout <- shuffle_markov_loop(k, seq.width, k, seqout, lets.uniq, trans,
                                trans.cols)

  seqout

}

shuffle_euler <- function(sequence, k) {
  # runtime increases with k
  # about 2X as slow and 2X as many mem allocs vs shuffle_markov_any()

  # no idea why, but for now this has to be done
  if (tolower(sequence) != sequence && toupper(sequence) != sequence)
    stop("lower and upper case letters cannot both be used for method = 'euler'")

  seq <- safeExplode(sequence)
  seqlen <- length(seq)
  alph <- sort(unique(seq))

  alph.i <- seq_along(alph)
  names(alph.i) <- alph

  first <- seq[seq_len(k - 1)]
  last <- collapse_cpp(seq[(seqlen - k + 2):seqlen])

  kletsm1 <- get_klets(alph, k - 1)
  # second slowest step
  klets <- letter_freqs(seq, k, "freqs", FALSE, alph)$counts

  edgematrix <- matrix(klets$counts, nrow = length(kletsm1), byrow = TRUE,
                       dimnames = list(kletsm1, alph))

  to.keep <- apply(edgematrix, 1, function(x) any(x > 0))
  kletsm1 <- kletsm1[to.keep]
  edgematrix <- edgematrix[to.keep, ]

  n <- length(kletsm1)

  lastlets <- get_lastlets(edgematrix, last, k, kletsm1, alph)

  for (i in seq_len(n)) {
    edgematrix[i, alph.i[lastlets[[i]]]] <- edgematrix[i, alph.i[lastlets[[i]]]] - 1
  }

  edgelist <- vector("list", n)
  names(edgelist) <- kletsm1
  for (i in seq_len(n)) {
    lets <- character()
    for (j in seq_along(alph)) {
      lets <- c(lets, rep(alph[j], edgematrix[i, j]))
    }
    edgelist[[i]] <- collapse_cpp(c(sample(lets), lastlets[[i]]))
  }
  edgelist <- unlist(edgelist)
  indices <- integer(n)
  names(indices) <- kletsm1

  # by far the slowest step
  out <- eulerian_walk_cpp(edgelist, first, seqlen, k, safeExplode(last), indices)

  collapse_cpp(out)

}

get_lastlets <- function(edgematrix, last, k, kletsm1, alph) {

  # Cycle-popping algorithm: RandomTreeWithRoot() (Propp and Wilson).
  # Only needs to be done with last letters of each vertex to obtain a
  # successfull Eulerian walk for entire sequence (Altschul and Erickson).

  n <- length(kletsm1)
  lastlets <- character(n)
  names(lastlets) <- kletsm1

  vertices <- logical(n)
  names(vertices) <- kletsm1

  vertices[last] <- TRUE  # the last (k-1)-let in the sequence is the tree root

  for (i in seq_len(n)) {
    u <- i
    while (!vertices[u]) {
      lastlets[u] <- sample(alph, 1, prob = edgematrix[u, ])
      u <- collapse_cpp(c(safeExplode(names(lastlets[u]))[-1], lastlets[u]))
    }
    u <- i
    while (!vertices[u]) {
      vertices[u] <- TRUE
      u <- collapse_cpp(c(safeExplode(names(lastlets[u]))[-1], lastlets[u]))
    }
  }

  as.list(lastlets)

}

#' Get k-let frequencies.
#'
#' @param seqs1 <CHAR> Split sequence, usually from strsplit(seq, "")[[1]].
#' @param k <INT> k-let size.
#' @param to.return <CHAR> Return k-let frequencies and/or transition matrix.
#' @param as.prob <BOOL> Return k-let counts or probabilities.
#' @param alph <CHAR> Alphabet letters, split.
#'
#' @return List, with entries 'transitions' (matrix), 'frequencies' (data.frame),
#'    and/or 'counts' (data.frame).
#'
#' @noRd
letter_freqs <- function(seqs1, k, to.return = c("freqs", "trans"),
                         as.prob = TRUE, alph = NULL) {
  # ~3 times slower than Biostrings::oligonucleotideTransitions

  if (is.null(alph)) lets.uniq <- sort(unique(seqs1))
  else lets.uniq <- sort(alph)

  possible.lets <- get_klets(lets.uniq, k)
  possible.lets <- data.frame(lets = possible.lets, stringsAsFactors = FALSE)

  seqs.let <- single_to_k(seqs1, k)

  seqs.counts <- table_cpp(seqs.let)
  seqs.counts <- data.frame(lets = names(seqs.counts), counts = as.numeric(seqs.counts),
                            stringsAsFactors = FALSE)
  final.table <- merge(possible.lets, seqs.counts, by = "lets", all.x = TRUE)
  final.table$counts[is.na(final.table$counts)] <- 0

  total.counts <- sum(final.table$counts)
  final.table$freqs <- final.table$counts / total.counts

  out <- list()

  if ("freqs" %in% to.return) {
    if (as.prob) out$frequencies <- final.table[, -2]
    else out$counts <- final.table[, -3]
  }

  if ("trans" %in% to.return) {
    trans <- t(matrix(final.table$counts, nrow = length(lets.uniq)))
    letskm1 <- get_klets(lets.uniq, k - 1)
    colnames(trans) <- lets.uniq
    rownames(trans) <- letskm1
    trans.sums <- rowSums(trans)
    if (as.prob)
      for (i in seq_len(nrow(trans))) trans[i, ] <- trans[i, ] / trans.sums[i]
    trans[is.nan(trans)] <- 0
    out$transitions <- trans
  }

  out

}
