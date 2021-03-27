#' Shuffle input sequences.
#'
#' Given a set of input sequences, shuffle the letters within those
#' sequences with any k-let size.
#'
#' @param sequences \code{\link{XStringSet}} Set of sequences to shuffle. Works
#'    with any set of characters.
#' @param k `numeric(1)` K-let size.
#' @param method `character(1)` One of `c('euler', 'markov', 'linear')`.
#'    Only relevant if `k > 1`. See details.
#' @param nthreads `numeric(1)` Run [shuffle_sequences()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'    Note that no speed up will occur for jobs with only a single sequence.
#' @param rng.seed `numeric(1)` Set random number generator seed. Since shuffling
#'    can occur simultaneously in multiple threads using C++, it cannot communicate
#'    with the regular `R` random number generator state and thus requires an
#'    independent seed. Each individual sequence in an \code{\link{XStringSet}} object will be
#'    given the following seed: `rng.seed * index`. The default is to pick a random
#'    number as chosen by [sample()], which effectively is making [shuffle_sequences()]
#'    dependent on the R RNG state.
#' @param window `logical(1)` Shuffle sequences iteratively over windows instead
#'    of all at once.
#' @param window.size `numeric(1)` Window size. Can be a fraction less than one, or
#'    an integer representing the actual window size.
#' @param window.overlap `numeric(1)` Overlap between windows. Can be a fraction less
#'    than one, or an integer representing the actual overlap size.
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
#'    Fitch (1983) for a discussion on the
#'    topic.
#'
#'    If `method = 'euler'`, then the sequence shuffling method proposed by
#'    Altschul and Erickson (1985) is used. As opposed
#'    to the 'markov' method, this one preserves exact k-let frequencies. This
#'    is done by creating a k-let edge graph, then following a
#'    random Eulerian walk through the graph. Not all walks will use up all
#'    available letters however, so the cycle-popping algorithm proposed by
#'    Propp and Wilson (1998) is used to find a
#'    random Eulerian path. A side effect of using this method is that the
#'    starting and ending sequence letters will remain unshuffled.
#'
#'    If `method = 'linear'`, then the input sequences are split linearly
#'    every `k` letters. For example, for `k = 3` 'ACAGATAGACCC' becomes
#'    'ACA GAT AGA CCC'; after which these `3`-lets are shuffled randomly.
#'
#'    Do note however, that the `method` parameter is only relevant for `k > 1`.
#'    For `k = 1`, a simple shuffling is performed using the `shuffle` function
#'    from the C++ standard library.
#'
#' @references
#'
#' Altschul SF, Erickson BW (1985). “Significance of Nucleotide
#' Sequence Alignments: A Method for Random Sequence Permutation That
#' Preserves Dinucleotide and Codon Usage.” _Molecular Biology and
#' Evolution_, *2*, 526-538.
#' 
#' Fitch WM (1983). “Random sequences.” _Journal of Molecular
#' Biology_, *163*, 171-176.
#' 
#' Propp JG, Wilson DW (1998). “How to get a perfectly random sample
#' from a generic markov chain and generate a random spanning tree of
#' a directed graph.” _Journal of Algorithms_, *27*, 170-217.
#'
#' @examples
#' if (R.Version()$arch != "i386") {
#' sequences <- create_sequences()
#' sequences.shuffled <- shuffle_sequences(sequences, k = 2)
#' }
#'
#' @seealso [create_sequences()], [scan_sequences()], [enrich_motifs()],
#'    [shuffle_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
shuffle_sequences <- function(sequences, k = 1, method = "euler",
  nthreads = 1, rng.seed = sample.int(1e4, 1), window = FALSE, window.size = 0.1,
  window.overlap = 0.01) {

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

  # NOTE: Use sort_unique_cpp() instead of sort() at ALL TIMES when sorting a
  #       character vector. Otherwise VERY annoying bugs occur when the sequences
  #       contain a mix of upper and lower case letters.
  #       Update: sort(..., method = "radix") performs the same sorting as
  #       Rcpp::sort_unique()! Turns out the issue is that R picks "shell" as
  #       the sorting method for character vectors by default. Might as well
  #       stick with sort_unique_cpp() though, since it's 10X faster than
  #       sort(..., method = "radix").
  #       Maybe consider adding this to utils-internal.R:
  #       sort <- function(x, decreasing = FALSE, method = "radix", ...) {
  #         base::sort(x, decreasing = decreasing, method = method, ...)
  #       }

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
  char_check <- check_fun_params(list(method = args$method),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(k = args$k,
                                     nthreads = args$nthreads,
                                     rng.seed = args$rng.seed),
                                     numeric(), logical(), TYPE_NUM)
  s4_check <- check_fun_params(list(sequences = args$sequences),
                               numeric(), logical(), TYPE_S4)
  all_checks <- c(all_checks, char_check, num_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  alph <- seqtype(sequences)

  seq.names <- names(sequences)

  sequences <- as.character(sequences)
  seed <- as.integer(abs(rng.seed))
  if (nthreads < 0) stop("'nthreads' cannot be less than 0")
  nthreads <- as.integer(nthreads)
  if (k < 1) stop("'k' must be greater than 0")
  k <- as.integer(k)

  if (window) {

    if (window.size <= 0)
      stop("`window.size` must be greater than zero")
    if (window.overlap < 0)
      stop("`window.overlap` must be greater than or equal to zero")

    window.size <- rep_len(window.size, length(sequences))
    window.overlap <- rep_len(window.overlap, length(sequences))

    window.size[window.size < 1] <- nchar(sequences)[window.size < 1] *
      window.size[window.size < 1]
    window.overlap[window.overlap < 1] <- nchar(sequences)[window.overlap < 1] *
      window.overlap[window.overlap < 1]

    window.size <- as.integer(window.size)
    window.overlap <- as.integer(window.overlap)

    sequences <- shuffle_local(sequences, k, method, nthreads, rng.seed,
      window.size, window.overlap)

  } else {

    if (k == 1) {
      sequences <- shuffle_k1_cpp(sequences, nthreads, seed)
    } else {
      sequences <- switch(method,
                           "euler" = shuffle_euler_cpp(sequences, k, nthreads, seed),
                           "markov" = shuffle_markov_cpp(sequences, k, nthreads, seed),
                           "linear" = shuffle_linear_cpp(sequences, k, nthreads, seed)
                         )
    }

  }

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

# this function won't leave any 'leftover' letters behind, but the k-lets are
# predictable; the sequence is split up linearly every 'k'-letters
shuffle_linear <- function(sequence, k, mode = 1) {
  # - benchmark timings: about 1.25 times slower than shuffle_k1
  # - runtime is independent of k!

  # somehow the safeExplode() + collapse_cpp() combo is faster than
  # utf8ToInt() + intToUtf8()!

  if (mode == 1) {
    seq1 <- safeExplode(sequence)
    seq.len <- length(seq1)
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

#' Get k-let frequencies.
#'
#' @param seqs1 CHAR Split sequence, usually from strsplit(seq, "")[[1]].
#' @param k INT k-let size.
#' @param to.return CHAR Return k-let frequencies and/or transition matrix.
#' @param as.prob BOOL Return k-let counts or probabilities.
#' @param alph CHAR Alphabet letters, split.
#'
#' @return List, with entries 'transitions' (matrix), 'frequencies' (data.frame),
#'    and/or 'counts' (data.frame).
#'
#' @noRd
letter_freqs <- function(string, k = 1, to.return = c("freqs", "trans"),
                          as.prob = TRUE, alph = NULL) {

  # The output for this is a mess, but a lot of functions depend on it
  # currently

  lets.inseq <- sort_unique_cpp(string)
  if (is.null(alph)) lets.uniq <- lets.inseq
  else lets.uniq <- sort_unique_cpp(alph)

  counts <- count_klets_cpp(collapse_cpp(string), k, 1)[[1]]
  klets1 <- get_klets_cpp(lets.inseq, k)
  freqs <- data.frame(lets = klets1, counts = counts, stringsAsFactors = FALSE)

  if (!is.null(alph)) {
    klets2 <- get_klets_cpp(alph, k)
    if (length(klets1) != length(klets2)) {
      freqs$order <- seq_len(nrow(freqs))
      freqs <- merge(data.frame(lets = klets2), freqs, all.x = TRUE,
                     by = "lets", sort = FALSE)
      freqs$counts[is.na(freqs$counts)] <- 0
      freqs <- freqs[order(freqs$order), ]
      freqs <- freqs[, -which(colnames(freqs) == "order")]
    }
  }

  out <- list()

  if ("freqs" %in% to.return) {
    if (as.prob) out$freqs <- data.frame(lets = freqs$lets,
                                         frequencies = freqs$counts / sum(freqs$counts),
                                         stringsAsFactors = FALSE)
    else out$counts <- freqs
  }

  if ("trans" %in% to.return) {
    trans <- t(matrix(freqs$counts, nrow = length(lets.uniq)))
    mlets <- get_klets_cpp(lets.uniq, k - 1)
    colnames(trans) <- lets.uniq
    rownames(trans) <- mlets
    if (as.prob) {
      trans.sums <- rowSums(trans)
      for (i in seq_len(nrow(trans))) trans[i, ] <- trans[i, ] / trans.sums[i]
    }
    trans[is.nan(trans)] <- 0
    out$transitions <- trans
  }

  out

}

shuffle_local <- function(seqs, k, method, nthreads, rng.seed, win, ov) {

  if (k == 1)
    method <- 4
  else
    method <- switch(method, euler = 1, markov = 2, linear = 3, stop("unknown method"))

  lens <- nchar(seqs)
  win <- rep_len(win, length(lens))
  ov <- rep_len(ov, length(lens))

  pos <- mapply(calc_wins, lens, win, ov, SIMPLIFY = FALSE)
  starts <- lapply(pos, function(x) x[[1]])
  stops <- lapply(pos, function(x) x[[2]])

  shuffle_seq_local_cpp(seqs, k, nthreads, rng.seed, starts, stops, method)

}
