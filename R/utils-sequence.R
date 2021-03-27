#' Sequence-related utility functions.
#'
#' @param alph `character(1)` A single character string with the desired
#'    sequence alphabet. If missing, finds the unique letters in the string.
#' @param lets `character` A character vector where each element will be
#'    considered a single unit.
#' @param k `integer(1)` K-let size.
#' @param letter `character(1)` Character to use for masking.
#' @param pattern `character(1)` Pattern to mask.
#' @param seqs `XStringSet` Sequences to mask. Cannot be `BStringSet`.
#' @param string `character(1)` A length one character vector.
#' @param method `character(1)` Shuffling method. One of `c("euler", "linear",
#'    "markov")`. See [shuffle_sequences()].
#' @param RC `logical(1)` Whether to mask the reverse complement of the pattern.
#' @param rng.seed `numeric(1)` Set random number generator seed. Since shuffling
#'    in [shuffle_sequences()] can occur simultaneously in multiple threads using C++,
#'    it cannot communicate
#'    with the regular `R` random number generator state and thus requires an
#'    independent seed. Since [shuffle_string()] uses the same underlying code
#'    as [shuffle_sequences()], it also requires a separate seed even if it is
#'    run in serial.
#'
#' @return
#'    For [count_klets()]: A `data.frame` with columns `lets` and `counts`.
#'
#'    For [get_klets()]: A `character` vector of k-lets.
#'
#'    For [mask_seqs()]: The masked `XStringSet` object.
#'
#'    For [shuffle_string()]: A single `character` string.
#'
#' @examples
#' #######################################################################
#' ## count_klets
#' ## Count k-lets for any string of characters
#' count_klets("GCAAATGTACGCAGGGCCGA", k = 2)
#' ## The default 'k' value (1) counts individual letters
#' count_klets("GCAAATGTACGCAGGGCCGA")
#'
#' #######################################################################
#' ## get_klets
#' ## Generate all possible k-lets for a set of characters
#' get_klets(c("A", "C", "G", "T"), 3)
#' ## Note that each element in 'lets' is considered a single unit;
#' ## see:
#' get_klets(c("AA", "B"), k = 2)
#'
#' #######################################################################
#' ## mask_seqs
#' ## Mask repetitive seqeuences
#' data(ArabidopsisPromoters)
#' mask_seqs(ArabidopsisPromoters, "AAAAAA")
#'
#' #######################################################################
#' ## shuffle_string
#' ## Shuffle any string of characters
#' shuffle_string("ASDADASDASDASD", k = 2)
#'
#' @seealso [create_sequences()], [shuffle_sequences()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @name utils-sequence
NULL

#' @rdname utils-sequence
#' @export
count_klets <- function(string, k = 1, alph) {

  if (k < 1) stop("k must be greater than 0")
  k <- as.integer(k)

  if (length(string) != 1) stop("'string' must be a length 1 character vector")
  if (nchar(string) < 1) stop("'string' cannot be empty")

  if (missing(alph)) {
    counts <- count_klets_cpp(string, k, 1)[[1]]
    klets <- get_klets_cpp(sort_unique_cpp(safeExplode(string)), k)
  } else {
    if (length(alph) > 1) stop("'alph' must be a single string")
    if (nchar(alph) < 1) stop("'alph' cannot be empty")
    counts <- count_klets_alph_cpp(string, alph, k, 1)
    klets <- get_klets_cpp(sort_unique_cpp(safeExplode(alph)), k)
  }

  data.frame(klets, counts, stringsAsFactors = FALSE)

}

# get_klets(lets, k = 1) --> see utils.cpp

#' @rdname utils-sequence
#' @export
get_klets <- function(lets, k = 1) {

  get_klets_cpp(lets, k)

}

#' @rdname utils-sequence
#' @export
mask_seqs <- function(seqs, pattern, RC = FALSE, letter = "-") {
  if (!is(seqs, "XStringSet"))
    stop("`seqs` must be an `XStringSet` object")
  if (length(pattern) > 1 || !is.character(pattern))
    stop("`pattern` must be a single character")
  alph <- seqtype(seqs)
  if (alph == "B")
    stop("`mask_seqs()` only works with DNA/RNA/AA sequences")
  fix_seqs <- function(seqs, pattern, letter) {
    seqs <- lapply(seqs, mask, pattern = pattern)
    seqs <- lapply(seqs, injectHardMask, letter = letter)
    switch(alph,
      DNA = DNAStringSet(seqs),
      RNA = RNAStringSet(seqs),
      AA = AAStringSet(seqs)
    )
  }
  seqs <- fix_seqs(seqs, pattern, letter)
  if (RC) {
    pattern <- as.character(switch(alph,
        DNA = reverseComplement(DNAString(pattern)),
        RNA = reverseComplement(RNAString(pattern)),
        stop("`RC = TRUE` is only valid for DNA/RNA sequences")
    ))
    seqs <- fix_seqs(seqs, pattern, letter)
  }
  seqs
}

#' @rdname utils-sequence
#' @export
shuffle_string <- function(string, k = 1, method = c("euler", "linear", "markov"),
                           rng.seed = sample.int(1e4, 1)) {

  method <- match.arg(method, c("euler", "linear", "markov"))

  if (length(string) != 1) stop("'string' must be a length 1 character vector")

  if (length(k) != 1) stop("'k' must be length 1")
  if (k < 1) stop("'k' must be greater than 0")
  k <- as.integer(k)

  seed <- as.integer(abs(rng.seed))[1]

  if (k == 1) {

    shuffle_k1_cpp(string, 1, seed)

  } else if (k > 1) {

    switch(method,
           "euler" = shuffle_euler_cpp(string, k, 1, seed),
           "linear" = shuffle_linear_cpp(string, k, 1, seed),
           "markov" = shuffle_markov_cpp(string, k, 1, seed),
           stop("'method' must be one of 'euler', 'linear', 'markov'"))

  } else {

    stop("k must be greater than 0")

  }

}
