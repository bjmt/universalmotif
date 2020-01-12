#' Create random sequences.
#'
#' Generate random sequences from any set of characters, represented as
#' \code{\link{XStringSet}} objects.
#'
#' @param alphabet `character(1)` One of `c('DNA', 'RNA', 'AA')`, or a string of
#'    characters to be used as the alphabet.
#' @param seqnum `numeric(1)` Number of sequences to generate.
#' @param seqlen `numeric(1)` Length of random sequences.
#' @param freqs `numeric` A named vector of probabilities. The length of the
#'    vector must be the power of the number of letters in the sequence alphabet.
#' @param nthreads `numeric(1)` Run [create_sequences()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'    Note that no speed up will occur for jobs with `seqnum = 1`.
#' @param rng.seed `numeric(1)` Set random number generator seed. Since sequence
#'    creation can occur simultaneously in multiple threads using C++, it cannot
#'    communicate with the regular `R` random number generator state and thus requires
#'    an independent seed. Each individual sequence creation instance is
#'    given the following seed: `rng.seed * index`. The default is to pick a random
#'    number as chosen by [sample()], which effectively is make [create_sequences()]
#'    dependent on the R RNG state.
#'
#' @return \code{\link{XStringSet}} The returned sequences are _unnamed_.
#'
#' @examples
#' ## create DNA sequences with slightly increased AT content:
#' sequences <- create_sequences(freqs = c(A=0.3, C=0.2, G=0.2, T=0.3))
#' ## create custom sequences:
#' sequences.QWER <- create_sequences("QWER")
#' ## you can include non-alphabet characters are well, even spaces:
#' sequences.custom <- create_sequences("!@#$ ")
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [create_motif()], [shuffle_sequences()]
#' @export
create_sequences <- function(alphabet = "DNA", seqnum = 100, seqlen = 100,
                             freqs, nthreads = 1,
                             rng.seed = sample.int(1e4, 1)) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(alphabet = args$alphabet),
                                 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(seqnum = args$seqnum,
                                     seqlen = args$seqlen,
                                     freqs = args$freqs,
                                     nthreads = args$nthreads,
                                     rng.seed = args$rng.seed),
                                c(1, 1, rep(0, 3)), c(FALSE, FALSE, rep(TRUE, 3)),
                                TYPE_NUM)
  all_checks <- c(char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  alph.letters <- switch(alphabet,
                         "DNA" = DNA_BASES,
                         "RNA" = RNA_BASES,
                         "AA"  = AA_STANDARD2,
                                 sort_unique_cpp(safeExplode(alphabet)))

  if (missing(freqs)) {
    freqs <- rep(1 / length(alph.letters), length(alph.letters))
    names(freqs) <- alph.letters
  } else {
    if (is.null(names(freqs))) stop("freqs must be NAMED vector")
  }

  freqs <- freqs[order(names(freqs))]
  k <- logb(length(freqs), length(alph.letters))
  if (k %% 1 != 0)
    stop(wmsg("The length of `freqs` must be the power of the number of letters ",
              "in the sequence alphabet"))

  trans <- if (k > 1) matrix(freqs, nrow = length(alph.letters)) else matrix()

  seqs <- create_sequences_cpp(seqlen, seqnum, alph.letters, k, freqs, nthreads,
                               rng.seed, trans)

  seqs <- switch(alphabet,
                 "DNA" = DNAStringSet(seqs),
                 "RNA" = RNAStringSet(seqs),
                 "AA"  = AAStringSet(seqs),
                         BStringSet(seqs))

  seqs

}

check_k_lets <- function(alph.letters, freqs, k) {
  lets1 <- names(freqs)
  lets2 <- get_klets(alph.letters, k)
  if (!isTRUE(all.equal(lets1, lets2, use.names = FALSE)))
    stop(wmsg("For a k-let size of ", k, ",",
              "probabilities should be provided for:\n",
              paste(lets2, collapse = " ")))
  invisible(NULL)
}
