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
#' @param monofreqs `numeric` Deprecated. Use `freqs` instead.
#' @param difreqs `numeric` Deprecated. Use `freqs` instead.
#' @param trifreqs `numeric` Deprecated. Use `freqs` instead.
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [create_sequences()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for large jobs (such as
#'    `create_sequences(seqlen=100000,seqnum=100000)`). Furthermore,
#'    the behaviour of `progress = TRUE` is
#'    changed if `BP = TRUE`; the default \pkg{BiocParallel} progress bar will
#'    be shown (which unfortunately is much less informative).
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
#' @references
#'    \insertRef{biostrings}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [create_motif()], [shuffle_sequences()]
#' @export
create_sequences <- function(alphabet = "DNA", seqnum = 100, seqlen = 100,
                             freqs, monofreqs, difreqs, trifreqs,
                             progress = FALSE, BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(alphabet = args$alphabet),
                                 1, FALSE, "character")
  num_check <- check_fun_params(list(seqnum = args$seqnum, seqlen = args$seqlen,
                                     freqs = args$freqs,
                                     monofreqs = args$monofreqs,
                                     difreqs = args$difreqs,
                                     trifreqs = args$trifreqs),
                                c(1, 1, rep(0, 4)), c(FALSE, FALSE, rep(TRUE, 4)),
                                "numeric")
  logi_check <- check_fun_params(list(progress = args$progress, BP = args$BP),
                                 numeric(), logical(), "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  alph.letters <- switch(alphabet,
                         "DNA" = DNA_BASES,
                         "RNA" = RNA_BASES,
                         "AA"  = AA_STANDARD,
                                 sort(safeExplode(alphabet)))

  if (missing(freqs)) {
    if (!missing(monofreqs)) {
      warning("`monofreqs` option is deprecated, please use `freqs`",
              immediate. = TRUE)
      freqs <- monofreqs
    }
    if (!missing(difreqs)) {
      warning("`difreqs` option is deprecated, please use `freqs`",
              immediate. = TRUE)
      freqs <- difreqs
    }
    if (!missing(trifreqs)) {
      warning("`trifreqs` option is deprecated, please use `freqs`",
              immediate. = TRUE)
      freqs <- trifreqs
    }
  }


  if (missing(monofreqs) &&
      missing(difreqs) &&
      missing(trifreqs) &&
      missing(freqs)) {
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
  seqs <- vector("list", seqnum)

  if (k == 1) {
    seqs <- lapply_(seq_len(seqnum),
                     function(x) create_k1(alph.letters = alph.letters,
                                           seqlen = seqlen, bkg = freqs),
                    BP = BP, PB = progress)
  } else  {
    check_k_lets(alph.letters, freqs, k)
    seqs <- lapply_(seq_len(seqnum),
                     function(x) create_kany(alph.letters = alph.letters,
                                             seqlen = seqlen, freqs = freqs,
                                             k = k),
                    BP = BP, PB = progress)
  }

  seqs <- unlist(seqs)

  seqs <- switch(alphabet,
                 "DNA" = DNAStringSet(seqs),
                 "RNA" = RNAStringSet(seqs),
                 "AA"  = AAStringSet(seqs),
                         BStringSet(seqs))

  seqs

}

create_k1 <- function(alph.letters, seqlen, bkg) {
  seqout <- sample(alph.letters, seqlen, replace = TRUE, prob = bkg)
  seqout <- collapse_cpp(seqout)
  seqout
}

check_k_lets <- function(alph.letters, freqs, k) {
  lets1 <- names(freqs)
  lets2 <- expand.grid(rep(list(alph.letters), k), stringsAsFactors = FALSE)
  lets2 <- sort(collapse_rows_df(lets2))
  if (!isTRUE(all.equal(lets1, lets2, use.names = FALSE)))
    stop(wmsg("For a k-let size of ", k, ",",
              "probabilities should be provided for:\n",
              paste(lets2, collapse = " ")))
  invisible(NULL)
}

create_kany <- function(alph.letters, seqlen, freqs, k) {

  seqout <- character(seqlen)
  lets <- names(freqs)
  first.k <- sample(lets, 1, prob = freqs)
  first.k <- safeExplode(first.k)
  seqout[seq_len(k)] <- first.k

  trans <- t(matrix(freqs, nrow = length(alph.letters)))
  lets2 <- expand.grid(rep(list(alph.letters), k - 1), stringsAsFactors = FALSE)
  lets2 <- collapse_rows_df(lets2)
  colnames(trans) <- alph.letters
  rownames(trans) <- lets2

  trans <- t(trans)
  trans[trans == 0] <- 1 / nrow(trans) / 1000
  names(lets2) <- lets2

  seqout <- shuffle_markov_loop(k, seqlen, k, seqout, alph.letters, trans, lets2)

  seqout

}
