#' Calculate sequence background.
#'
#' For a set of input sequences, calculate the overall sequence background for
#' any k-let size. Only recommended for non-DNA/RNA sequences: otherwise use
#' the much faster and more efficient [Biostrings::oligonucleotideFrequency()].
#'
#' @param sequences \code{\link{XStringSet}} Input sequences. Note that if
#'    multiple sequences are present, they will be combined into one.
#' @param k `integer` Size of k-let. Background can be calculated for any
#'    k-let size.
#' @param as.prob `logical(1)` Whether to return k-let counts or probabilities.
#' @param pseudocount `integer(1)` Add a count to each possible k-let. Prevents
#'    any k-let from having 0 or 1 probabilities.
#' @param alphabet `character` Provide a custom alphabet to calculate a
#'    background for. If `NULL`, then standard letters will be assumed for
#'    DNA, RNA and AA sequences, and all unique letters found will be used
#'    for `BStringSet` type sequences.
#' @param to.meme If not `NULL`, then [get_bkg()] will return the sequence
#'    background in MEME Markov Background Model format. Input for this argument
#'    will be used for `cat(..., file = to.meme)` within [get_bkg()]. See
#'    \url{http://meme-suite.org/doc/bfile-format.html} for a description of
#'    the format.
#' @param RC `logical(1)` Calculate the background of the reverse complement
#'    of the input sequences as well. Only valid for DNA/RNA.
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [shuffle_sequences()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for large jobs. Furthermore,
#'    the behaviour of `progress = TRUE` is changed if `BP = TRUE`; the default
#'    \pkg{BiocParallel} progress bar will be shown (which unfortunately is much
#'    less informative).
#'
#' @return
#'    If `to.meme = NULL`: a list with each entry being a named numeric vector
#'    for every element in `k`. Otherwise: `NULL`, invisibly.
#'
#' @examples
#' ## Compare to Biostrings version
#' library(Biostrings)
#' seqs.DNA <- create_sequences()
#' bkg.DNA <- get_bkg(seqs.DNA, k = 3, as.prob = FALSE)[[1]]
#' bkg.DNA2 <- oligonucleotideFrequency(seqs.DNA, 3, 1, as.prob = FALSE)
#' bkg.DNA2 <- colSums(bkg.DNA2)
#' all(bkg.DNA == bkg.DNA2)
#'
#' ## Create a MEME background file
#' get_bkg(seqs.DNA, k = 1:3, to.meme = stdout(), pseudocount = 1)
#'
#' ## Non-DNA/RNA/AA alphabets
#' seqs.QWERTY <- create_sequences("QWERTY")
#' bkg.QWERTY <- get_bkg(seqs.QWERTY, k = 1:2)
#'
#' @references
#'    \insertRef{meme3}{universalmotif}
#'
#' @seealso [create_sequences()], [scan_sequences()], [shuffle_sequences()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
get_bkg <- function(sequences, k = 1:3, as.prob = TRUE, pseudocount = 0,
                    alphabet = NULL, to.meme = NULL, RC = FALSE,
                    progress = FALSE, BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  s4_check <- check_fun_params(list(sequences = args$sequences),
                               numeric(), logical(), "S4")
  num_check <- check_fun_params(list(k = args$k, pseudocount = args$pseudocount),
                                c(0, 1), c(FALSE, FALSE), "numeric")
  char_check <- check_fun_params(list(alphabet = args$alphabet), numeric(),
                                 TRUE, "character")
  logi_check <- check_fun_params(list(as.prob = args$as.prob, RC = args$RC,
                                      progress = args$progress, BP = args$BP),
                                 numeric(), logical(), "logical")
  all_checks <- c(all_checks, char_check, num_check, s4_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  k <- as.integer(k)
  if (RC && (is(sequences, "DNAStringSet") || is(sequences, "RNAStringSet")))
    sequences <- c(sequences, reverseComplement(sequences))

  if (!is.null(to.meme)) {
    if (!all(k == seq_len(k[length(k)])))
      stop(wmsg("To create a MEME background file, `all(k == seq_along(k))` must be true"))
  }

  no.alph <- FALSE
  if (is.null(alphabet)) {
    if (is(sequences, "DNAStringSet")) alphabet <- DNA_BASES
    else if (is(sequences, "RNAStringSet")) alphabet <- RNA_BASES
    else if (is(sequences, "AAStringSet")) alphabet <- AA_STANDARD
    else if (is(sequences, "BStringSet")) {
      no.alph <- TRUE
    } else if (!is(sequences, "XStringSet")) {
      stop("`sequences` must be an `XStringSet` object")
    }
  } else alphabet <- sort(alphabet)

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- as.character(seq_len(length(sequences)))
  seqs <- as.character(sequences)
  seqs <- lapply(seqs, safeExplode)

  # This function can be made many times faster for processing multiple
  # sequences by combining them into. However, this introduces a big
  # problem: combining sequences means creating new k-lets between them.
  # (this only matters for k > 1)
  #
  # It's probably faster to find out what these new k-lets are, then subtracting,
  # rather than going back to repeating the letter_freqs() function for every
  # sequence.
  # Update: turns out this approach is actually slower. Best to go back to
  # working on sequences individually.

  if (no.alph) alphabet <- sort(unique(do.call(c, lapply(seqs, unique))))

  check.k1 <- FALSE
  if (1 %in% k) {
    check.k1 <- TRUE
    k <- k[k != 1]
    k1 <- vector("list", length(seqs))
    for (i in seq_along(k1)) {
      k1[[i]] <- as.numeric(table(seqs[[i]]))
    }
    k1 <- do.call(cbind, k1)
    k1 <- rowSums(k1)
    names(k1) <- alphabet
  }

  counts <- list()

  if (length(k) > 0) {
    for (i in seq_along(k)) {
      counts[[as.character(k[i])]] <- lapply_(seqs, get_counts, k = k[i],
                                              alph = alphabet, BP = BP,
                                              PB = progress)
      counts[[as.character(k[i])]] <- do.call(cbind, counts[[as.character(k[i])]])
      counts[[as.character(k[i])]] <- rowSums(counts[[as.character(k[i])]])
      names(counts[[as.character(k[i])]]) <- get_lets(alphabet, k[i])
    }
  }

  if (check.k1) {
    counts[["1"]] <- k1
  }

  if (pseudocount > 0) counts <- lapply(counts, function(x) x + 1)
  if (as.prob) counts <- lapply(counts, function(x) x / sum(x))
  counts <- counts[sort(names(counts))]

  if (!is.null(to.meme)) {

    if (pseudocount < 1) {
      zero.check <- vapply(counts, function(x) any(x == 0), logical(1))
      one.check <- vapply(counts, function(x) any(x == 1), logical(1))
      if (any(zero.check) || any(one.check))
        stop(wmsg("MEME background files do not allow 0 or 1 values, ",
                  "please try again with a `pseudocount` higher than 0"))
    }
    if (!as.prob) counts <- lapply(counts, function(x) x / sum(x))
    out <- character(0)
    if (check.k1) k <- c(1, k)
    out <- to_meme_bkg(counts)
    cat(out, sep = "\n", file = to.meme)

  } else {

    counts

  }

}

get_counts <- function(x, k, alph) {
  a <- letter_freqs(x, k, to.return = "freqs", as.prob = FALSE, alph = alph)
  a$counts$counts
}

get_lets <- function(alph, k) {
  lets <- expand.grid(rep(list(alph), k), stringsAsFactors = FALSE)
  sort(collapse_rows_df(lets))
}

to_meme_bkg <- function(counts) {
  count.names <- lapply(counts, names)
  meme.order <- seq(0, length(counts) - 1)
  out <- vector("list", length(counts))
  for (i in seq_along(meme.order)) {
    out[[i]] <- to_meme_bkg_single(counts[[i]], meme.order[i],
                                   names(counts[[i]]))
  }
  do.call(c, out)
}

to_meme_bkg_single <- function(counts, meme.order, count.names) {
  out <- character(length(counts) + 1)
  out[1] <- paste0("#   order ", meme.order)
  for (i in seq_along(counts)) {
    out[i + 1] <- paste0(format(count.names[i], width = 8),
                         formatC(counts[i], format = "e", digits = 3),
                         " ")
  }
  out
}
