#' Calculate sequence background.
#'
#' For a set of input sequences, calculate the overall sequence background for
#' any k-let size. Only recommended for non-DNA/RNA sequences: otherwise use
#' the much faster and more efficient
#' \code{\link[Biostrings:nucleotideFrequency]{Biostrings::oligonucleotideFrequency()}}.
#'
#' @param sequences \code{\link{XStringSet}} Input sequences. Note that if
#'    multiple sequences are present, they will be combined into one.
#' @param k `integer` Size of k-let. Background can be calculated for any
#'    k-let size.
#' @param as.prob `logical(1)` Whether to return k-let counts or probabilities.
#' @param pseudocount `integer(1)` Add a count to each possible k-let. Prevents
#'    any k-let from having 0 or 1 probabilities.
#' @param alphabet `character(1)` Provide a custom alphabet to calculate a
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
#' @param list.out `logical(1)` Return background frequencies as list, with an
#'    entry for each `k`. If `FALSE`, return a single vector.
#' @param progress `logical(1)` Deprecated. Does nothing.
#' @param BP `logical(1)` Deprecated. See `nthreads`.
#' @param nthreads `numeric(1)` Run [get_bkg()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'    Note that no speed up will occur for jobs with only a single sequence.
#'
#' @return
#'    If `to.meme = NULL` and `list.out = TRUE`: a list with each entry being a
#'    named numeric vector for every element in `k`. If `to.meme = NULL` and
#'    `list.out = FALSE`: a named numeric vector. Otherwise: `NULL`, invisibly.
#'
#' @examples
#' ## Compare to Biostrings version
#' library(Biostrings)
#' seqs.DNA <- create_sequences()
#' bkg.DNA <- get_bkg(seqs.DNA, k = 3, as.prob = FALSE, list.out = FALSE)
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
                    list.out = TRUE, progress = FALSE, BP = FALSE,
                    nthreads = 1) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (any(k < 1)) {
    k_check <- paste0(" * Incorrect 'k': values below 1 are not allowed; found `",
                      paste0(k[k < 1], collapse = ", "), "`")
    all_checks <- c(all_checks, k_check)
  }
  s4_check <- check_fun_params(list(sequences = args$sequences),
                               numeric(), logical(), TYPE_S4)
  num_check <- check_fun_params(list(k = args$k, pseudocount = args$pseudocount),
                                c(0, 1), c(FALSE, FALSE), TYPE_NUM)
  char_check <- check_fun_params(list(alphabet = args$alphabet), 1,
                                 TRUE, TYPE_CHAR)
  logi_check <- check_fun_params(list(as.prob = args$as.prob, RC = args$RC,
                                      progress = args$progress, BP = args$BP,
                                      list.out = args$list.out),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(all_checks, char_check, num_check, s4_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (progress)
    warning("'progress' is deprecated and does nothing", immediate. = TRUE)
  if (BP)
    warning("'BP' is deprecated, see 'nthreads'", immediate. = TRUE)

  k <- as.integer(k)
  if (RC && seqtype(sequences) %in% c("DNA", "RNA"))
    sequences <- c(sequences, reverseComplement(sequences))

  if (!is.null(to.meme)) {
    if (!all(k == seq_len(k[length(k)])))
      stop(wmsg("To create a MEME background file, `all(k == seq_along(k))` must be true"))
  }

  no.alph <- FALSE
  if (is.null(alphabet)) {
    if (!is(sequences, "XStringSet"))
      stop("`sequences` must be an `XStringSet` object")
    switch(seqtype(sequences),
           "DNA" = alphabet <- DNA_BASES,
           "RNA" = alphabet <- RNA_BASES,
           "AA" = alphabet <- AA_STANDARD,
           no.alph <- TRUE)
  } else {
    if (length(alphabet) == 1) alphabet <- sort_unique_cpp(safeExplode(alphabet))
    else alphabet <- sort_unique_cpp(alphabet)
  }

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- as.character(seq_len(length(sequences)))
  seqs1 <- as.character(sequences)
  seqs <- lapply(seqs1, safeExplode)

  # This function can be made many times faster for processing multiple
  # sequences by combining them into one. However, this introduces a big
  # problem: combining sequences means creating new k-lets between them.
  # (this only matters for k > 1)
  #
  # It's probably faster to find out what these new k-lets are, then subtracting,
  # rather than going back to repeating the letter_freqs() function for every
  # sequence.
  # Update: turns out this approach is actually slower. Best to go back to
  # working on sequences individually.

  if (no.alph) alphabet <- sort_unique_cpp(do.call(c, lapply(seqs, unique)))
  alph <- collapse_cpp(alphabet)

  counts <- vector("list", length(k))
  names(counts) <- as.character(k)
  for (i in seq_along(k)) {
    counts[[as.character(k[i])]] <- count_klets_alph_cpp(seqs1, alph, k[i], nthreads)
    counts[[as.character(k[i])]] <- do.call(data.frame, counts[[as.character(k[i])]])
    counts[[as.character(k[i])]] <- rowSums(counts[[as.character(k[i])]])
    names(counts[[as.character(k[i])]]) <- get_klets(alphabet, k[i])
  }

  if (pseudocount > 0) counts <- lapply(counts, function(x) x + 1)
  if (as.prob) counts <- lapply(counts, function(x) x / sum(x))
  counts <- counts[sort_unique_cpp(names(counts))]

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
    out <- to_meme_bkg(counts)
    cat(out, sep = "\n", file = to.meme)

  } else {

    if (!list.out) {
      counts <- unlist(counts)
      names(counts) <- vapply(names(counts),
                              function(x) strsplit(x, ".", fixed = TRUE)[[1]][2],
                              character(1))
    }

    counts

  }

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
