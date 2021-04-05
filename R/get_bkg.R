#' Calculate sequence background.
#'
#' For a set of input sequences, calculate the overall sequence background for
#' any k-let size. For very large sequences DNA and RNA sequences (in the billions of bases),
#' please be aware of the much faster and more efficient
#' \code{\link[Biostrings:nucleotideFrequency]{Biostrings::oligonucleotideFrequency()}}.
#' [get_bkg()] can still be used in these cases, though it may take several seconds or
#' minutes to calculate the results (depending on requested k-let sizes).
#'
#' @param sequences \code{\link{XStringSet}} Input sequences. Note that if
#'    multiple sequences are present, the results will be combined into one.
#' @param k `integer` Size of k-let. Background can be calculated for any
#'    k-let size.
#' @param as.prob Deprecated.
#' @param pseudocount `integer(1)` Add a count to each possible k-let. Prevents
#'    any k-let from having 0 or 1 probabilities.
#' @param alphabet `character(1)` Provide a custom alphabet to calculate a
#'    background for. If `NULL`, then standard letters will be assumed for
#'    DNA, RNA and AA sequences, and all unique letters found will be used
#'    for `BStringSet` type sequences. Note that letters which are not a part
#'    of the standard DNA/RNA/AA alphabets or in the provided alphabet will
#'    not be counted in the totals during probability calculations.
#' @param to.meme If not `NULL`, then [get_bkg()] will return the sequence
#'    background in MEME Markov Background Model format. Input for this argument
#'    will be used for `cat(..., file = to.meme)` within [get_bkg()]. See
#'    \url{http://meme-suite.org/doc/bfile-format.html} for a description of
#'    the format.
#' @param RC `logical(1)` Calculate the background of the reverse complement
#'    of the input sequences as well. Only valid for DNA/RNA.
#' @param list.out Deprecated.
#' @param nthreads `numeric(1)` Run [get_bkg()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'    Note that no speed up will occur for jobs with only a single sequence.
#' @param merge.res `logical(1)` Whether to merge results from all sequences
#'    or return background data for individual sequences.
#' @param window `logical(1)` Determine background in windows.
#' @param window.size `numeric` Window size. If a number between 0 and 1 is
#'    provided, the value is calculated as the number multiplied by the sequence
#'    length.
#' @param window.overlap `numeric` Overlap between windows. If a number
#'    between 0 and 1 is provided, the value is calculated as the number
#'    multiplied by the sequence length.
#'
#' @return
#'    If `to.meme = NULL`, a `DataFrame` with columns `klet`, `count`,
#'    and `probability`. If `merge.res = FALSE`, there will be an additional
#'    `sequence` column. If `window = TRUE`, there will be an additional `start`
#'    and `stop` columns.
#'
#'    If `to.meme` is not `NULL`, then `NULL` is returned, invisibly.
#'
#' @examples
#' ## Compare to Biostrings version
#' library(Biostrings)
#' seqs.DNA <- create_sequences()
#' bkg.DNA <- get_bkg(seqs.DNA, k = 3)
#' bkg.DNA2 <- oligonucleotideFrequency(seqs.DNA, 3, 1, as.prob = FALSE)
#' bkg.DNA2 <- colSums(bkg.DNA2)
#' all(bkg.DNA$count == bkg.DNA2)
#'
#' ## Create a MEME background file
#' get_bkg(seqs.DNA, k = 1:3, to.meme = stdout(), pseudocount = 1)
#'
#' ## Non-DNA/RNA/AA alphabets
#' seqs.QWERTY <- create_sequences("QWERTY")
#' bkg.QWERTY <- get_bkg(seqs.QWERTY, k = 1:2)
#'
#' @references
#'
#' Bailey TL, Elkan C (1994). “Fitting a mixture model by expectation
#' maximization to discover motifs in biopolymers.” *Proceedings of
#' the Second International Conference on Intelligent Systems for
#' Molecular Biology*, **2**, 28-36.
#'
#' @seealso [create_sequences()], [scan_sequences()], [shuffle_sequences()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
get_bkg <- function(sequences, k = 1:3, as.prob = NULL, pseudocount = 0,
  alphabet = NULL, to.meme = NULL, RC = FALSE, list.out = NULL, nthreads = 1,
  merge.res = TRUE, window = FALSE, window.size = 0.1, window.overlap = 0) {

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
  logi_check <- check_fun_params(list(RC = args$RC),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(all_checks, char_check, num_check, s4_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (!is.null(as.prob))
    warning(wmsg("`as.prob` has been disabled and now does nothing ",
        "(both counts and probabilities are now shown)"), immediate. = TRUE)
  if (!is.null(list.out))
    warning(wmsg("`list.out` has been disabled and now does nothing ",
        "(all posisble outputs are now returned as DataFrames)"), immediate. = TRUE)

  if (!is.null(to.meme) && window)
    stop("`window = TRUE` is not valid if `to.meme` is not `NULL`")
  if (!is.null(to.meme) && !merge.res)
    stop("`merge.res = FALSE` is not valid if `to.meme` is not `NULL`")

  pseudocount <- as.integer(pseudocount)
  k <- as.integer(k)
  if (RC && seqtype(sequences) %in% c("DNA", "RNA"))
    sequences <- c(sequences, reverseComplement(sequences))

  if (!is.null(to.meme)) {
    if (!all(k == seq_len(k[length(k)])))
      stop(wmsg("To create a MEME background file, `all(k == seq_along(k))` must be true"))
  }

  seq_type <- seqtype(sequences)

  no.alph <- FALSE
  if (is.null(alphabet)) {
    if (!is(sequences, "XStringSet"))
      stop("`sequences` must be an `XStringSet` object")
    switch(seq_type,
           "DNA" = alphabet <- DNA_BASES,
           "RNA" = alphabet <- RNA_BASES,
           "AA" = alphabet <- AA_STANDARD2,
           no.alph <- TRUE)
  } else {
    if (length(alphabet) == 1) alphabet <- sort_unique_cpp(safeExplode(alphabet))
    else alphabet <- sort_unique_cpp(alphabet)
  }

  seqlens <- width(sequences)
  if (min(seqlens) < max(k))
    stop("`k` must be less than the length of the shortest sequence")

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- as.character(seq_len(length(sequences)))
  seqs1 <- as.character(sequences)
  seqs <- lapply(seqs1, safeExplode)

  if (no.alph) alphabet <- sort_unique_cpp(do.call(c, lapply(seqs, unique)))
  alph <- collapse_cpp(alphabet)

  if (!window) {

    counts <- vector("list", length(k))
    names(counts) <- as.character(k)
    for (i in seq_along(k)) {
      counts[[as.character(k[i])]] <- count_klets_alph_cpp(seqs1, alph, k[i], nthreads)
      counts[[as.character(k[i])]] <- do.call(data.frame, counts[[as.character(k[i])]])
      if (merge.res) {
        counts[[as.character(k[i])]] <- rowSums(counts[[as.character(k[i])]])
        names(counts[[as.character(k[i])]]) <- get_klets(alphabet, k[i])
      } else {
        colnames(counts[[as.character(k[i])]]) <- seq.names
        rownames(counts[[as.character(k[i])]]) <- get_klets(alphabet, k[i])
        counts[[as.character(k[i])]] <- as.matrix(t(counts[[as.character(k[i])]]))
      }
    }

    if (pseudocount > 0) counts <- lapply(counts, function(x) x + pseudocount)
    if (merge.res) {
      probs <- lapply(counts, function(x) x / sum(x))
    } else {
      probs <- lapply(counts, function(x) t(apply(x, 1, function(y) y / sum(y))))
    }
    counts <- counts[sort_unique_cpp(names(counts))]
    probs <- probs[sort_unique_cpp(names(probs))]

    if (!is.null(to.meme)) {

      if (pseudocount < 1) {
        zero.check <- vapply(probs, function(x) any(x == 0), logical(1))
        one.check <- vapply(probs, function(x) any(x == 1), logical(1))
        if (any(zero.check) || any(one.check))
          stop(wmsg("MEME background files do not allow 0 or 1 values, ",
                    "please try again with a `pseudocount` higher than 0"))
      }
      out <- character(0)
      out <- to_meme_bkg(probs)
      cat(out, sep = "\n", file = to.meme)
      return(invisible())

    }

    if (merge.res) {
      counts <- unlist(counts)
      names(counts) <- vapply(strsplit(names(counts), ".", fixed = TRUE),
        function(x) x[2], character(1))
      probs <- unlist(probs)
      names(probs) <- vapply(strsplit(names(probs), ".", fixed = TRUE),
        function(x) x[2], character(1))
      res <- data.frame(names(counts), unname(counts), unname(probs), row.names = NULL)
      colnames(res) <- c("klet", "count", "probability")
      res <- as(res, "DataFrame")
    } else {
      for (i in seq_along(counts)) {
        counts[[i]] <- stack(counts[[i]])
        probs[[i]] <- stack(probs[[i]])
      }
      probs <- do.call(rbind, probs)
      counts <- do.call(rbind, counts)
      counts[[1]] <- as.character(counts[[1]])
      counts[[2]] <- as.character(counts[[2]])
      res <- counts
      colnames(res) <- c("sequence", "klet", "count")
      res$probability <- probs[[3]]
    }

  } else {

    window.size <- rep_len(window.size, length(sequences))
    window.overlap <- rep_len(window.overlap, length(sequences))

    if (any(window.size == 0))
      stop("`window.size` must be greater than 0")

    window.size[window.size < 1] <- as.integer(window.size[window.size < 1] *
      seqlens[window.size < 1])
    window.overlap[window.overlap < 1 & window.overlap > 0] <- as.integer(
      window.overlap[window.overlap < 1 & window.overlap > 0] *
      seqlens[window.overlap < 1 & window.overlap > 0])

    if (any(window.size <= window.overlap))
      stop("`window.overlap` cannot be larger than or equal to `window.size`")

    wins <- mapply(calc_wins, seqlens, window.size, window.overlap,
      SIMPLIFY = FALSE)
    starts <- lapply(wins, function(x) x$starts)
    stops <- lapply(wins, function(x) x$stops)

    seqs.split <- mapply(split_seq_by_win, seqs1, starts, stops,
      SIMPLIFY = FALSE)

    split_n <- vapply(seqs.split, length, integer(1))
    seqnames_n <- mapply(function(x, y) rep(x, y), seq.names, split_n, SIMPLIFY = FALSE)
    seqnames_n <- do.call(c, seqnames_n)

    seqs.split2 <- do.call(c, unname(seqs.split))
    if (min(nchar(seqs.split2)) < max(k))
      stop("`k` must be less than the smallest `window.size`")

    starts <- do.call(c, starts)
    stops <- do.call(c, stops)

    res <- vector("list", length(k))
    for (i in seq_along(res)) {
      res[[i]] <- count_klets_alph_cpp(seqs.split2, alph, k[i], nthreads)
      res[[i]] <- do.call(data.frame, res[[i]])
      rownames(res[[i]]) <- get_klets(alphabet, k[i])
      colnames(res[[i]]) <- NULL
      res[[i]] <- stack(as.matrix(res[[i]]))
      colnames(res[[i]]) <- c("klet", "col", "count")
      res[[i]]$klet <- as.character(res[[i]]$klet)
      res[[i]]$sequence <- seqnames_n[as.integer(res[[i]]$col)]
      res[[i]]$start <- starts[as.integer(res[[i]]$col)]
      res[[i]]$stop <- stops[as.integer(res[[i]]$col)]
    }
    res <- do.call(rbind, res)
    res <- res[, c("sequence", "start", "stop", "klet", "count")]

    if (merge.res) {
      if (length(unique(seqlens)) != 1)
        stop(wmsg("`merge.res = TRUE` and `window = TRUE` is only valid if ",
            "all sequences are of equal length"))
      res <- aggregate(count ~ start + stop + klet, data = res[, -1], sum)
      if (pseudocount > 0) res$count <- res$count + pseudocount
      res <- by(res, INDICES = list(res$start, nchar(res$klet)), FUN = function(x) {
        x$probability <- x$count / sum(x$count) ; x
      })
      res <- do.call(rbind, res)
      rownames(res) <- NULL
      res <- res[, c("start", "stop", "klet", "count", "probability")]
    } else {
      if (pseudocount > 0) res$count <- res$count + pseudocount
      res <- by(res, INDICES = list(res$sequence, res$start, nchar(res$klet)),
        FUN = function(x) {
          x$probability <- x$count / sum(x$count) ; x
        }
      )
      res <- do.call(rbind, res)
      rownames(res) <- NULL
      res <- res[, c("sequence", "start", "stop", "klet", "count", "probability")]
    }

    if (!merge.res)
      res <- res[order(res$sequence, res$klet, method = "radix"), ]
    else
      res <- res[order(res$klet, method = "radix"), ]

    res$start <- as.integer(res$start)
    res$stop <- as.integer(res$stop)

    res <- as(res, "DataFrame")

  }

  res$probability[is.nan(res$probability)] <- 0

  res

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

calc_wins <- function(seqlen, winsize, ovrlp) {
  starts <- seq(from = 1, to = seqlen, by = winsize)
  if (ovrlp == 0) {
    stops <- c(c(1, starts[c(-1, -length(starts))]) + winsize - 1, seqlen)
    list(starts = starts, stops = stops)
  } else {
    starts[-1] <- starts[-1] - ovrlp
    starts <- c(starts, starts[length(starts)] + winsize)
    if (starts[length(starts)] >= seqlen) starts <- starts[-length(starts)]
    stops <- c((starts[-length(starts)] + winsize) - 1, seqlen)
    list(starts = starts, stops = stops)
  }
}
