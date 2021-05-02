#' Scan sequences for matches to input motifs.
#'
#' For sequences of any alphabet, scan them using the PWM matrices of
#' a set of input motifs.
#'
#' @param motifs See `convert_motifs()` for acceptable motif formats.
#' @param sequences \code{\link{XStringSet}} Sequences to scan. Alphabet
#'    should match motif.
#' @param threshold `numeric(1)` See details.
#' @param threshold.type `character(1)` One of `c('logodds', 'logodds.abs',
#'    'pvalue')`. See details.
#' @param RC `logical(1)` If `TRUE`, check reverse complement of input
#'    sequences.
#' @param use.freq `numeric(1)` The default, 1, uses the motif matrix (from
#'    the `motif['motif']` slot) to search for sequences. If a higher
#'    number is used, then the matching k-let matrix from the
#'    `motif['multifreq']` slot is used. See [add_multifreq()].
#' @param verbose `numeric(1)` Describe progress, from none (`0`) to
#'    verbose (`3`).
#' @param nthreads `numeric(1)` Run [scan_sequences()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'    Note that no speed up will occur for jobs with only a single motif and
#'    sequence.
#' @param motif_pvalue.k `numeric(1)` Control [motif_pvalue()] approximation.
#'    See [motif_pvalue()].
#' @param use.gaps `logical(1)` Set this to `FALSE` to ignore motif gaps, if
#'    present.
#' @param allow.nonfinite `logical(1)` If `FALSE`, then apply a pseudocount if
#'    non-finite values are found in the PWM. Note that if the motif has a
#'    pseudocount greater than zero and the motif is not currently of type PWM,
#'    then this parameter has no effect as the pseudocount will be
#'    applied automatically when the motif is converted to a PWM internally. This
#'    value is set to `FALSE` by default in order to stay consistent with
#'    pre-version 1.8.0 behaviour.
#' @param warn.NA `logical(1)` Whether to warn about the presence of non-standard
#'    letters in the input sequence, such as those in masked sequences.
#' @param calc.pvals `logical(1)` Calculate P-values for each hit. This is a
#'    convenience option which simply gives `motif_pvalue()` the input motifs
#'    and the scores of each hit. Be careful about setting this to `TRUE` if
#'    you anticipate getting thousands of hits: expect to wait a few seconds or
#'    minutes for the calculations to finish. Increasing the `nthreads` value
#'    can help greatly here. See Details for more information on P-value
#'    calculation.
#' @param return.granges `logical(1)` Return the results as a `GRanges` object.
#'    Requires the `GenomicRanges` package to be installed.
#' @param no.overlaps `logical(1)` Remove overlapping hits from the same motifs.
#'    Overlapping hits from different motifs are preserved. Please note that the
#'    current implementation of this feature can add significantly to the run
#'    time for large inputs.
#' @param no.overlaps.by.strand `logical(1)` Whether to discard overlapping hits
#'    from the opposite strand, or to only discard overlapping hits on the same
#'    strand.
#' @param no.overlaps.strat `character(1)` One of `c("score", "order")`.
#'    The former option keeps the highest scoring overlapping hit (and the first
#'    of these within ties), and the latter simply keeps the first overlapping hit.
#'    keeps the highest scoring 
#' @param respect.strand `logical(1)` If `RC = TRUE` and motifs are DNA/RNA,
#'    then setting this option to `TRUE` will make `scan_sequences()` only
#'    scan the strands of the input sequences as indicated in the motif
#'    `strand` slot.
#'
#' @return `DataFrame` with each row representing one hit. If the input
#'    sequences are \code{\link{DNAStringSet}} or
#'    \code{\link{RNAStringSet}}, then an
#'    additional column with the strand is included. Function args are stored
#'    in the `metadata` slot.
#'
#' @details
#'    Similar to [Biostrings::matchPWM()], the scanning method uses
#'    logodds scoring. (To see the scoring matrix for any motif, simply
#'    run `convert_type(motif, "PWM")`. For a `multifreq` scoring
#'    matrix: `apply(motif["multifreq"][["2"]], 2, ppm_to_pwm)`). In order
#'    to score a sequence, at each position within a sequence of length equal
#'    to the length of the motif, the scores for each base are summed. If the
#'    score sum is above the desired threshold, it is kept.
#'
#'    If `threshold.type = 'logodds'`, then the `threshold` value is multiplied
#'    by the maximum possible motif scores. To calculate the
#'    maximum possible scores a motif (of type PWM) manually, run
#'    `motif_score(motif, 1)`. If \code{threshold.type = 'pvalue'},
#'    then threshold logodds scores are generated using [motif_pvalue()].
#'    Finally, if \code{threshold.type = 'logodds.abs'}, then the exact values
#'    provided will be used as thresholds.
#'
#'    Non-standard letters (such as "N", "+", "-", ".", etc in \code{\link{DNAString}}
#'    objects) will be safely ignored, resulting only in a warning and a very
#'    minor performance cost. This can used to scan
#'    masked sequences. See \code{\link[Biostrings:maskMotif]{Biostrings::mask()}}
#'    for masking sequences
#'    (generating \code{\link{MaskedXString}} objects), and [Biostrings::injectHardMask()]
#'    to recover masked \code{\link{XStringSet}} objects for use with [scan_sequences()].
#'    There is also a provided wrapper function which performs both steps: [mask_seqs()].
#'
#'    When `calc.pvals = TRUE`, [motif_pvalue()] will calculate the probabilities
#'    of getting the input scores or higher, which is why it can take time to
#'    calculate the P-values. If you simply wish to calculate the
#'    probabilities of getting individual matches based on background frequencies,
#'    then the following code can be used to achieve
#'    this (using the list of input motifs and [scan_sequences()] results):
#'    `mapply(prob_match, motifs[scanRes$motif.i], scanRes$match)`. Of course
#'    this only matters if you do not have uniform background frequencies, or
#'    else the probability of each match is simply `(1 / nrow(motif))^ncol(motif)`.
#'
#' @examples
#' ## any alphabet can be used
#' \dontrun{
#' set.seed(1)
#' alphabet <- paste(c(letters), collapse = "")
#' motif <- create_motif("hello", alphabet = alphabet)
#' sequences <- create_sequences(alphabet, seqnum = 1000, seqlen = 100000)
#' scan_sequences(motif, sequences)
#' }
#'
#' ## Sequence masking:
#' if (R.Version()$arch != "i386") {
#' library(Biostrings)
#' data(ArabidopsisMotif)
#' data(ArabidopsisPromoters)
#' seq <- mask_seqs(ArabidopsisPromoters, "AAAAA")
#' scan_sequences(ArabidopsisMotif, seq)
#' # A warning regarding the presence of non-standard letters will be given,
#' # but can be safely ignored in this case.
#' }
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @seealso [add_multifreq()], [Biostrings::matchPWM()],
#'    [enrich_motifs()], [motif_pvalue()]
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.001,
  threshold.type = "pvalue", RC = FALSE, use.freq = 1, verbose = 0,
  nthreads = 1, motif_pvalue.k = 8, use.gaps = TRUE, allow.nonfinite = FALSE,
  warn.NA = TRUE, calc.pvals = FALSE, return.granges = FALSE,
  no.overlaps = FALSE, no.overlaps.by.strand = FALSE, no.overlaps.strat = "score",
  respect.strand = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!threshold.type %in% c("logodds", "pvalue", "logodds.abs")) {
    threshold.type_check <- wmsg2(paste0(" * Incorrect 'threshold.type': expected ",
                                         "`logodds`, `logodds.abs` or `pvalue`; got `",
                                         threshold.type, "`"), exdent = 3, indent = 1)
    all_checks <- c(all_checks, threshold.type_check)
  }
  char_check <- check_fun_params(list(threshold.type = args$threshold.type,
                                      no.overlaps.strat = args$no.overlaps.strat),
                                 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(threshold = args$threshold,
                                     use.freq = args$use.freq,
                                     verbose = args$verbose,
                                     nthreads = args$nthreads,
                                     motif_pvalue.k = args$motif_pvalue.k),
                                c(0, 1, 1, 1, 1), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(RC = args$RC, use.gaps = args$use.gaps,
                                      return.granges = args$return.granges,
                                      no.overlaps = args$no.overlaps,
                                      no.overlaps.by.strand = args$no.overlaps.by.strand),
                                 numeric(), logical(), TYPE_LOGI)
  s4_check <- check_fun_params(list(sequences = args$sequences), numeric(),
                               logical(), TYPE_S4)
  all_checks <- c(all_checks, char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (verbose > 2) {
    message(" * Input parameters")
    message("   * motifs:              ", deparse(substitute(motifs)))
    message("   * sequences:           ", deparse(substitute(sequences)))
    message("   * threshold:           ", ifelse(length(threshold) > 1, "...",
                                                 threshold))
    message("   * threshold.type:      ", threshold.type)
    message("   * RC:                  ", RC)
    message("   * use.freq:            ", use.freq)
    message("   * verbose:             ", verbose)
  }

  if (!no.overlaps.strat %in% c("score", "order"))
    stop("`no.overlaps.strat` must be \"score\" or \"order\"", call. = FALSE)

  if (missing(motifs) || missing(sequences)) {
    stop("need both motifs and sequences")
  }

  if (!RC && respect.strand)
    stop("`respect.strand` cannot be TRUE if `RC = FALSE`")

  if (verbose > 0) message(" * Processing motifs")

  if (verbose > 1) message("   * Scanning ", length(motifs),
                           ifelse(length(motifs) > 1, " motifs", " motif"))

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  motifs <- convert_type_internal(motifs, "PWM")
  needsfix <- vapply(motifs, function(x) any(is.infinite(x@motif)), logical(1))
  if (any(needsfix) && !allow.nonfinite) {
    message(wmsg("Note: found -Inf values in motif PWM(s), adding a pseudocount. ",
      "Set `allow.nonfinite = TRUE` to prevent this behaviour."))
    for (i in which(needsfix)) {
      motifs[[i]] <- suppressMessages(normalize(motifs[[i]]))
    }
  }

  mot.names <- vapply(motifs, function(x) x@name, character(1))

  mot.gaps <- lapply(motifs, function(x) x@gapinfo)
  mot.hasgap <- vapply(mot.gaps, function(x) x@isgapped, logical(1))
  if (any(mot.hasgap) && use.gaps) {
    gapdat <- process_gapped_motifs(motifs, mot.hasgap)
  }

  mot.pwms <- lapply(motifs, function(x) x@motif)
  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) != 1) stop("can only scan using one alphabet")
  mot.alphs <- unique(mot.alphs)
  if (verbose > 1) message("   * Motif alphabet: ", mot.alphs)

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- as.character(seq_len(length(sequences)))

  seq.alph <- seqtype(sequences)
  if (seq.alph != "B" && seq.alph != mot.alphs)
    stop("Motif and Sequence alphabets do not match")
  else if (seq.alph == "B")
    seq.alph <- mot.alphs
  if (RC && !seq.alph %in% c("DNA", "RNA"))
    stop("`RC = TRUE` is only valid for DNA/RNA motifs")

  if (use.freq > 1) {
    if (any(mot.hasgap) && use.gaps)
      stop("use.freq > 1 cannot be used with gapped motifs")
    if (any(vapply(motifs, function(x) length(x@multifreq) == 0, logical(1))))
      stop("missing multifreq slots")
    check_multi <- vapply(motifs,
                          function(x) any(names(x@multifreq) %in%
                                          as.character(use.freq)),
                          logical(1))
    if (!any(check_multi)) stop("not all motifs have correct multifreqs")
  }

  if (use.freq == 1) {
    score.mats <- mot.pwms
  } else {
    score.mats <- lapply(motifs,
                         function(x) x@multifreq[[as.character(use.freq)]])
    for (i in seq_along(score.mats)) {
      score.mats[[i]] <- MATRIX_ppm_to_pwm(score.mats[[i]],
                                           nsites = motifs[[i]]@nsites,
                                           pseudocount = motifs[[i]]@pseudocount,
                                           bkg = motifs[[i]]@bkg[rownames(score.mats[[i]])])
    }
  }

  max.scores <- vapply(motifs, function(x)
    suppressMessages(motif_score(x, 1, use.freq, threshold.type = "fromzero",
        allow.nonfinite = allow.nonfinite)),
    numeric(1))
  if (!allow.nonfinite)
    min.scores <- vapply(motifs, function(x)
      suppressMessages(motif_score(x, 0, use.freq)), numeric(1))
  else
    min.scores <- vapply(motifs, function(x) motif_score_min(x, use.freq), numeric(1))

  switch(threshold.type,

    "logodds" = {

      thresholds <- max.scores * threshold

    },

    "logodds.abs" = {

      if (!length(threshold) %in% c(length(motifs), 1))
        stop(wmsg("for threshold.type = 'logodds.abs', a threshold must be provided for
                  every single motif or one threshold recycled for all motifs"))

      if (length(threshold) == 1) threshold <- rep(threshold, length(motifs))
      thresholds <- threshold

    },

    "pvalue" = {

      if (verbose > 0)
        message(" * Converting P-values to logodds thresholds")
      thresholds <- vector("numeric", length(motifs))
      thresholds <- motif_pvalue(motifs, pvalue = threshold, use.freq = use.freq,
                                 k = motif_pvalue.k, allow.nonfinite = allow.nonfinite)
      for (i in seq_along(thresholds)) {
        if (thresholds[i] > max.scores[i]) thresholds[i] <- max.scores[i]
      }
      if (verbose > 3) {
        for (i in seq_along(thresholds)) {
          message("   * Motif ", mot.names[i], ": max.score = ", max.scores[i],
                  ", threshold = ", thresholds[i])
        }
      }
      thresholds <- unlist(thresholds)

    },

    stop("unknown 'threshold.type'")

  )

  for (i in seq_along(threshold)) {
    if (threshold[i] > max.scores[i])
      warning(wmsg("Threshold [", threshold[i], "] for motif ", i,
          " is higher than the max possible threshold [", max.scores[i], "]"),
        immediate. = TRUE, call. = FALSE)
  }

  alph <- switch(seq.alph, "DNA" = "ACGT", "RNA" = "ACGU",
                 "AA" = collapse_cpp(AA_STANDARD2), seq.alph)
  sequences <- as.character(sequences)
  strands <- rep("+", length(score.mats))

  if (any(mot.hasgap) && use.gaps) {
    strands <- strands[gapdat$IDs]
    mot.names <- mot.names[gapdat$IDs]
    score.mats <- lapply(gapdat$motifs, function(x) x@motif)
    thresholds <- thresholds[gapdat$IDs]
    min.scores <- min.scores[gapdat$IDs]
    max.scores <- max.scores[gapdat$IDs]
  }

  if (RC) {
    if (respect.strand) {
      mot.strands <- vapply(motifs, function(x) x@strand, character(1))
      keep.pos <- rep(TRUE, length(motifs))
      keep.neg <- rep(TRUE, length(motifs))
      keep.pos[mot.strands == "-"] <- FALSE
      keep.neg[mot.strands == "+"] <- FALSE
    } else {
      keep.pos <- rep(TRUE, length(motifs))
      keep.neg <- rep(TRUE, length(motifs))
    }
    strands <- c(strands[keep.pos], rep("-", length(score.mats))[keep.neg])
    mot.names <- c(mot.names[keep.pos], mot.names[keep.neg])
    thresholds <- c(thresholds[keep.pos], thresholds[keep.neg])
    score.mats.rc <- lapply(score.mats,
                            function(x) matrix(rev(as.numeric(x)), nrow = nrow(x)))
    score.mats <- c(score.mats[keep.pos], score.mats.rc[keep.neg])
    min.scores <- c(min.scores[keep.pos], min.scores[keep.neg])
    max.scores <- c(max.scores[keep.pos], max.scores[keep.neg])
    mot.indices <- c(seq_along(motifs)[keep.pos], seq_along(motifs)[keep.neg])
    motifs <- c(motifs[keep.pos], motifs[keep.neg])
  }

  thresholds[thresholds == Inf] <- min_max_ints()$max / 1000
  thresholds[thresholds == -Inf] <- min_max_ints()$min / 1000

  if (allow.nonfinite) {
    for (i in seq_along(score.mats)) {
      if (any(is.infinite(score.mats[[i]]))) {
        min_val1 <- min_max_ints()$min / ncol(score.mats[[i]])
        min_val2 <- as.integer(log2(nrow(score.mats[[i]])) * ncol(score.mats[[i]])) * 1000
        min_val <- (min_val1 + min_val2) / 1000
        score.mats[[i]][is.infinite(score.mats[[i]])] <- min_val
      }
    }
  }

  if (verbose > 0) message(" * Scanning")

  res <- scan_sequences_cpp(score.mats, sequences, use.freq, alph, thresholds,
    nthreads, allow.nonfinite, warn.NA)

  if (verbose > 1) message("   * Number of matches: ", nrow(res))
  if (verbose > 0) message(" * Processing results")

  thresholds[thresholds <= min_max_ints()$min / 1000] <- -Inf
  thresholds[thresholds >= min_max_ints()$max / 1000] <- Inf

  res$thresh.score <- thresholds[res$motif]
  res$min.score <- min.scores[res$motif]
  res$max.score <- max.scores[res$motif]
  res$score.pct <- res$score / res$max.score * 100
  if (seq.alph %in% c("DNA", "RNA")) res$strand <- strands[res$motif]
  res$motif <- mot.names[res$motif]
  res$sequence <- seq.names[res$sequence]

  if (nrow(res) == 0) message("No hits found.")

  if (RC && nrow(res) > 0) res <- adjust_rc_hits(res, seq.alph)

  if (RC) res$motif.i <- mot.indices[res$motif.i]

  out <- as(res, "DataFrame")
  out@metadata <- list(
    args = args[-c(1:2)],
    seqlengths = structure(width(sequences), names = names(sequences))
  )

  if (any(mot.hasgap) && use.gaps) {
    out$match <- add_gap_dots_cpp(out$match, gapdat$gaplocs[out$motif.i])
    out$motif.i <- gapdat$IDs[out$motif.i]
  }

  if (calc.pvals) {
    out$pvalue <- motif_pvalue(motifs[out$motif.i], out$score, use.freq = use.freq,
      nthreads = nthreads, allow.nonfinite = allow.nonfinite, k = motif_pvalue.k)
  }

  if (no.overlaps && nrow(out) > 1) {
    # TODO: The no.overlaps code is rather slow, not too happy.
    if (RC && no.overlaps.by.strand) {
      row.indices.plus <- which(out$strand == "+")
      row.indices.minus <- which(out$strand == "-")
      row.indices.plus <- remove_masked_hits(out, row.indices.plus, no.overlaps.strat)
      row.indices.minus <- remove_masked_hits(switch_antisense_coords_cpp(out),
        row.indices.minus, no.overlaps.strat)
      row.indices <- c(row.indices.plus, row.indices.minus)
    } else if (RC) {
      row.indices <- remove_masked_hits(switch_antisense_coords_cpp(out),
        seq_len(nrow(out)), no.overlaps.strat)
    } else {
      row.indices <- seq_len(nrow(out))
      row.indices <- remove_masked_hits(out, seq_len(nrow(out)), no.overlaps.strat)
    }
    out <- out[row.indices, ]
  }

  if (return.granges) {
    colnames(out)[3] <- "seqname"
    if (RC) {
      out <- switch_antisense_coords_cpp(out)
    }
    out <- granges_fun(GenomicRanges::GRanges(out,
        seqlengths = structure(width(sequences), names = names(sequences))))
  }

  out

}

# What about this kind of situation?
#  seq:   CAAAAACCAAAACCAAAACC
#  hit 1: ++++++++
#  hit 2:       ++++++++
#  hit 3:             ++++++++
# What should happen here? Right now only one of the three is kept.

remove_masked_hits <- function(x, i = seq_len(nrow(x)), strat = "score") {
  if (!length(i)) return(i)
  y <- x[i, ]
  y$index.tokeep <- i
  switch(strat, score = remove_masked_hits_by_score(y),
    order = remove_masked_hits_by_order(y))
}

remove_masked_hits_by_order <- function(y) {
  sort(do.call(c, by(y, list(y$sequence, y$motif.i), function(z) {
    dedup_by_order(z, flatten_group_matrix(get_overlap_groups(z)))
  }, simplify = FALSE)))
}

remove_masked_hits_by_score <- function(y) {
  sort(do.call(c, by(y, list(y$sequence, y$motif.i), function(z) {
    dedup_by_score(z, flatten_group_matrix(get_overlap_groups(z)))
  }, simplify = FALSE)))
}

get_overlap_groups <- function(x) {
  y <- as.matrix(findOverlaps(IRanges(x$start, x$stop)))
  y <- xtabs(~queryHits + subjectHits, y)
  matrix(as.integer(y), nrow = nrow(x))
}

flatten_group_matrix <- function(x) {
  if (all(x == 1)) return(seq_len(nrow(x)))
  cutree(hclust(as.dist(1 - x)), h = 0.5)
}

dedup_by_order <- function(x, i) {
  x$index.tokeep[!duplicated(i)]
}

dedup_by_score <- function(x, i) {
  do.call(c, by(x, i, function(y) {
    y$index.tokeep[which.max(y$score)]
  }, simplify = FALSE))
}

adjust_rc_hits <- function(res, alph) {
  rev.strand <- res$strand == "-"
  if (any(rev.strand)) {
    start <- res$stop[rev.strand]
    stop <- res$start[rev.strand]
    res$stop[rev.strand] <- stop
    res$start[rev.strand] <- start
    matches <- res$match[rev.strand]
    if (alph == "DNA")
      matches <- as.character(reverseComplement(DNAStringSet(matches)))
    else if (alph == "RNA")
      matches <- as.character(reverseComplement(RNAStringSet(matches)))
    res$match[rev.strand] <- matches
  }
  res
}

# Note: It's probably a lot faster to scan the individual submotifs and then
# process the gapped motifs afterwards, versus scanning all possible gapped
# motif combinations. Would need to think about how to score the submotifs
# though; so for now, go with the dumb and slow brute force option.

process_gapped_motifs <- function(motifs, hasgap) {
  motifs[hasgap] <- lapply(motifs[hasgap], ungap_single)
  motifs_gapped <- mapply(function(x, y) rep(x, length(y)), hasgap, motifs,
    SIMPLIFY = FALSE)
  IDs <- mapply(function(x, y) rep(x, length(y)), seq_along(motifs), motifs,
    SIMPLIFY = FALSE)
  out <- list(
    motifs = do.call(c, motifs),
    gapped = do.call(c, motifs_gapped),
    IDs = do.call(c, IDs)
  )
  out$gaplocs <- lapply(seq_along(out$motifs), function(x) integer())
  out$gaplocs[out$gapped] <- lapply(out$motifs[out$gapped], get_gaplocs)
  out
}

get_gaplocs <- function(x) {
  y <- strsplit(x@name, "/", fixed = TRUE)[[1]]
  npos <- seq_len(ncol(x@motif))
  lens <- vapply(y, function(x) strsplit(x, "_L", fixed = TRUE)[[1]][2], character(1))
  lens <- as.numeric(lens)
  lens <- lapply(lens, function(x) seq(1, x))
  lenslens <- cumsum(vapply(lens[-length(lens)], length, integer(1)))
  lens <- mapply(function(x, y) x + y, lens, c(0, lenslens), SIMPLIFY = FALSE)
  gapped <- grepl("BLANK", y)
  do.call(c, lens[gapped])
}

get_submotifs <- function(m) {
  n <- length(m@gapinfo@gaploc)
  mname <- m@name
  submotifs <- vector("list", n + 1)
  submotifs[[1]] <- subset(m, seq(1, m@gapinfo@gaploc[1]))
  submotifs[[length(submotifs)]] <- subset(
    m, seq(m@gapinfo@gaploc[n] + 1, ncol(m))
  )
  if (length(submotifs) > 2) {
    for (i in seq_along(submotifs)[-c(1, length(submotifs))]) {
      submotifs[[i]] <- subset(
        m, seq(m@gapinfo@gaploc[i - 1] + 1, m@gapinfo@gaploc[i])
      )
    }
  }
  for (i in seq_along(submotifs)) {
    submotifs[[i]]@name <- paste0("SUB_N", i, "_L", ncol(submotifs[[i]]@motif))
  }
  submotifs
}

make_blank_motif <- function(n, N, alph) {
  alphlen <- switch(alph, DNA = 4, RNA = 4, AA = 20, nchar(alph))
  mot <- matrix(0, nrow = alphlen, ncol = n)
  create_motif(mot, type = "PWM", alphabet = alph,
    name = paste0("BLANK_N", N, "_L", n))
}

ungap_single <- function(m) {
  gaplens <- mapply(
    seq, m@gapinfo@mingap, m@gapinfo@maxgap, SIMPLIFY = FALSE
  )
  gaplens <- expand.grid(gaplens)
  out <- vector("list", nrow(gaplens))
  submotifs <- get_submotifs(m)
  for (i in seq_along(out)) {
    tmp <- list(submotifs[[1]])
    for (j in seq_len(ncol(gaplens))) {
      if (gaplens[[j]][i] == 0) {
        tmp <- c(tmp, list(submotifs[[j + 1]]))
      } else {
        tmp <- c(tmp, list(make_blank_motif(gaplens[[j]][i], j, m@alphabet),
            submotifs[[j + 1]]))
      }
    }
    out[[i]] <- do.call(cbind, tmp)
  }
  out
}

motif_score_min <- function(x, use.freq) {
  if (any(is.infinite(x@motif)))
    -Inf
  else
    suppressMessages(motif_score(x, 0, use.freq))
}

granges_fun <- function(FUN, env = parent.frame()) {
  if (requireNamespace("GenomicRanges", quietly = TRUE)) {
    eval(substitute(FUN), envir = env)
  } else {
    stop(wmsg("The 'GenomicRanges' package must be installed for `return.granges=TRUE`. ",
        "[BiocManager::install(\"GenomicRanges\")]"), call. = FALSE)
  }
}
