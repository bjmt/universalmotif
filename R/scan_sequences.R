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
#'    If `threshold.type = 'logodds'`, then to calculate the minimum
#'    allowed score the max possible score for a motif is multiplied
#'    by the value set by `threshold`. To determine the maximum 
#'    possible scores a motif (of type PWM), run
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
#' seq <- ArabidopsisPromoters[[1]]  # Only works for XString, not XStringSet
#' seq <- mask(seq, pattern = "AAAA")  # MaskedDNAString class
#' seq <- injectHardMask(seq, letter = "+")  # Recover XString
#' seq <- DNAStringSet(seq)  # scan_sequences() needs XStringSet
#' scan_sequences(ArabidopsisMotif, seq)
#' # A warning regarding the presence of non-standard letters will be given,
#' # but can be safely ignored in this case.
#'
#' # If you'd like to do it to a whole DNAStringSet object:
#' seq <- ArabidopsisPromoters
#' for (i in 1:length(ArabidopsisPromoters)) {
#'   seq[[i]] <- injectHardMask(mask(seq[[i]], pattern = "AAAA"), letter = "+")
#' }
#' scan_sequences(ArabidopsisMotif, seq)
#' }
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [add_multifreq()], [Biostrings::matchPWM()],
#'    [enrich_motifs()], [motif_pvalue()]
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.001,
                           threshold.type = "pvalue", RC = FALSE,
                           use.freq = 1, verbose = 0, nthreads = 1,
                           motif_pvalue.k = 8, use.gaps = TRUE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!threshold.type %in% c("logodds", "pvalue", "logodds.abs")) {
    threshold.type_check <- wmsg2(paste0(" * Incorrect 'threshold.type': expected ",
                                         "`logodds`, `logodds.abs` or `pvalue`; got `",
                                         threshold.type, "`"), exdent = 3, indent = 1)
    all_checks <- c(all_checks, threshold.type_check)
  }
  char_check <- check_fun_params(list(threshold.type = args$threshold.type),
                                 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(threshold = args$threshold,
                                     use.freq = args$use.freq,
                                     verbose = args$verbose,
                                     nthreads = args$nthreads,
                                     motif_pvalue.k = args$motif_pvalue.k),
                                c(0, 1, 1, 1, 1), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(RC = args$RC, use.gaps = args$use.gaps),
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

  if (missing(motifs) || missing(sequences)) {
    stop("need both motifs and sequences")
  }

  if (verbose > 0) message(" * Processing motifs")

  if (verbose > 1) message("   * Scanning ", length(motifs),
                           ifelse(length(motifs) > 1, " motifs", " motif"))

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  motifs <- convert_type_internal(motifs, "PWM")
  motifs <- lapply(motifs, fix_mots)

  mot.names <- vapply(motifs, function(x) x@name, character(1))

  mot.gaps <- lapply(motifs, function(x) x@gapinfo)
  mot.hasgap <- vapply(mot.gaps, function(x) x@isgapped, logical(1))
  if (any(mot.hasgap) && use.gaps) {
    gapdat <- process_gapped_motifs(motifs, mot.hasgap)
    # motifs.ungapped <- gapdat$motifs
    # mot.names.ungapped <- mot.names[gapdat$IDs]
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

  max.scores <- vapply(motifs, function(x) motif_score(x, 1, use.freq), numeric(1))
  min.scores <- vapply(motifs, function(x) motif_score(x, 0, use.freq), numeric(1))

  switch(threshold.type,

    "logodds" = {

      # thresholds <- ((abs(max.scores) + abs(min.scores)) * threshold) - abs(min.scores)
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
                                 k = motif_pvalue.k)
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
    strands <- c(strands, rep("-", length(score.mats)))
    mot.names <- c(mot.names, mot.names)
    thresholds <- c(thresholds, thresholds)
    score.mats.rc <- lapply(score.mats,
                            function(x) matrix(rev(as.numeric(x)), nrow = nrow(x)))
    score.mats <- c(score.mats, score.mats.rc)
    min.scores <- c(min.scores, min.scores)
    max.scores <- c(max.scores, max.scores)
  }

  if (verbose > 0) message(" * Scanning")

  res <- scan_sequences_cpp(score.mats, sequences, use.freq, alph, thresholds, nthreads)

  if (verbose > 1) message("   * Number of matches: ", nrow(res))
  if (verbose > 0) message(" * Processing results")

  res$thresh.score <- thresholds[res$motif]
  res$min.score <- min.scores[res$motif]
  res$max.score <- max.scores[res$motif]
  res$score.pct <- res$score / res$max.score * 100
  if (seq.alph %in% c("DNA", "RNA")) res$strand <- strands[res$motif]
  res$motif <- mot.names[res$motif]
  res$sequence <- seq.names[res$sequence]

  if (nrow(res) == 0) message("No hits found.")

  if (RC && nrow(res) > 0) res <- adjust_rc_hits(res, seq.alph)

  out <- as(res, "DataFrame")
  out@metadata <- args[-c(1:2)]

  if (any(mot.hasgap) && use.gaps) {
    out$match <- add_gap_dots_cpp(out$match, gapdat$gaplocs[out$motif.i])
    out$motif.i <- gapdat$IDs[out$motif.i]
  }

  out

}

fix_mots <- function(x) {
  if (any(is.infinite(x@motif))) {
    message("Note: detected non-finite positional weights, normalising")
    normalize(x)
  } else {
    x
  }
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
