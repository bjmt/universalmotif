#' Scan sequences for matches to input motifs.
#'
#' For sequences of any alphabet, scan them using the PWM matrices of
#' a set of input motifs.
#'
#' @param motifs See `convert_motifs()` for acceptable motif formats.
#' @param sequences \code{\link{XStringSet}} Sequences to scan. Alphabet
#'    should match motif.
#' @param threshold `numeric(1)` Between 0 and 1. See details.
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
#' @param progress `logical(1)` Deprecated. Does nothing.
#' @param BP `logical(1)` Deprecated. See `ncores`.
#' @param nthreads `numeric(1)` Run [scan_sequences()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'    The work is split by motif, so no speed up will occur for jobs with a
#'    single motif.
#'
#' @return `data.frame` with each row representing one hit; if the input
#'    sequences are \code{\link{DNAStringSet}} or
#'    \code{\link{RNAStringSet}}, then an
#'    additional column with the strand is included.
#'
#' @details
#'    Similar to [Biostrings::matchPWM()], the scanning method uses
#'    logodds scoring. (To see the scoring matrix for any motif, simply
#'    run `convert_type(motif, "PWM")`; for a `multifreq` scoring
#'    matrix: `apply(motif["multifreq"]$`2`, 2, ppm_to_pwm)`). In order
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
#'    Finaly, if \code{threshold.type = 'logodds.abs'}, then the exact values
#'    provided will be used as thresholds.
#'
#'    Non-standard letters (such as "N", "+", "-", ".", etc in `DNAString`
#'    objects) will be safely ignored, resulting only in a warning and a very
#'    minor performance cost. This can used to scan
#'    masked sequences. See \code{\link[Biostrings:maskMotif]{Biostrings::mask()}}
#'    for masking sequences
#'    (generating `MaskedXString` objects), and [Biostrings::injectHardMask()]
#'    to recover masked `XStringSet` objects for use with [scan_sequences()].
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
#' library(Biostrings)
#' data(ArabidopsisMotif)
#' data(ArabidopsisPromoters)
#' seq <- ArabidopsisPromoters[[1]]  # Only works for XString, not XStringSet
#' seq <- mask(seq, pattern = "AAAA")  # MaskedDNAString class
#' seq <- injectHardMask(seq, letter = "+")  # Recover XString
#' seq <- DNAStringSet(seq)  # scan_sequences() needs XStringSet
#' scan_sequences(ArabidopsisMotif, seq, verbose = 0, progress = FALSE)
#' # A warning regarding the presence of non-standard letters will be given,
#' # but can be safely ignored in this case.
#'
#' @references
#'    \insertRef{biostrings}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [add_multifreq()], [Biostrings::matchPWM()],
#'    [enrich_motifs()], [motif_pvalue()]
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.001,
                           threshold.type = "pvalue", RC = FALSE,
                           use.freq = 1, verbose = 0, progress = FALSE,
                           BP = FALSE, nthreads = 1) {

  # TODO: Work with Masked*String objects. Masked letters show up as "#" after
  #       as.character() calls, which should just cause scan_sequences() to
  #       ignore these and work as intended. For now, just using "-", "." or
  #       "+" within DNA/RNA/AAStringSet objects will work the same.

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
                                     nthreads = args$nthreads),
                                c(0, 1, 1, 1), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(RC = args$RC),
                                 numeric(), logical(), TYPE_LOGI)
  s4_check <- check_fun_params(list(sequences = args$sequences), numeric(),
                               logical(), TYPE_S4)
  all_checks <- c(all_checks, char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (progress)
    warning("'progress' is deprecated and does nothing")
  if (BP)
    warning("'BP' is deprecated; use 'ncores' instead")

  if (verbose <= 0) progress <- FALSE

  if (verbose > 2) {
    cat(" * Input parameters\n")
    cat("   * motifs:              ", deparse(substitute(motifs)), "\n")
    cat("   * sequences:           ", deparse(substitute(sequences)), "\n")
    cat("   * threshold:           ", ifelse(length(threshold) > 1, "...",
                                             threshold), "\n")
    cat("   * threshold.type:      ", threshold.type, "\n")
    cat("   * RC:                  ", RC, "\n")
    cat("   * use.freq:            ", use.freq, "\n")
    cat("   * verbose:             ", verbose, "\n")
  }

  if (missing(motifs) || missing(sequences)) {
    stop("need both motifs and sequences")
  }

  if (verbose > 0) cat(" * Processing motifs\n")

  if (verbose > 1) cat("   * Scanning", length(motifs),
                       ifelse(length(motifs) > 1, "motifs\n", "motif\n"))

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  motifs <- convert_type_internal(motifs, "PWM")
  motifs <- lapply(motifs, function(x) if (any(is.infinite(x@motif)))
                                         normalize(x) else x)

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  mot.pwms <- lapply(motifs, function(x) x@motif)
  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) != 1) stop("can only scan using one alphabet")
  mot.alphs <- unique(mot.alphs)
  if (verbose > 1) cat("   * Motif alphabet:", mot.alphs, "\n")

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))

  seq.alph <- seqtype(sequences)
  if (seq.alph != "B" && seq.alph != mot.alphs)
    stop("Motif and Sequence alphabets do not match")
  else if (seq.alph == "B")
    seq.alph <- mot.alphs

  if (use.freq > 1) {
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
      score.mats[[i]] <- apply(score.mats[[i]], 2, ppm_to_pwmC,
                               nsites = motifs[[i]]@nsites,
                               pseudocount = motifs[[i]]@pseudocount)
    }
  }

  max.scores <- vapply(motifs, function(x) motif_score(x, 1), numeric(1))
  min.scores <- vapply(motifs, function(x) motif_score(x, 0), numeric(1))

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

      if (progress && !BP && verbose > 0)
        cat(" * Converting P-values to logodds thresholds ...")
      else if ((progress && BP && verbose > 0) || verbose > 0)
        cat(" * Converting P-values to logodds thresholds\n")
      thresholds <- vector("numeric", length(motifs))
      thresholds <- motif_pvalue(motifs, pvalue = threshold, use.freq = use.freq,
                                 k = 6, progress = progress, BP = BP)
      for (i in seq_along(thresholds)) {
        if (thresholds[i] > max.scores[i]) thresholds[i] <- max.scores[i]
      }
      if (verbose > 3) {
        for (i in seq_along(thresholds)) {
          cat("   * Motif ", mot.names[i], ": max.score = ", max.scores[i],
              ", threshold = ", thresholds[i], "\n", sep = "")
        }
      }
      thresholds <- unlist(thresholds)

    },

    stop("unknown 'threshold.type'")

  )

  alph <- switch(seq.alph, "DNA" = "ACGT", "RNA" = "ACGU",
                 "AA" = collapse_cpp(AA_STANDARD), seq.alph)
  sequences <- as.character(sequences)
  strands <- rep("+", length(score.mats))

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

  if (verbose > 0) cat(" * Scanning\n")

  res <- scan_sequences_cpp(score.mats, sequences, use.freq, alph, thresholds, nthreads)

  if (verbose > 1) cat ("   * Number of matches: ", nrow(res), "\n", sep = "")
  if (verbose > 0) cat(" * Processing results\n")

  res$thresh.score <- thresholds[res$motif]
  res$min.score <- min.scores[res$motif]
  res$max.score <- max.scores[res$motif]
  res$score.pct <- res$score / res$max.score * 100
  if (seq.alph %in% c("DNA", "RNA")) res$strand <- strands[res$motif]
  res$motif <- mot.names[res$motif]
  res$sequence <- seq.names[res$sequence]

  if (nrow(res) == 0) message("No hits found.")

  if (RC && nrow(res) > 0) res <- adjust_rc_hits(res, seq.alph)

  res[, c(1:5, 6, 11, 7:10)]

}

adjust_rc_hits <- function(res, alph) {
  rev.strand <- res$strand == "-"
  if (any(rev.strand)) {
    start <- res$stop[rev.strand]
    stop <- res$start[rev.strand]
    res$stop <- stop
    res$start <- start
    matches <- res$match[rev.strand]
    if (alph == "DNA")
      matches <- as.character(reverseComplement(DNAStringSet(matches)))
    else if (alph == "RNA")
      matches <- as.character(reverseComplement(RNAStringSet(matches)))
    res$match[rev.strand] <- matches
  }
  res
}
