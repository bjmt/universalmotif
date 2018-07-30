#' Find motif binding sites in a set of sequences.
#'
#' Find matches to motifs in a set of input sequences.
#'
#' @param motifs \linkS4class{universalmotif} objects.
#' @param sequences XStringSet object. Sequences to scan. Alphabet should
#'    match motif.
#' @param threshold Numeric. Between 0 and 1. See details.
#' @param threshold.type Character. One of 'logodds' and 'pvalue'. See details.
#' @param RC Logical. If \code{TRUE}, check reverse complement of input
#'    sequences.
#' @param use.freq Numeric. The default, 1, uses the motif matrix (from
#'    the \code{motif["motif"]} slot) to search for sequences. If a higher
#'    number is used, then the matching k-let matrix from the
#'    \code{motif["multifreq"]} slot is used. See \code{\link{add_multifreq}}.
#' @param progress_bar Logical. Shoe progress bar.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Results as a data.frame object.
#'
#' @details
#'    If \code{use.freq = 1} and the sequences are DNAStringSet or RNAStringSet
#'    class, then the \code{\link[Biostrings]{matchPWM}} function from the
#'    Biostrings package is used. Otherwise, a different (and less efficient)
#'    scanning function is used.
#'    
#'    Similar to \code{\link[Biostrings]{matchPWM}}, the scanning method uses
#'    logodds scoring. (To see the scoring matrix for any motif, simply
#'    run \code{convert_type(motif, "PWM")}; for a \code{multifreq} scoring
#'    matrix: \code{apply(motif["multifreq"]$`2`, 2, ppm_to_pwm)}). In order
#'    to score a sequence, at each position within a sequence of length equal
#'    to the length of the motif, the scores for each base are summed. If the
#'    score sum is above the desired threshold, it is kept.
#'
#'    If \code{threshold.type = 'logodds'}, then to calculate the minimum
#'    allowed score the maximum possible score for a motif is multiplied
#'    by the value set by \code{threshold}. To determine the maximum
#'    possible score a motif (of type PWM), run
#'    \code{sum(apply(motif['motif'], 2, max))}.
#'
#'    Note: memory usage increases exponentially with k.
#'
#' @examples
#' # any alphabet can be used
#' \dontrun{
#' set.seed(1)
#' alphabet <- paste(c(letters), collapse = "")
#' motif <- create_motif("hello", alphabet = alphabet)
#' sequences <- create_sequences(alphabet, seqnum = 1000, seqlen = 100000)
#' scan_sequences(motif, sequences)
#' }
#'
#' @references
#'    \insertRef{biostrings}{universalmotif}
#'
#'    \insertRef{tfmpvalue}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso \code{\link{add_multifreq}}, \code{\link[Biostrings]{matchPWM}},
#'    \code{\link{enrich_motifs}}
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.01,
                            threshold.type = "pvalue", RC = FALSE,
                            use.freq = 1, verbose = TRUE,
                            progress_bar = FALSE, BPPARAM = SerialParam()) {

  if (missing(motifs) || missing(sequences)) {
    stop("need both motifs and sequences")
  }

  if (progress_bar) pb_prev <- BPPARAM$progressbar

  if (verbose) cat(" * Processing motifs\n")

  if (!is.list(motifs)) motifs <- list(motifs)

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)

  mot.names <- vapply(motifs, function(x) x["name"], character(1))
  mot.pwms <- convert_type(motifs, "PWM", BPPARAM = BPPARAM)
  mot.pwms <- lapply(mot.pwms, function(x) x["motif"])
  mot.lens <- vapply(mot.pwms, ncol, numeric(1))
  mot.alphs <- vapply(motifs, function(x) x["alphabet"], character(1))
  if (length(unique(mot.alphs)) != 1) stop("can only scan using one alphabet")
  mot.alphs <- unique(mot.alphs)
  alph <- unique(mot.alphs)
  mot.bkgs <- lapply(motifs, function(x) x["bkg"])
  seq.lens <- width(sequences)

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))

  if (use.freq > 1) {
    if (any(vapply(motifs, function(x) length(x["multifreq"]) == 0, logical(1))))
      stop("missing multifreq slots")
    check_multi <- vapply(motifs,
                          function(x) any(names(x["multifreq"]) %in%
                                          as.character(use.freq)),
                          logical(1))
    if (!any(check_multi)) stop("not all motifs have correct multifreqs")
  }

  for (i in seq_along(motifs)) {
    if (motifs[[i]]["pseudocount"] == 0) motifs[[i]]["pseudocount"] <- 0.0001
    if (length(motifs[[i]]["nsites"]) == 0) motifs[[i]]["nsites"] <- 100
  }
  
  if (use.freq == 1) {
    score.mats <- mot.pwms
  } else {
    score.mats <- lapply(motifs,
                         function(x) x["multifreq"][[as.character(use.freq)]]) 
    for (i in seq_along(score.mats)) {
      score.mats[[i]] <- apply(score.mats[[i]], 2, ppm_to_pwmC,
                               nsites = motifs[[i]]["nsites"],
                               pseudocount = motifs[[i]]["pseudocount"])
    }
  }

  max.scores <- vapply(score.mats, function(x) sum(apply(x, 2, max)), numeric(1))
  min.scores <- vapply(score.mats, function(x) sum(apply(x, 2, min)), numeric(1))

  if (threshold.type == "logodds") {
    thresholds <- max.scores * threshold
  } else if (threshold.type == "pvalue") {
    if (verbose) cat(" * Converting P-values to logodds thresholds\n")
    thresholds <- vector("numeric", length(motifs))
    if (progress_bar) BPPARAM$progressbar <- TRUE
    thresholds <- motif_pvalue(motifs, pvalue = threshold, use.freq = use.freq,
                               BPPARAM = BPPARAM)
    if (progress_bar) BPPARAM$progressbar <- FALSE
    thresholds <- unlist(thresholds)
  } else stop("unknown 'threshold.type'")

  if (verbose) cat(" * Processing sequences\n")

  seqs.aschar <- as.character(sequences)
  seqs.aschar <- bplapply(seqs.aschar, function(x) strsplit(x, "")[[1]],
                          BPPARAM = BPPARAM)

  seq.lens <- width(sequences)
  seq.matrices <- lapply(seq.lens, function(x) matrix(ncol = x - use.freq + 1,
                                                      nrow = use.freq))

  seq.matrices <- bpmapply(.process_seqs, seq.matrices, seqs.aschar,
                           MoreArgs = list(k = use.freq),
                           BPPARAM = BPPARAM, SIMPLIFY = FALSE)

  if (mot.alphs == "DNA") {
    mot.alphs <- DNA_BASES
  } else if (mot.alphs == "RNA") {
    mot.alphs <- RNA_BASES
  } else if (mot.alphs == "AA") {
    mot.alphs <- AA_STANDARD
  } else if (alph == "custom") {
    if (RC) stop("RC search is only available for DNA/RNA")
    mot.alphs <- lapply(seqs.aschar, unique)
    mot.alphs <- unique(do.call(c, mot.alphs))
  } else {
    if (RC) stop("RC search is only available for DNA/RNA")
    mot.alphs <- strsplit(mot.alphs, "")[[1]]
  }

  alph.int <- as.integer(seq_len(length(mot.alphs)))

  seq.matrices <- bplapply(seq.matrices,
                           function(x) string_to_factor(x, mot.alphs),
                           BPPARAM = BPPARAM)

  seq.ints <- bplapply(seq.matrices,
                       function(x) LETTER_to_int(as.integer(x) - 1,
                                                 use.freq, alph.int),
                       BPPARAM = BPPARAM)

  if (verbose) cat(" * Scanning sequences for motifs\n")

  if (progress_bar) BPPARAM$progressbar <- TRUE

  if (RC && verbose) cat("   * Forward strand\n")
  to.keep <- bplapply(seq_along(score.mats),
                      function(x) .score_motif(seq.ints, score.mats[[x]],
                                               thresholds[x]),
                      BPPARAM = BPPARAM)

  if (RC) {
    score.mats.rc <- lapply(score.mats,
                            function(x) matrix(rev(as.numeric(x)),
                                               ncol = ncol(x)))
    if (verbose) cat("   * Reverse strand\n")
    to.keep.rc <- bplapply(seq_along(score.mats.rc),
                           function(x) .score_motif(seq.ints, score.mats.rc[[x]],
                                                    thresholds[x]),
                           BPPARAM = BPPARAM)
  }

  if (progress_bar) BPPARAM$progressbar <- FALSE

  if (verbose) cat(" * Processing results\n")

  if (RC && verbose) cat("   * Forward strand\n")

  if (progress_bar) pb <- txtProgressBar(max = length(to.keep), style = 3)
  res <- vector("list", length(to.keep))
  for (i in seq_along(res)) {
    to.keep[[i]] <- lapply(to.keep[[i]], res_to_index)
    res[[i]] <- get_res_cpp(to.keep[[i]], seqs.aschar, seq.ints, mot.lens[i],
                         min.scores[i], max.scores[i], mot.names[i],
                         seq.names, score.mats[[i]], "+",
                         seq.lens, use.freq)
    if (progress_bar) setTxtProgressBar(pb, i)
  }
  if (progress_bar) close(pb)

  if (RC) {
    if (verbose) cat("   * Reverse strand\n")

    if (progress_bar) pb <- txtProgressBar(max = length(to.keep), style = 3)
    res.rc <- vector("list", length(to.keep.rc))
    for (i in seq_along(res.rc)) {
      to.keep.rc[[i]] <- lapply(to.keep.rc[[i]], res_to_index)
      res.rc[[i]] <- get_res_cpp(to.keep.rc[[i]], seqs.aschar, seq.ints, mot.lens[i],
                           min.scores[i], max.scores[i], mot.names[i],
                           seq.names, score.mats.rc[[i]], "-",
                           seq.lens, use.freq)
      if (progress_bar) setTxtProgressBar(pb, i)
    }
    if (progress_bar) close(pb)

    res <- c(res, res.rc)
  }

  if (verbose) cat(" * Generating output\n")
  res <- res_list_to_df_cpp(res)

  if (nrow(res) == 0) {
    message("no matches found")
    return(invisible(NULL))
  }

  if (!alph %in% c("DNA", "RNA")) res <- res[, -9]

  if (progress_bar) BPPARAM$progressbar <- pb_prev

  rownames(res) <- NULL

  res

}

.score_motif <- function(seqs, score.1, thresh) {
  lapply(seqs, function(x) scan_seq_internal(x, score.1, thresh))
}

.process_seqs <- function(seq.matrix, seq.aschar, k) {
  for (i in seq_len(k)) {
    to.remove <- k - 1
    if (to.remove == 0) {
      seq.matrix[i, ] <- seq.aschar[seq_len(ncol(seq.matrix))]
      next
    }
    seq.i <- seq.aschar[-seq_len(to.remove)]
    seq.matrix[i, ] <- seq.aschar[seq_len(ncol(seq.matrix))]
  }
  seq.matrix
}
