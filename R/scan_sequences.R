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
#'    If \code{threshold.type = 'pvalue'}, then the P-value as set by
#'    \code{threshold} is converted to a minimum logodds score using
#'    the \code{\link[TFMPvalue]{TFMpv2sc}} function. Note that this
#'    is available for DNA motifs only. Accordingly, if DNA motifs are
#'    used then the P-value for each result is also reported using the
#'    \code{\link[TFMPvalue]{TFMsc2pv}} function.
#'
#'    Note: the memory and processing costs for this function increase
#'    exponentially with increasing k.
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
scan_sequences <- function(motifs, sequences, threshold = 0.6,
                           threshold.type = "logodds", RC = FALSE,
                           use.freq = 1, progress_bar = FALSE,
                           BPPARAM = SerialParam()) {

  if (is.list(motifs)) {
    results <- lapply(motifs, scan_sequences, threshold = threshold,
                      threshold.type = threshold.type,
                      sequences = sequences, RC = RC, 
                      use.freq = use.freq, BPPARAM = BPPARAM)
    results <- do.call(rbind, results)
    rownames(results) <- seq_len(nrow(results))
    return(results)
  }

  if (missing(motifs) || missing(sequences)) {
    stop("missing 'motifs' and/or 'sequences'")
  }

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  mot.name <- motifs["name"]
  mot.mat <- convert_type(motifs, "PWM")["motif"]
  mot.pfm <- convert_type(motifs, "PPM")["motif"]
  mot.len <- ncol(mot.mat)
  max.score <- sum(apply(mot.mat, 2, max))

  bkg <- motifs["bkg"]
  names(bkg) <- DNA_BASES

  pb_prev <- BPPARAM$progressbar

  if (threshold.type == "pvalue" && motifs["alphabet"] == "DNA") {
    threshold <- TFMpv2sc(mot.pfm, threshold, bkg, type = "PFM")
    threshold <- threshold / sum(apply(mot.mat, 2, max))
  } else if (threshold.type == "pvalue") {
    stop("'threshold.type = pvalue' is only valid for DNA")
  }
  if (threshold < 0) stop("cannot have negative threshold")

  min.score <- threshold * 100
  min.score <- paste0(min.score, "%")

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))


  if (use.freq == 1 &&
      motifs["alphabet"] %in% c("DNA", "RNA") &&
      sequences@elementType %in% c("DNAString", "RNAString")) {

    if (motifs["alphabet"] == "RNA") motifs <- switch_alph(motifs)
    if (sequences@elementType == "RNAString") {
      sequences <- DNAStringSet(sequences)
      RNA <- TRUE
    } else RNA <- FALSE

    motif <- motifs
    motif <- convert_type(motif, "PWM")
    motif.rc <- suppressWarnings(motif_rc(motif))
    motif <- motif["motif"]
    motif.rc <- motif.rc["motif"]

    parse_hits <- function(x, y, strand = "+") {
      if (length(x) == 0) {
        data.frame(motif = NULL, sequence = NULL, start = NULL, stop = NULL,
                   max.score = NULL, score.pct = NULL, match = NULL,
                   strand = NULL)
      } else {
        motif <- rep(mot.name, length(x))
        sequence = rep(y, length(x))
        start <- x@ranges@start
        stop <- start + (mot.len - 1)
        if (strand == "-") {
          tmp <- start
          start <- stop
          stop <- tmp
        }
        score <- x@elementMetadata$score 
        max.score <- rep(max.score, length(x))
        score.pct <- (score / max.score) * 100
        match <- as.character(x)

        data.frame(motif = motif, sequence = sequence, start = start, stop = stop,
                   score = score, max.score = max.score, score.pct = score.pct,
                   match = match, strand = rep(strand, length(x)))
      }
    }

    if (progress_bar) BPPARAM$progressbar <- TRUE
    if (progress_bar) cat("Foward strand:\n")
    sequence.hits <- bplapply(seq_len(length(sequences)),
                              function(x) matchPWM(motif, sequences[[x]],
                                                   min.score = min.score,
                                                   with.score = TRUE),
                              BPPARAM = BPPARAM)
    if (progress_bar) BPPARAM$progressbar <- FALSE

    sequence.hits <- bpmapply(parse_hits, sequence.hits, seq.names,
                              BPPARAM = BPPARAM, SIMPLIFY = FALSE)

    sequence.hits <- do.call(rbind, sequence.hits)
    if (nrow(sequence.hits) > 0) sequence.hits$strand <- "+"

    if (RC) {

      if (progress_bar) BPPARAM$progressbar <- TRUE
      if (progress_bar) cat("Reverse strand:\n")
      sequence.hits.rc <- bplapply(seq_len(length(sequences)),
                                   function(x) matchPWM(motif.rc, sequences[[x]],
                                                        min.score = min.score,
                                                        with.score = TRUE),
                                   BPPARAM = BPPARAM)
      if (progress_bar) BPPARAM$progressbar <- FALSE

      sequence.hits.rc <- bpmapply(parse_hits, sequence.hits.rc, seq.names,
                                   "-", BPPARAM = BPPARAM, SIMPLIFY = FALSE)

      sequence.hits.rc <- do.call(rbind, sequence.hits.rc)
      
      if (nrow(sequence.hits.rc) > 0) sequence.hits.rc$strand <- "-"
      match <- sequence.hits.rc$match
      match <- as.character(reverseComplement(DNAStringSet(match)))
      sequence.hits.rc$match <- match

      sequence.hits <- rbind(sequence.hits, sequence.hits.rc)

    }

    if (RNA) {
      sequence.hits$match <- gsub("T", "U", sequence.hits$match)
    }


  } else {

    if (!as.character(use.freq) %in% names(motifs@multifreq) && use.freq != 1) {
      stop("no ", use.freq, "-letter frequencies found in motif ", motifs["name"])
    }

    if (progress_bar) BPPARAM$progressbar <- TRUE
    if (progress_bar) cat("Foward strand:\n")
    results <- bpmapply(function(x, y) get_res_k(motifs, x, y, "+", use.freq,
                                                 threshold),
                        sequences, seq.names, SIMPLIFY = FALSE, BPPARAM = BPPARAM)
    if (progress_bar) BPPARAM$progressbar <- FALSE

    if (RC && motifs["alphabet"] %in% c("DNA", "RNA") &&
        sequences@elementType %in% c("DNAString", "RNAString")) {

      if (progress_bar) BPPARAM$progressbar <- TRUE
      if (progress_bar) cat("Reverse strand:\n")
      results.rc <- bpmapply(function(x, y) get_res_k(motifs, x, y, "-", use.freq,
                                                      threshold),
                             reverseComplement(sequences), seq.names,
                             SIMPLIFY = FALSE, BPPARAM = BPPARAM)
      if (progress_bar) BPPARAM$progressbar <- FALSE

      results <- c(results, results.rc)

    } else if (RC) warning("'RC = TRUE' is only supported for DNA/RNA")

    sequence.hits <- do.call(rbind, results)

  } 

  if (is.null(sequence.hits)) {
    message("no matches found using current threshold for motif ", mot.name)
    return(NULL)
  } else if (nrow(sequence.hits) == 0) {
    message("no matches found using current threshold for motif ", mot.name)
    return(NULL)
  }

  sequence.hits <- sequence.hits[order(sequence.hits$score.pct,
                                       decreasing = TRUE), ]
  rownames(sequence.hits) <- seq_len(nrow(sequence.hits))

  BPPARAM$progressbar <- pb_prev

  if (motifs["alphabet"] == "DNA") {
    sequence.hits$p.value <- vapply(sequence.hits$score,
                                    function(x) TFMsc2pv(mot.pfm, x, bkg, "PFM"),
                                    numeric(1))
    sequence.hits[, c(1:7, 10, 8:9)]
  } else if (motifs["alphabet"] == "RNA") {
    sequence.hits
  } else {
    sequence.hits[, -9]
  }


}

get_res_k <- function(motif, seq, seq.name, seqstrand, k, threshold) {

  sequence <- as.character(seq)
  sequence <- strsplit(sequence, "")[[1]]

  if (motif["pseudocount"] == 0) motif["pseudocount"] <- 0.0001
  if (length(motif["nsites"]) == 0) motif["nsites"] <- 100

  if (k == 1) {
    scores <- convert_type(motif, "PWM")["motif"]
  } else {
    scores <- apply(motif@multifreq[[as.character(k)]], 2, ppm_to_pwmC,
                    nsites = motif["nsites"],
                    pseudocount = motif["pseudocount"])
  }
  max.score <- sum(apply(scores, 2, max))
  min.score <- max.score * threshold

  max.len <- length(sequence)
  seq.mat <- matrix(ncol = max.len - k + 1, nrow = k)

  for (i in seq_len(k)) {
    to.remove <- i - 1
    if (to.remove == 0) {
      seq.mat[i, ] <- sequence[seq_len(ncol(seq.mat))]
      next
    }
    sequence.i <- sequence[-seq_len(to.remove)]
    seq.mat[i, ] <- sequence.i[seq_len(ncol(seq.mat))]
  }

  alph <- motif["alphabet"]
  if (alph == "DNA") {
    alph <- DNA_BASES
  } else if (alph == "RNA") {
    alph <- RNA_BASES
  } else if (alph == "AA") {
    alph <- AA_STANDARD
  } else if (alph == "custom") {
    alph <- unique(sequence)
  } else {
    alph <- strsplit(alph, "")[[1]]
  }
  alph.int <- as.integer(seq_len(length(alph)))

  ############
  # cpp
  seq.mat <- string_to_factor(seq.mat, alph)
  seqs <- LETTER_to_int(as.integer(seq.mat) - 1, k, alph.int)
  to_keep <- scan_seq_internal(seqs, scores, min.score)
  ############

  to_keep <- which(as.logical(to_keep))

  mot_len <- ncol(motif["motif"]) - k + 1
  name <- motif["name"]

  res <- lapply(to_keep, function(x) parse_k_res(x, sequence, seqs,
                mot_len, min.score, max.score, name, seq.name, scores,
                seqstrand, max.len, k))

  do.call(rbind, res)

}

parse_k_res <- function(to_keep, sequence, seqs, mot_len, min.score,
                        max.score, mot.name, seq.name, score.mat,
                        strand, seq.length, k) {

  if (strand == "+") {
    hit <- seqs[to_keep:(to_keep + mot_len - 1)]
    score <- score_seq(hit, score.mat)
    data.frame(motif = mot.name, sequence = seq.name, start = to_keep,
               stop = to_keep + mot_len + k - 2, score = score,
               max.score = max.score, score.pct = score / max.score * 100,
               match = paste(sequence[to_keep:(to_keep + mot_len + k - 2)],
                             collapse = ""), strand = strand)
  } else if (strand == "-") {
    hit <- seqs[to_keep:(to_keep + mot_len - 1)]
    score <- score_seq(hit, score.mat)
    start <- seq.length - to_keep + 1
    stop <- to_keep + mot_len
    stop <- seq.length - stop - k + 3
    match <- paste(sequence[to_keep:(to_keep + mot_len + k - 2)],
                   collapse = "")
    data.frame(motif = mot.name, sequence = seq.name, start = start,
               stop = stop, score = score, max.score = max.score,
               score.pct = score / max.score * 100,
               match = as.character(match), strand = strand)
  }
  
}
