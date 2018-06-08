#' Find motif binding sites in a set of sequences.
#'
#' @param motifs List of motifs or a single motif.
#' @param sequences XStringSet object. List of sequences to scan
#' @param threshold Numeric. Logodds threshold.
#' @param RC Logical. Check reverse strand.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Site search results as a data.frame object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.6,
                           RC = TRUE, BPPARAM = bpparam()) {

  if (is.list(motifs)) {
    results <- lapply(motifs, scan_sequences, threshold = threshold,
                      sequences = sequences, RC = RC, BPPARAM = BPPARAM)
    results <- do.call(rbind, results)
    return(results)
  }

  motif <- motifs
  motif <- convert_type(motif,"PWM")
  motif.rc <- motif_rc(motif)
  motif <- motif["motif"]
  motif.rc <- motif.rc["motif"]

  min.score <- threshold * 100
  min.score <- paste0(min.score, "%")

  mot.name <- motifs["name"]
  mot.mat <- convert_type(motifs, "PWM")["motif"]
  mot.len <- ncol(mot.mat)
  max.score <- sum(apply(mot.mat, 2, max))

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))

  parse_hits <- function(x, y) {
    if (length(x) == 0) {
      data.frame(motif = NULL, sequence = NULL, start = NULL, stop = NULL,
                 max.score = NULL, score.pct = NULL, match = NULL,
                 strand = NULL)
    } else {
      motif <- rep(mot.name, length(x))
      sequence = rep(y, length(x))
      start <- x@ranges@start
      stop <- start + (mot.len - 1)
      score <- x@elementMetadata$score 
      max.score <- rep(max.score, length(x))
      score.pct <- (score / max.score) * 100
      match <- as.character(x)
      data.frame(motif = motif, sequence = sequence, start = start, stop = stop,
                 score = score, max.score = max.score, score.pct = score.pct,
                 match = match, strand = rep(NA, length(x)))
    }
  }

  ## + strand

  sequence.hits <- bplapply(seq_len(length(sequences)),
                            function(x) matchPWM(motif, sequences[[x]],
                                                 min.score = min.score,
                                                 with.score = TRUE),
                            BPPARAM = BPPARAM)

  sequence.hits <- bpmapply(parse_hits, sequence.hits, seq.names,
                            BPPARAM = BPPARAM, SIMPLIFY = FALSE)

  sequence.hits <- do.call(rbind, sequence.hits)
  if (nrow(sequence.hits) > 0) sequence.hits$strand <- "+"

  ## - strand

  if (RC) {
    sequence.hits.rc <- bplapply(seq_len(length(sequences)),
                                 function(x) matchPWM(motif.rc, sequences[[x]],
                                                      min.score = min.score,
                                                      with.score = TRUE),
                                 BPPARAM = BPPARAM)

    sequence.hits.rc <- bpmapply(parse_hits, sequence.hits.rc, seq.names,
                                 BPPARAM = BPPARAM, SIMPLIFY = FALSE)

    sequence.hits.rc <- do.call(rbind, sequence.hits.rc)
    
    if (nrow(sequence.hits.rc) > 0) sequence.hits.rc$strand <- "-"

    sequence.hits <- rbind(sequence.hits, sequence.hits.rc)
  }

  if (nrow(sequence.hits) == 0) {
    message("no matches found using current threshold for motif ", mot.name)
    return(NULL)
  }

  rownames(sequence.hits) <- seq_len(nrow(sequence.hits))

  sequence.hits

}
