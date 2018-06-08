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

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, "PWM", BPPARAM = BPPARAM)

  if (is.list(motifs)) {
    results <- lapply(motifs, scan_sequences, threshold = threshold,
                      BPPARAM = BPPARAM)
    results <- do.call(rbind, results)
    return(results)
  }

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))

  motif <- motifs

  mot.name <- motif["name"]
  mot.mat <- motif["motif"]
  max.score <- sum(apply(mot.mat, 2, max))
  min.score <- max.score * threshold
  mot.len <- ncol(mot.mat)

  mot.mat.rc <- NULL
  if (RC) mot.mat.rc <- motif_rc(motif)["motif"]

  get_score <- function(x) {
    vapply(1:mot.len, function(y) mot.mat[rownames(mot.mat) == x[y], y],
           numeric(1))
  }

  get_scoreRC <- function(x) {
    vapply(1:mot.len, function(y) mot.mat.rc[rownames(mot.mat.rc) == x[y], y],
           numeric(1))
  }

  sequences <- as.matrix(sequences)

  sequences.split <- lapply(1:nrow(sequences), function(x) sequences[x, ])

  parse_sequences <- function(x) {
    toadd <- mot.len - 1
    out <- matrix(ncol = mot.len, nrow = length(x) - toadd)
    for (i in 1:(length(x) - toadd)) {
      out[i, ] <- x[i:(i + toadd)]
    }
    out
  }

  ## + strand

  sequences.split <- bplapply(sequences.split, parse_sequences,
                              BPPARAM = BPPARAM)

  if (RC) sequences.split.rc <- sequences.split

  sequence.scores <- bplapply(sequences.split,
                              function(x) t(apply(x, 1, get_score)),
                              BPPARAM = BPPARAM)

  sequence.scores <- bplapply(sequence.scores, rowSums, BPPARAM = BPPARAM)

  sequence.matches <- bplapply(sequence.scores, function(x) x >= min.score,
                               BPPARAM = BPPARAM)

  for (i in seq_along(sequences.split)) {
    sequences.split[[i]] <- sequences.split[[i]][sequence.matches[[i]], ]
  }

  for (i in seq_along(sequence.scores)) {
    sequence.scores[[i]] <- sequence.scores[[i]][sequence.matches[[i]]]
  }

  sequence.locations <- bplapply(sequence.matches, which,
                                 BPPARAM = BPPARAM)

  sequences.split <- bplapply(sequences.split, 
                              function(x) apply(x, 1, paste, collapse = ""),
                              BPPARAM = BPPARAM)

  collate_results <- function(n, x, y, z) {
    if (length(y) == 0) {
      data.frame(motif = NULL, sequence = NULL, start = NULL, stop = NULL,
                 max.score = NULL, score.pct = NULL, match = NULL)
    } else {
      data.frame(motif = rep(mot.name, length(y)), sequence = n, start = x,
                 stop = x + (mot.len - 1),
                 score = y, max.score = max.score,
                 score.pct = (y / max.score) * 100, match = z)  
    }
  }

  results.table <- mapply(collate_results,
                          seq.names, sequence.locations,
                          sequence.scores, sequences.split,
                          SIMPLIFY = FALSE)

  results.table <- do.call(rbind, results.table)
  if (nrow(results.table) > 0) results.table$strand <- "+"

  ## - strand

  if (RC) {

    sequence.scores.rc <- bplapply(sequences.split.rc,
                                   function(x) t(apply(x, 1, get_scoreRC)),
                                   BPPARAM = BPPARAM)
    sequence.scores.rc <- bplapply(sequence.scores.rc, rowSums, BPPARAM = BPPARAM)

    sequence.matches.rc <- bplapply(sequence.scores.rc, function(x) x >= min.score,
                                    BPPARAM = BPPARAM)

    for (i in seq_along(sequences.split.rc)) {
      sequences.split.rc[[i]] <- sequences.split.rc[[i]][sequence.matches.rc[[i]], ]
    }

    for (i in seq_along(sequence.scores.rc)) {
      sequence.scores.rc[[i]] <- sequence.scores.rc[[i]][sequence.matches.rc[[i]]]
    }

    sequence.locations.rc <- bplapply(sequence.matches.rc, which,
                                      BPPARAM = BPPARAM)

    sequences.split.rc <- bplapply(sequences.split.rc, 
                                   function(x) apply(x, 1, paste, collapse = ""),
                                   BPPARAM = BPPARAM)

    results.table.rc <- mapply(collate_results,
                               seq.names, sequence.locations.rc,
                               sequence.scores.rc, sequences.split.rc,
                               SIMPLIFY = FALSE)

    results.table.rc <- do.call(rbind, results.table.rc)
    if (nrow(results.table.rc) > 0) results.table.rc$strand <- "-"
  
  }

  if (RC) results.table <- rbind(results.table, results.table.rc)

  if (nrow(results.table) == 0) {
    message("no matches found using current threshold")
    return(NULL)
  }

  rownames(results.table) <- seq_len(nrow(results.table))

  results.table

}
