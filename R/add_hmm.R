#' Add first and/or second order HMM data to a motif.
#'
#' @param motif universalmotif object.
#' @param sequences DNAStringSet.
#' @param add.first Logical.
#' @param add.second Logical.
#' @param threshold Numeric.
#' @param RC Logical.
#' @param motifs.perseq Numeric.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
add_hmm <- function(motif, sequences, add.first = TRUE, add.second = TRUE,
                    threshold = 0.01, RC = FALSE, motifs.perseq = 1,
                    BPPARAM = bpparam()) {

  motif <- convert_motifs(motif, BPPARAM = BPPARAM)
  
  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))

  seq.res <- scan_sequences(motif, sequences, threshold = threshold, RC = RC,
                            BPPARAM = BPPARAM)
  seqs.out <- vector("list", length(sequences))

  for (i in seq_len(length(sequences))) {
    seq.out <- seq.res[seq.res$sequence == seq.names[i], ]
    seq.out <- seq.out[order(seq.out$score, decreasing = TRUE), ]
    seq.out <- seq.out[seq_len(motifs.perseq), ]
    seqs.out[[i]] <- seq.out
  }

  seqs.out <- do.call(cbind, seqs.out)
  seqs.out <- seqs.out$match
  seqs.out <- levels(seqs.out)
  sequences <- DNAStringSet(seqs.out)

  if (add.first) {
    motif@hmmfirst <- .my_create_first(motif["bkg"], sequences)
  }
  if (add.second) {
    motif@hmmsecond <- .my_create_second(motif["bkg"], sequences)
  }

  motif

}
