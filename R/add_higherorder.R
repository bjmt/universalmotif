#' Add multi-nucleotide information to a motif.
#'
#' @param motif universalmotif object.
#' @param sequences DNAStringSet.
#' @param add.k Numeric.
#' @param threshold Numeric.
#' @param RC Logical.
#' @param motifs.perseq Numeric.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
add_higherorder <- function(motif, sequences, add.k = 2:3,
                            threshold = 0.01, RC = FALSE,
                            motifs.perseq = 1, BPPARAM = bpparam()) {

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

  for (i in seq_along(add.k)) {
    motif@multifreq[[i]] <- add_multi(motif["bkg"], sequences, add.k[i])
  }
  names(motif@multifreq) <- add.k

  motif

}
