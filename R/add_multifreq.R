#' Add multi-nucleotide information to a motif.
#'
#' @param motif universalmotif object.
#' @param sequences DNAStringSet.
#' @param add.k Numeric.
#' @param threshold Numeric.
#' @param threshold.type Character.
#' @param RC Logical.
#' @param motifs.perseq Numeric.
#' @param BPPARAM See \code{\link[BiocParallel]{SerialParam}}.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
add_multifreq <- function(motif, sequences, add.k = 2:3, RC = FALSE,
                          threshold = 0.01, threshold.type = "logodds",
                          motifs.perseq = 1, BPPARAM = SerialParam()) {

  motif <- convert_motifs(motif, BPPARAM = BPPARAM)
  
  if (all(ncol(motif["motif"])) != unique(width(sequences))) {

    seq.names <- names(sequences)
    if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
    seq.res <- scan_sequences(motif, sequences, threshold = threshold, RC = RC,
                              threshold.type = threshold.type, BPPARAM = BPPARAM)

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

  }

  if (sequences@elementType == "DNAString") {
    sequences <- DNAStringSet(seqs.out)
    for (i in seq_along(add.k)) {
      motif@multifreq[[i]] <- add_multi(motif["bkg"], sequences, add.k[i])
    }
  } else {
    alph <- rownames(motif["motif"])
    for (i in seq_along(add.k)) {
      motif@multifreq[[i]] <- add_multi_ANY(sequences, add.k[i], alph)
    }
  }

  names(motif@multifreq) <- add.k

  motif

}
