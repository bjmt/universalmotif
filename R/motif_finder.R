motif_finder <- function(sequences, bkg.sequences = NULL, nmotifs = 3,
  min.size = 5, max.size = 12, max.p = 1e-3, min.edge.ic = 0.5,
  min.nsites = as.integer(width(sequences) * 0.2), seed.k = 6:8, RC = TRUE,
  pseudocount = 1, nthreads = 1) {

  alph <- switch(seqtype(sequences),
    "DNA" = DNA_BASES,
    "RNA" = RNA_BASES,
    "AA" = AA_STANDARD2,
    stop("incomplete")
  )
  alph <- collapse_cpp(alph)

  k1 <- get_bkg(sequences, k = c(1, seed.k), nthreads = nthreads)
  seqsk <- k1[nchar(k1$klet) > 1, ]
  k1 <- k1[nchar(k1$klet) == 1, ]
  k1 <- structure(k1$probability, names = k1$klet)

  countsums <- tapply(seqsk$count, nchar(seqsk$klet), sum)
  countsums <- structure(as.vector(countsums), names = names(countsums))

  if (is.null(bkg.sequences)) {
    seqsk$nullProb <- calc_seq_probs_cpp(seqsk$klet, k1, alph, nthreads)
    seqsk$nullCount <- seqsk$nullProb * countsums[as.character(nchar(seqsk$klet))]
  } else {
    bkgCounts <- get_bkg(bkg.sequences, k = seed.k, nthreads = nthreads)
    seqsk$nullCount <- bkgCounts$count
    seqsk$nullProb <- bkgCounts$probability
  }

  seqsk$Ratio <- (seqsk$count + pseudocount) / (seqsk$nullCount + pseudocount)
  seqsk$Zscore <- scale(seqsk$Ratio)
  # seqsk$Zscore <- scale(log2(seqsk$Ratio))
  seqsk$log10GaussP <- as.vector(pnorm(seqsk$Zscore, lower.tail = FALSE,
      log.p = TRUE)) / log(10)

  seqsk$log10BinomP <- as.vector(pbinom(seqsk$count + pseudocount - 1,
      countsums + pseudocount, seqsk$nullProb + pseudocount / countsums,
      lower.tail = FALSE, log.p = TRUE)) / log(10)

  seqsk <- seqsk[seqsk$log10GaussP < log10(max.p), ]

  # Binom vs Gauss: binom favours cases with high counts in target seqs, even
  # if that results in a lower ratio. Gauss ranks purely based on the ratio
  # of counts between target and bkg.

  seqsk[order(seqsk$log10BinomP), ]

  # next: merge overlapping klets?
  # nsites cutoff: do this after clustering?

}
