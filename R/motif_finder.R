motif_finder <- function(sequences, bkg.sequences = NULL, nmotifs = 5,
  max.p = 1e-6, min.nsites = as.integer(length(sequences) * 0.2),
  starting.sizes = c(8, 10, 12), RC = TRUE, extend.motifs = TRUE,
  trim.motifs = TRUE, min.edge.ic = 0.5,
  pseudocount = 1, nthreads = 1) {

  # add option for min motif ambiguity based on average per position IC

  if (!seqtype(sequences) %in% c("DNA", "RNA"))
    stop("Only DNA/RNA alphabets are currently supported.", call. = FALSE)

  alph <- switch(seqtype(sequences),
    "DNA" = DNA_BASES,
    "RNA" = RNA_BASES
  )
  alph <- collapse_cpp(alph)

  if (RC && !seqtype(sequences) %in% c("DNA", "RNA"))
    stop("`RC=TRUE` can only be used with DNA/RNA alphabets", call. = FALSE)

  message("Calculating sequence [", paste0(starting.sizes, collapse = ", "),
    "]-mer content...")

  k1 <- get_bkg(sequences, k = c(1, starting.sizes), nthreads = nthreads, RC = RC)
  seqsk <- k1[nchar(k1$klet) > 1, ]
  k1 <- k1[nchar(k1$klet) == 1, ]
  k1 <- structure(k1$probability, names = k1$klet)

  countsums <- tapply(seqsk$count, nchar(seqsk$klet), sum)
  countsums <- structure(as.vector(countsums), names = names(countsums))

  message("Calculating background [", paste0(starting.sizes, collapse = ", "),
    "]-mer content...")

  if (is.null(bkg.sequences)) {
    seqsk$nullProb <- calc_seq_probs_cpp(seqsk$klet, k1, alph, nthreads)
    seqsk$nullCount <- seqsk$nullProb * countsums[as.character(nchar(seqsk$klet))]
  } else {
    bkgCounts <- get_bkg(bkg.sequences, k = starting.sizes, nthreads = nthreads, RC = RC)
    seqsk$nullCount <- bkgCounts$count
    seqsk$nullProb <- bkgCounts$probability
  }

  message("Looking for over-represented [", paste0(starting.sizes, collapse = ", "),
    "]-mers...")

  seqsk$Ratio <- (seqsk$count + pseudocount) / (seqsk$nullCount + pseudocount)
  # seqsk$Zscore <- scale(seqsk$Ratio)
  # seqsk$log10GaussP <- as.vector(pnorm(seqsk$Zscore, lower.tail = FALSE,
  #     log.p = TRUE)) / log(10)

  seqsk$log10BinomP <- as.vector(pbinom(seqsk$count + pseudocount - 1,
      countsums[as.character(nchar(seqsk$klet))] + pseudocount,
      seqsk$nullProb + pseudocount / countsums[as.character(nchar(seqsk$klet))],
      lower.tail = FALSE, log.p = TRUE)) / log(10)

  if (RC) {
    # get rid of duplicate rows
    rows_rc <- as.character(reverseComplement(DNAStringSet(seqsk$klet)))
    rows_uq <- pmin(seqsk$klet, rows_rc)
    seqsk <- seqsk[!duplicated(rows_uq), ]
  }

  # seqsk <- seqsk[seqsk$log10GaussP < log10(max.p), ]
  seqsk <- seqsk[seqsk$log10BinomP < log10(max.p), ]
  # seqsk <- seqsk[seqsk$Ratio >= 2, ]  # delete this?
  seqsk <- seqsk[seqsk$count >= min.nsites, ]

  if (!nrow(seqsk)) {
    message("No over-represented [", paste0(starting.sizes, collapse = ", "),
      "]-mers found.")
    return(invisible())
  }

  message("Found ", nrow(seqsk), " over-represented [",
    paste0(starting.sizes, collapse = ", "), "]-mers.")

  seqsk <- seqsk[order(seqsk$log10BinomP), ]

  # next: merge overlapping klets?

  # Next step is to find motifs within each k-mer size.

  print(seqsk)
  suppressMessages(scan_sequences(
    lapply(seqsk$klet, function(x) create_motif(x, name = x)),
    sequences, threshold = 1, RC = RC, nthreads = nthreads, warn.NA = FALSE,
    threshold.type = "logodds", return.granges = TRUE)[, "motif"])

  # Maybe: starting from the top of the list, try merging with all subsequent
  # k-mers and testing whether that increases the enrichment p-value. Once
  # finished, remove those k-mers from the list and repeat until desired # of
  # motifs is achieved.

  # Maybe go through this whole process one size k-mer at a time.

  # Give more weight to top motif when merging?

  # motif_pvalue(create_motif("AAAAAA"), , 10^-6) --> no mismatches allowed

  # What to do about k-mer hits which are just 1 bp shifts? These artificially
  # inflate the counts. Maybe after creating enriched k-mer table, go through
  # and clean up inflated counts. If bkg != NULL then clean up bkg as well.
  # (Keep in mind strand?)

}
