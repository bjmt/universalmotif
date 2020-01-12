motif_clusters <- function(motifs, target.sequences, bkg.sequences,
                           min.dist.inner = 0, min.dist.outer = 20,
                           max.dist.inner = 20, max.dist.outer = 100,
                           min.size = 2, max.size = 5, allow.dup = TRUE,
                           allow.disorder = TRUE, allow.overlap = FALSE,
                           RC = FALSE, scan.threshold = 0.001,
                           scan.threshold.type = "pvalue",
                           scan.use.freq = 1, shuffle.k = 2,
                           shuffle.method = "euler",
                           motif_pvalue.k = 8, max.p = 0.001,
                           nthreads = 1, rng.seed = sample.int(1e4, 1)) {

  scan.res <- scan_sequences(motifs, target.sequences, threshold = scan.threshold,
                             threshold.type = scan.threshold.type, RC = RC,
                             use.freq = use.freq, nthreads = nthreads,
                             motif_pvalue.k = motif_pvalue.k)

  if (missing(background.sequences)) {
    bkg.sequences <- shuffle_sequences(target.sequences, k = shuffle.k,
                                       method = shuffle.method, nthreads = nthreads,
                                       rng.seed = rng.seed)
  }

  bkg.res <- scan_sequences(motifs, bkg.sequences, threshold = scan.threshold,
                            threshold.type = scan.threshold.type, RC = RC,
                            use.freq = use.freq, nthreads = nthreads,
                            motif_pvalue.k = motif_pvalue.k)

  # For every sequence, check for clusters in sliding windows max.dist.outer long

}
