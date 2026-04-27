context("scan_sequences()")

test_that("Results are accurate", {

  motif <- create_motif("AAAA", pseudocount = 1, nsites = 100)
  seq <- Biostrings::DNAStringSet("GGGAAAAGGGTTTTGGG")  # width 17
  res <- scan_sequences(motif, seq, RC = TRUE, verbose = 0)

  expect_equal(res$start[1], 4)
  expect_equal(res$motif[1], "motif")
  expect_equal(res$sequence[1], "1")
  expect_equal(res$stop[1], 7)
  expect_equal(res$match[1], "AAAA")
  expect_equal(res$strand[1], "+")
  expect_equal(res$score[1], 7.956, tolerance = 0.001)
  expect_equal(res$`max.score`[1], 7.957, tolerance = 0.001)
  expect_equal(res$`score.pct`[1], 100, tolerance = 0.1)

  expect_equal(res$score[1], res$score[2])
  expect_equal(res$`max.score`[1], res$`max.score`[2])
  expect_equal(res$strand[2], "-")
  expect_equal(res$start[2], 14)
  expect_equal(res$stop[2], 11)

  m <- create_motif(create_sequences(seqlen = 10), add.multifreq = 2, pseudocount = 1)
  s <- create_sequences()
  r <- scan_sequences(m, s, RC = TRUE, use.freq = 2, threshold = 0.8,
                      threshold.type = "logodds", verbose = 0)

  expect_true(is(r, "DataFrame"))

})

test_that("scan_sequences() streaming kernel doesn't OOM: no-hit path returns 0-row DataFrame (regression: full O(nmotifs*nseqs*seqlen) grid)", {

  set.seed(1)
  motifs <- lapply(1:30, function(i) create_motif(create_sequences(seqlen = 8, seqnum = 5)))
  seq <- create_sequences(seqlen = 5000, seqnum = 1)
  # Threshold above any achievable score => zero hits, but the old kernel
  # would still build the full 30 x 5000 score grid before filtering.
  r <- suppressWarnings(
    scan_sequences(motifs, seq, threshold = 1e9,
                   threshold.type = "logodds.abs",
                   calc.pvals = FALSE, verbose = 0)
  )
  expect_equal(nrow(r), 0)

})

test_that("scan_sequences() streaming kernel: hits pass threshold invariant (score >= thresh.score)", {

  set.seed(1)
  motifs <- lapply(1:10, function(i) create_motif(create_sequences(seqlen = 6, seqnum = 5)))
  seq <- create_sequences(seqlen = 500, seqnum = 2)
  r <- scan_sequences(motifs, seq, threshold = 0.8, threshold.type = "logodds",
                      calc.pvals = FALSE, verbose = 0)
  expect_true(nrow(r) > 0)
  expect_true(all(r$score >= r$thresh.score))

})

test_that("scan_sequences() with use.freq=3 doesn't underflow for sequences exactly as wide as the motif (regression: unsigned underflow in scan loop)", {

  # With motif width W and k=3, a sequence of length W satisfies the existing
  # width check (W >= W) but triggers the unsigned underflow:
  # W - 3 + 1 - W + 1 = -1 wraps to UINT_MAX without the guard.
  seqs <- create_sequences(seqlen = 10)
  motif <- suppressMessages(create_motif(seqs, pseudocount = 1, add.multifreq = 3))
  exact_seq <- Biostrings::DNAStringSet("ACGTACGTAC")  # exactly 10 bp == motif width
  r <- scan_sequences(motif, exact_seq, use.freq = 3, verbose = 0)
  expect_true(is(r, "DataFrame"))

})

test_that("calc.qvals.method = 'BH' matches the textbook BH formula (regression: inverted ratio + 100x bias)", {

  motif <- create_motif("ACGT", pseudocount = 1, nsites = 100)
  set.seed(1)
  seq <- create_sequences(seqlen = 200, seqnum = 5)
  r <- suppressWarnings(scan_sequences(motif, seq, threshold = 1e-3,
                      threshold.type = "pvalue", calc.pvals = TRUE,
                      calc.qvals = TRUE, calc.qvals.method = "BH",
                      verbose = 0))
  if (nrow(r) > 0) {
    mLen <- ncol(motif@motif)
    mMax <- sum(Biostrings::width(seq) - mLen + 1)
    expected <- pmin(r$pvalue * mMax / rank(r$pvalue), 1)
    expect_equal(r$qvalue, expected, tolerance = 1e-12)
  }

})
