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

test_that("scan_sequences() with multifreq motif on a short sequence gives a friendly error (regression: size_t underflow)", {

  m <- create_motif("ACGT", nsites = 100, pseudocount = 1)
  # Sequences the same length as the motif are used directly (no scan step)
  seqs_match <- Biostrings::DNAStringSet(rep("ACGT", 20))
  m <- add_multifreq(m, seqs_match, add.k = 2)
  short <- Biostrings::DNAStringSet("AC")
  expect_error(suppressWarnings(scan_sequences(m, short, use.freq = 2,
      verbose = 0)), regexp = "shorter than the width", fixed = TRUE)

})

test_that("scan_sequences() suggests scan_sequences2() when arguments are compatible", {

  m <- create_motif("ACGTAC")
  seqs <- create_sequences(seqnum = 3, seqlen = 100, rng.seed = 1)

  old_opt <- getOption("universalmotif.suggest.scan_sequences2")
  on.exit(options(universalmotif.suggest.scan_sequences2 = old_opt), add = TRUE)

  # Default call -- should fire.
  options(universalmotif.suggest.scan_sequences2 = TRUE)
  expect_message(
    suppressWarnings(scan_sequences(m, seqs, threshold = 1e-2,
                                    threshold.type = "pvalue", RC = TRUE,
                                    calc.pvals = FALSE, calc.qvals = FALSE,
                                    verbose = 0)),
    "scan_sequences2()", fixed = TRUE
  )

  # Each disqualifying argument should suppress the hint.
  hint_silent <- function(...) {
    msgs <- capture_messages(
      suppressWarnings(scan_sequences(m, seqs, threshold = 1e-2,
                                      threshold.type = "pvalue", RC = TRUE,
                                      calc.pvals = FALSE, calc.qvals = FALSE,
                                      verbose = 0, ...))
    )
    expect_false(any(grepl("scan_sequences2", msgs, fixed = TRUE)))
  }
  hint_silent(no.overlaps         = TRUE)
  hint_silent(respect.strand      = TRUE)
  hint_silent(motif_pvalue.method = "exhaustive")
  # allow.nonfinite is incompatible with threshold.type = "pvalue" + dynamic,
  # so verify suppression via a logodds-threshold path instead.
  msgs <- capture_messages(
    suppressWarnings(scan_sequences(m, seqs, threshold = 0.6,
                                    threshold.type = "logodds.abs", RC = TRUE,
                                    allow.nonfinite = TRUE,
                                    calc.pvals = FALSE, calc.qvals = FALSE,
                                    verbose = 0))
  )
  expect_false(any(grepl("scan_sequences2", msgs, fixed = TRUE)))

  # threshold.type other than pvalue: also suppress.
  msgs <- capture_messages(
    suppressWarnings(scan_sequences(m, seqs, threshold = 0.6,
                                    threshold.type = "logodds", RC = TRUE,
                                    calc.pvals = FALSE, calc.qvals = FALSE,
                                    verbose = 0))
  )
  expect_false(any(grepl("scan_sequences2", msgs, fixed = TRUE)))

  # The option silences the hint even when conditions are met.
  options(universalmotif.suggest.scan_sequences2 = FALSE)
  msgs <- capture_messages(
    suppressWarnings(scan_sequences(m, seqs, threshold = 1e-2,
                                    threshold.type = "pvalue", RC = TRUE,
                                    calc.pvals = FALSE, calc.qvals = FALSE,
                                    verbose = 0))
  )
  expect_false(any(grepl("scan_sequences2", msgs, fixed = TRUE)))

})

test_that("scan_sequences runs on RNA motifs and sequences", {
  m <- create_motif("ACGU", alphabet = "RNA", nsites = 10)
  seqs <- Biostrings::RNAStringSet(c("GGGACGUAAA", "UUUACGUUUU"))
  res <- scan_sequences(m, seqs, RC = FALSE, verbose = 0,
                        threshold = 0.5, threshold.type = "logodds.abs")
  expect_true(nrow(res) >= 2L)
  expect_true(all(c("start", "stop", "match") %in% colnames(res)))
})

test_that("scan_sequences runs on AA motifs and sequences", {
  m <- create_motif("YYAA", alphabet = "AA", nsites = 10)
  seqs <- Biostrings::AAStringSet(c("GGGYYAAGGG", "AAAYYAAAAA"))
  res <- scan_sequences(m, seqs, RC = FALSE, verbose = 0,
                        threshold = 0.5, threshold.type = "logodds.abs",
                        calc.pvals = FALSE)
  expect_true(nrow(res) >= 1L)
  expect_true(all(c("start", "stop", "match") %in% colnames(res)))
})

test_that("scan_sequences gives the same hits at nthreads = 1 and 2", {
  suppressMessages({
    set.seed(1)
    m <- create_motif("CACGTG", nsites = 100)
    seqs <- Biostrings::DNAStringSet(vapply(seq_len(8), function(i) {
      paste(sample(c("A","C","G","T"), 200, replace = TRUE), collapse = "")
    }, character(1)))
    a <- scan_sequences(m, seqs, nthreads = 1, verbose = 0,
                        threshold = 0.6, threshold.type = "logodds",
                        calc.pvals = FALSE)
    b <- scan_sequences(m, seqs, nthreads = 2, verbose = 0,
                        threshold = 0.6, threshold.type = "logodds",
                        calc.pvals = FALSE)
    expect_equal(as.data.frame(a), as.data.frame(b))
  })
})
