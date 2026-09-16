context("inclusive score and P-value thresholds")

threshold_fixture <- function(width = 4L) {
  pwm <- matrix(rep(c(1, 0, -1, -1), width), nrow = 4L,
                dimnames = list(c("A", "C", "G", "T"), NULL))
  create_motif(pwm, type = "PWM", alphabet = "DNA", pseudocount = 0)
}

test_that("dynamic thresholds select exactly the eligible integer-score support", {
  m <- threshold_fixture()
  words <- as.matrix(expand.grid(rep(list(seq_len(4L)), 4L)))
  scores <- rowSums(matrix(c(1, 0, -1, -1)[words], ncol = 4L))
  support <- sort(unique(scores))
  tails <- vapply(support, function(s) mean(scores >= s), numeric(1))
  cutoffs <- unique(c(1e-4, 1/256, .01, 5/256, .1, tails, 1))
  for (threads in c(1L, 2L)) {
    thresholds <- motif_pvalue(m, pvalue = cutoffs, nthreads = threads)
    for (i in seq_along(cutoffs)) {
      expect_identical(support[support >= thresholds[i]], support[tails <= cutoffs[i]])
      if (is.finite(thresholds[i]))
        expect_lte(motif_pvalue(m, score = thresholds[i]), cutoffs[i])
      else expect_identical(thresholds[i], Inf)
    }
  }
})

test_that("dynamic score queries respect inclusive decimal boundaries", {
  m <- threshold_fixture(2L)
  queries <- c(-1.001, -1, -.9999, 0, .0001, .9999, 1, 1.001)
  words <- as.matrix(expand.grid(rep(list(seq_len(4L)), 2L)))
  scores <- rowSums(matrix(c(1, 0, -1, -1)[words], ncol = 2L))
  expected <- vapply(queries, function(q) mean(scores >= q), numeric(1))
  expect_equal(motif_pvalue(m, score = queries), expected, tolerance = 1e-14)
})

test_that("both scanners enforce P-value cutoffs including impossible and equal cutoffs", {
  m <- threshold_fixture()
  seqs <- Biostrings::DNAStringSet(c("AAAA", "AAAC", "TTTT", "GTTT"))
  scan <- function(engine, cutoff, threads, motifs = list(m)) {
    if (engine == "lite")
      scan_sequences_lite(motifs, seqs, pvalue = cutoff, RC = TRUE,
        return.granges = FALSE, nthreads = threads)
    else as.data.frame(suppressMessages(scan_sequences(motifs, seqs,
      threshold = cutoff, RC = TRUE, calc.pvals = TRUE, calc.qvals = FALSE,
      return.granges = FALSE, nthreads = threads, verbose = 0)))
  }
  for (engine in c("standard", "lite")) for (threads in c(1L, 2L)) {
    hits <- scan(engine, .01, threads)
    expect_equal(nrow(hits), 2L)
    expect_true(all(hits$match == "AAAA"))
    expect_true(all(hits$pvalue <= .01))
    expect_setequal(hits$strand, c("+", "-"))
    expect_equal(nrow(scan(engine, 5/256, threads)), 4L)
    expect_equal(nrow(scan(engine, 1/256, threads)), 2L)
    empty <- scan(engine, 1e-4, threads)
    expect_equal(nrow(empty), 0L)
    expect_true(all(c("score", "pvalue", "motif.i", "sequence.i") %in% names(empty)))
    # Skipping an impossible motif must preserve indices of the others.
    mixed <- scan(engine, .001, threads,
      motifs = list(m, create_motif("AAAA", bkg = c(.1, .4, .4, .1),
                                    pseudocount = 1, nsites = 100)))
    expect_true(nrow(mixed) > 0L)
    expect_true(all(mixed$motif.i == 2L))
    expect_true(all(mixed$pvalue <= .001))
  }
})

test_that("absolute score thresholds never admit a lower score", {
  m <- threshold_fixture(2L)
  seqs <- Biostrings::DNAStringSet(c("AA", "AC", "CC", "CG", "GG"))
  for (cutoff in c(-1.001, -1, -.9999, 0, .0001, .9999, 1, 1.001)) {
    hits <- suppressMessages(scan_sequences(m, seqs, threshold = cutoff,
      threshold.type = "logodds.abs", RC = FALSE, calc.pvals = FALSE, verbose = 0))
    expected <- c(2, 1, 0, -1, -2)
    expect_equal(sort(hits$score), sort(expected[expected >= cutoff]))
  }
})

test_that("integer-grid round trips preserve equality on both sides of zero", {
  # These integer scores exercise decimal-to-binary representation boundaries.
  pwm <- matrix(c(.501, 0, -.501, -1, .5, 0, -.5, -1), 4L,
                dimnames = list(c("A", "C", "G", "T"), NULL))
  m <- create_motif(pwm, type = "PWM", alphabet = "DNA", pseudocount = 0)
  quantised <- as.vector(outer(trunc(pwm[, 1] * 1000),
                               trunc(pwm[, 2] * 1000), "+")) / 1000
  for (score in quantised) {
    expected <- mean(quantised >= score)
    expect_equal(motif_pvalue(m, score = score), expected, tolerance = 1e-14)
  }
  hits <- suppressMessages(scan_sequences(m, Biostrings::DNAStringSet("AA"),
    threshold = 1.001, threshold.type = "logodds.abs", RC = FALSE,
    calc.pvals = FALSE, verbose = 0))
  expect_equal(hits$score, 1.001)
})

test_that("empty GRanges retain metadata with an impossible P-value cutoff", {
  skip_if_not_installed("GenomicRanges")
  m <- threshold_fixture()
  seqs <- Biostrings::DNAStringSet(c(a = "AAAA", b = "TTTT"))
  a <- suppressMessages(scan_sequences(m, seqs, threshold = 1e-4,
    return.granges = TRUE, calc.pvals = TRUE, verbose = 0))
  b <- scan_sequences_lite(m, seqs, pvalue = 1e-4, return.granges = TRUE)
  for (result in list(a, b)) {
    expect_length(result, 0L)
    expect_true(all(c("score", "pvalue") %in% colnames(S4Vectors::mcols(result))))
    expect_equal(unname(GenomeInfoDb::seqlengths(result)), c(4L, 4L))
  }
})

test_that("an explicit changed background still rebases PWM scores", {
  m <- threshold_fixture(2L)
  bg <- c(A = .1, C = .4, G = .4, T = .1)
  pwm <- log2(matrix(rep(c(.5, .25, .125, .125), 2L), nrow = 4L) / bg)
  words <- as.matrix(expand.grid(rep(list(seq_len(4L)), 2L)))
  scores <- trunc(pwm[cbind(words[, 1], 1L)] * 1000) +
    trunc(pwm[cbind(words[, 2], 2L)] * 1000)
  weights <- bg[words[, 1]] * bg[words[, 2]]
  queries <- c(-2.5, 0, 3)
  expected <- vapply(queries, function(q) sum(weights[scores >= ceiling(q * 1000)]),
                     numeric(1))
  expect_equal(motif_pvalue(m, score = queries, bkg.probs = bg), expected,
               tolerance = 1e-12)
})
