context("scan_sequences_lite()")

suppressPackageStartupMessages({
  library(Biostrings)
})

## Small reusable fixture: two DNA motifs + a handful of random sequences.
make_fixture <- function(seed = 1) {
  set.seed(seed)
  motifs <- list(
    create_motif("TATAAA", name = "M1"),
    create_motif("CACGTG", name = "M2")
  )
  seqs <- create_sequences(seqnum = 5, seqlen = 400, rng.seed = seed)
  list(motifs = motifs, seqs = seqs)
}

suppressMessages({  ## suppress "Added a pseudocount" notes from PWM coercion

test_that("scan_sequences_lite is deterministic", {
  fx <- make_fixture()
  a <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 1e-2,
                       return.granges = FALSE)
  b <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 1e-2,
                       return.granges = FALSE)
  expect_equal(a, b)
})

test_that("coordinates are always sequential", {
  fx <- make_fixture()
  hits <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 5e-2,
                          return.granges = FALSE)
  expect_true(all(hits$start <= hits$end))
  expect_true(all(hits$start >= 1L))
})

test_that("match on '-' strand reverse-complements back to the sequence", {
  fx <- make_fixture()
  hits <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 5e-2,
                          return.granges = FALSE)
  neg <- hits[hits$strand == "-", , drop = FALSE]
  skip_if(nrow(neg) == 0L, "no '-' strand hits in fixture")
  rc_match <- as.character(reverseComplement(DNAStringSet(neg$match)))
  expected <- as.character(subseq(fx$seqs[neg$sequence.i], neg$start, neg$end))
  expect_equal(rc_match, expected)
})

test_that("match on '+' strand equals the literal substring", {
  fx <- make_fixture()
  hits <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 5e-2,
                          return.granges = FALSE)
  pos <- hits[hits$strand == "+", , drop = FALSE]
  expected <- as.character(subseq(fx$seqs[pos$sequence.i], pos$start, pos$end))
  expect_equal(pos$match, expected)
})

test_that("score is in (-Inf, max.score] and matches score.pct", {
  fx <- make_fixture()
  hits <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 5e-2,
                          return.granges = FALSE)
  ## Mirror the function: normalize any motif whose PWM has -Inf entries
  ## (consensus-derived motifs), then take the sum of column maxes.
  ms <- vapply(fx$motifs, function(m) {
    pwm <- convert_type(m, "PWM")@motif
    if (any(is.infinite(pwm)))
      pwm <- convert_type(suppressMessages(normalize(m)), "PWM")@motif
    sum(apply(pwm, 2, max))
  }, numeric(1))
  expect_true(all(hits$score <= ms[hits$motif.i] + 1e-9))
  expect_true(all(hits$score.pct <= 100 + 1e-9))
  expect_equal(hits$score.pct,
               100 * hits$score / ms[hits$motif.i],
               tolerance = 1e-6)
})

test_that("pvalue is in [0, 1] for every hit", {
  fx <- make_fixture()
  hits <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 5e-2,
                          return.granges = FALSE)
  ## Underflow to 0 is legitimate at very high scores (CDF tail).
  expect_true(all(hits$pvalue >= 0))
  expect_true(all(hits$pvalue <= 1))
})

test_that("data.frame output has documented columns in order", {
  fx <- make_fixture()
  hits <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 1e-2,
                          return.granges = FALSE)
  expect_s3_class(hits, "data.frame")
  expect_identical(
    colnames(hits),
    c("motif", "motif.i", "sequence", "sequence.i",
      "start", "end", "strand",
      "score", "score.pct", "match", "pvalue")
  )
})

test_that("GRanges output carries seqlengths and mcols", {
  skip_if_not_installed("GenomicRanges")
  fx <- make_fixture()
  gr <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 1e-2,
                        return.granges = TRUE)
  expect_s4_class(gr, "GRanges")
  expect_setequal(
    colnames(GenomicRanges::mcols(gr)),
    c("motif", "motif.i", "sequence.i",
      "score", "score.pct", "match", "pvalue")
  )
  sl <- GenomeInfoDb::seqlengths(gr)
  expect_equal(unname(sl), width(fx$seqs))
})

test_that("RC = FALSE returns only '+' hits", {
  fx <- make_fixture()
  hits <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 5e-2, RC = FALSE,
                          return.granges = FALSE)
  expect_true(all(hits$strand == "+"))
})

test_that("agrees with scan_sequences() on shared hits", {
  ## scan_sequences swaps start/stop on '-' strand; harmonise before joining.
  fx <- make_fixture()
  a <- scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 1e-3,
                       return.granges = FALSE)
  b <- as.data.frame(scan_sequences(fx$motifs, fx$seqs,
                                    threshold = 1e-3,
                                    threshold.type = "pvalue",
                                    RC = TRUE, verbose = 0,
                                    calc.pvals = FALSE,
                                    calc.qvals = FALSE))
  ## Bring scan_sequences output into the same coord convention.
  neg <- b$strand == "-"
  if (any(neg)) {
    swap <- b$start[neg]
    b$start[neg] <- b$stop[neg]
    b$stop[neg]  <- swap
  }
  b$end <- b$stop
  ## Join on (motif.i, sequence.i, start, end, strand).
  ka <- paste(a$motif.i, a$sequence.i, a$start, a$end, a$strand, sep = "|")
  kb <- paste(b$motif.i, b$sequence.i, b$start, b$end, b$strand, sep = "|")
  expect_true(length(intersect(ka, kb)) > 0)
  expect_setequal(ka, kb)
})

test_that("AA motifs are rejected", {
  m <- create_motif("LLNN")
  expect_error(scan_sequences_lite(m, AAStringSet("LLNNQQAA")),
               "DNA/RNA")
})

test_that("alphabet mismatch between motif and sequence errors", {
  m <- create_motif("ACGT")
  expect_error(scan_sequences_lite(m, RNAStringSet("ACGUACGU")),
               "do not match")
})

test_that("invalid pvalue values are rejected", {
  fx <- make_fixture()
  expect_error(scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 0),  "in \\(0, 1\\)")
  expect_error(scan_sequences_lite(fx$motifs, fx$seqs, pvalue = 1),  "in \\(0, 1\\)")
  expect_error(scan_sequences_lite(fx$motifs, fx$seqs, pvalue = NA), "in \\(0, 1\\)")
  expect_error(scan_sequences_lite(fx$motifs, fx$seqs, pvalue = c(0.1, 0.2)),
               "in \\(0, 1\\)")
})

test_that("empty result returns the documented shape", {
  ## Force a no-hits situation by scanning a single sequence that does
  ## not contain either motif at any reasonable score.
  m <- list(create_motif("AAAAAA", name = "polyA"))
  seqs <- DNAStringSet(c(only = paste(rep("G", 100), collapse = "")))
  hits <- scan_sequences_lite(m, seqs, pvalue = 1e-4, return.granges = FALSE)
  expect_s3_class(hits, "data.frame")
  expect_equal(nrow(hits), 0L)
  expect_identical(
    colnames(hits),
    c("motif", "motif.i", "sequence", "sequence.i",
      "start", "end", "strand",
      "score", "score.pct", "match", "pvalue")
  )
})

test_that("empty GRanges result retains the documented metadata columns", {
  skip_if_not_installed("GenomicRanges")
  m <- list(create_motif("AAAAAA", name = "polyA"))
  seqs <- DNAStringSet(c(only = paste(rep("G", 100), collapse = "")))
  hits <- scan_sequences_lite(m, seqs, pvalue = 1e-4,
                              return.granges = TRUE)
  expect_s4_class(hits, "GRanges")
  expect_equal(length(hits), 0L)
  expect_identical(
    colnames(S4Vectors::mcols(hits)),
    c("motif", "motif.i", "sequence.i", "score", "score.pct", "match",
      "pvalue")
  )
})

})  ## end suppressMessages

test_that("scan_sequences_lite runs on RNA motifs and sequences", {
  suppressMessages({
    m <- create_motif("ACGU", alphabet = "RNA", nsites = 10)
    seqs <- Biostrings::RNAStringSet(c("GGGACGUAAA", "UUUACGUUUU"))
    res <- scan_sequences_lite(list(m), seqs, pvalue = 0.5, return.granges = FALSE)
    expect_s3_class(res, "data.frame")
    expect_true(nrow(res) >= 1L)
  })
})

test_that("scan_sequences and scan_sequences_lite both find an implanted motif (v1/v2 parity smoke)", {
  ## Different scoring stacks; this verifies both find at least one hit for an
  ## obvious implanted match. Strict numerical parity is intentionally not
  ## asserted.
  suppressMessages({
    m <- create_motif("CACGTG", name = "test", nsites = 100)
    seq <- Biostrings::DNAStringSet(paste0(paste(rep("A", 50), collapse = ""),
                                            "CACGTG",
                                            paste(rep("A", 50), collapse = "")))
    h1 <- scan_sequences(m, seq, RC = FALSE, verbose = 0,
                         threshold = 0.9, threshold.type = "logodds")
    h2 <- scan_sequences_lite(list(m), seq, pvalue = 1e-2, return.granges = FALSE)
    expect_gte(nrow(h1), 1L)
    expect_gte(nrow(h2), 1L)
    ## Both should find a hit overlapping position 51.
    expect_true(any(h1$start <= 51 & h1$stop >= 51))
    expect_true(any(h2$start <= 51 & h2$end >= 51))
  })
})

test_that("scan_sequences_lite gives the same hits at nthreads = 1 and 2", {
  suppressMessages({
    set.seed(1)
    m <- create_motif("CACGTG", nsites = 100)
    seqs <- Biostrings::DNAStringSet(vapply(seq_len(8), function(i) {
      paste(sample(c("A","C","G","T"), 200, replace = TRUE), collapse = "")
    }, character(1)))
    a <- scan_sequences_lite(list(m), seqs, pvalue = 1e-2, nthreads = 1,
                         return.granges = FALSE)
    b <- scan_sequences_lite(list(m), seqs, pvalue = 1e-2, nthreads = 2,
                         return.granges = FALSE)
    expect_equal(a, b)
  })
})

## Gapped motifs ----------------------------------------------------------

## gapped_fixture() and read_hits() are in helper-gapped.R.

## Reference values from the Rust/Python port (universalmotif-py 0.7.0,
## scan_sequences() with pvalue = 0.003), in 1-based inclusive coordinates.
port_hits_corrected <- read_hits("
u  s1 33 38 + 6.522 CACGTG       0.000244140625
u  s1 33 38 - 6.522 CACGTG       0.000244140625
g  s1 11 20 + 8.696 TGAC..GTCA   4.57763671875e-05
g  s1 11 20 - 8.696 TGAC..GTCA   4.57763671875e-05
g  s1 37 46 + 8.696 TGAC..GTCA   4.57763671875e-05
g  s1 37 46 - 8.696 TGAC..GTCA   4.57763671875e-05
g  s1 49 56 + 8.696 TGACGTCA     4.57763671875e-05
g  s1 49 56 - 8.696 TGACGTCA     4.57763671875e-05
g  s3 13 21 + 8.696 TGAC.GTCA    4.57763671875e-05
g  s3 13 21 - 8.696 TGAC.GTCA    4.57763671875e-05
g  s3 32 41 + 8.696 TGAC..GTCA   4.57763671875e-05
g  s3 32 41 - 8.696 TGAC..GTCA   4.57763671875e-05
g2 s2  9 20 + 9.783 AGG..TCA.GGA 1.52587890625e-05
g2 s2 36 45 - 9.783 AGG.TCAGGA   1.52587890625e-05
")
port_hits_uncorrected <- read_hits("
u  s1 33 38 + 6.522 CACGTG       0.000244140625
u  s1 33 38 - 6.522 CACGTG       0.000244140625
g  s1 11 20 + 8.696 TGAC..GTCA   1.52587890625e-05
g  s1 11 20 - 8.696 TGAC..GTCA   1.52587890625e-05
g  s1 37 46 + 8.696 TGAC..GTCA   1.52587890625e-05
g  s1 37 46 - 8.696 TGAC..GTCA   1.52587890625e-05
g  s1 49 56 + 8.696 TGACGTCA     1.52587890625e-05
g  s1 49 56 - 8.696 TGACGTCA     1.52587890625e-05
g  s3 13 21 + 8.696 TGAC.GTCA    1.52587890625e-05
g  s3 13 21 - 8.696 TGAC.GTCA    1.52587890625e-05
g  s3 32 41 + 8.696 TGAC..GTCA   1.52587890625e-05
g  s3 32 41 - 8.696 TGAC..GTCA   1.52587890625e-05
g2 s2  9 19 + 6.253 AGG..TCACGG  0.0013427734375
g2 s2  9 20 + 9.783 AGG..TCA.GGA 3.814697265625e-06
g2 s2 10 20 + 6.253 GGT.TCA.GGA  0.0013427734375
g2 s2 28 39 - 6.253 AGG..TAA.CGA 0.0013427734375
g2 s2 35 45 - 6.253 AGG.TCA.GAT  0.0013427734375
g2 s2 36 45 - 9.783 AGG.TCAGGA   3.814697265625e-06
g2 s2 36 46 - 6.253 CAG..TCAGGA  0.0013427734375
")

expect_hits <- function(res, ref) {
  res <- as.data.frame(res)[, colnames(ref)]
  rownames(res) <- rownames(ref) <- NULL
  cols <- c("motif", "sequence", "start", "end", "strand", "match")
  expect_identical(res[, cols], ref[, cols])
  expect_equal(res$score, ref$score, tolerance = 1e-9)
  expect_equal(res$pvalue, ref$pvalue, tolerance = 1e-9)
}

scan_gapped_fixture <- function(...) {
  f <- gapped_fixture()
  suppressMessages(scan_sequences_lite(f$motifs, f$seqs, return.granges = FALSE,
                                       ...))
}

test_that("gapped scans agree with the port's reference values", {
  expect_hits(scan_gapped_fixture(pvalue = 0.003), port_hits_corrected)
  expect_hits(scan_gapped_fixture(pvalue = 0.003, correct.gaps = FALSE),
              port_hits_uncorrected)
  ## The forward strand alone.
  expect_hits(scan_gapped_fixture(pvalue = 0.003, correct.gaps = FALSE,
                                  RC = FALSE),
              port_hits_uncorrected[port_hits_uncorrected$strand == "+", ])
})

test_that("correct.gaps scales gapped P-values and cutoffs by the number of combinations", {
  combos <- c(u = 1, g = 3, g2 = 4)
  raw <- scan_gapped_fixture(pvalue = 0.01, correct.gaps = FALSE)
  cor <- scan_gapped_fixture(pvalue = 0.01)
  expect_lt(nrow(cor), nrow(raw))
  key <- function(x) paste(x$motif, x$sequence, x$start, x$end, x$strand,
                           x$match)
  m <- match(key(cor), key(raw))
  expect_false(anyNA(m))
  expect_equal(cor$score, raw$score[m])
  expect_equal(cor$pvalue, unname(pmin(raw$pvalue[m] * combos[cor$motif], 1)))
  keep <- raw$pvalue * combos[raw$motif] <= 0.01
  expect_setequal(key(cor), key(raw)[keep])
})

test_that("use.gaps = FALSE scans gapped motifs without their gaps", {
  res <- scan_gapped_fixture(pvalue = 0.003, use.gaps = FALSE)
  expect_hits(res, read_hits("
u  s1 33 38 + 6.522 CACGTG    0.000244140625
u  s1 33 38 - 6.522 CACGTG    0.000244140625
g  s1 49 56 + 8.696 TGACGTCA  1.52587890625e-05
g  s1 49 56 - 8.696 TGACGTCA  1.52587890625e-05
g2 s2 36 44 - 6.253 GGTTCAGGA 0.0013427734375
"))
  f <- gapped_fixture()
  ungapped <- lapply(f$motifs, function(m) {
    m@gapinfo@isgapped <- FALSE
    m
  })
  ref <- suppressMessages(scan_sequences_lite(ungapped, f$seqs, pvalue = 0.003,
                                              return.granges = FALSE))
  expect_identical(res, ref)
})

test_that("ungapped motifs scan the same next to gapped motifs", {
  f <- gapped_fixture()
  run <- function(motifs, ...) suppressMessages(scan_sequences_lite(motifs,
    f$seqs, pvalue = 0.01, return.granges = FALSE, ...))
  alone <- run(f$motifs[1])
  for (correct in c(TRUE, FALSE)) {
    mixed <- run(f$motifs, correct.gaps = correct)
    mixed <- mixed[mixed$motif == "u", ]
    rownames(mixed) <- NULL
    expect_identical(mixed, alone)
  }
  ## Gapped hits are the same whatever comes before them.
  both <- run(f$motifs)
  both <- both[both$motif != "u", ]
  both$motif.i <- both$motif.i - 1L
  rownames(both) <- NULL
  expect_identical(both, run(f$motifs[2:3]))
})

test_that("gapped hits at the same start are listed by strand and then by combination", {
  f <- gapped_fixture()
  res <- suppressMessages(scan_sequences_lite(f$motifs[2], f$seqs,
                                              pvalue = 0.01, correct.gaps = FALSE,
                                              return.granges = FALSE))
  at <- res[res$sequence == "s1" & res$start == 32, ]
  expect_identical(at$strand, c("+", "+", "-", "-"))
  expect_identical(at$end, c(39L, 41L, 39L, 41L))
  expect_identical(at$match,
                   c("TCACGTGA", "TCAC..GACA", "TCACGTGA", "TGTC..GTGA"))
})

test_that("gapped GRanges results span the gaps", {
  f <- gapped_fixture()
  gr <- suppressMessages(scan_sequences_lite(f$motifs, f$seqs, pvalue = 0.003,
                                             return.granges = TRUE))
  df <- scan_gapped_fixture(pvalue = 0.003)
  expect_s4_class(gr, "GRanges")
  expect_identical(BiocGenerics::end(gr), df$end)
  expect_identical(S4Vectors::mcols(gr)$match, df$match)
})
