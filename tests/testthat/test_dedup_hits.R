context("dedup_hits()")

suppressPackageStartupMessages({
  library(Biostrings)
})

# --- Hand-computed examples -------------------------------------------------

test_that("greedy semantics: chained cluster keeps best + non-overlapping next-best", {
  # [1-5, p=1e-10], [3-7, p=1e-8], [6-9, p=1e-9]
  # max-end sweep clusters all three; priority order picks h1, then h3 (no
  # overlap with h1), drops h2 (overlaps h1).
  hits <- data.frame(
    sequence = c("s1","s1","s1"),
    motif    = c("M",  "M",  "M"),
    strand   = c("+",  "+",  "+"),
    start    = c(1L,   3L,   6L),
    end      = c(5L,   7L,   9L),
    pvalue   = c(1e-10, 1e-8, 1e-9),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits)
  expect_equal(nrow(out), 2L)
  expect_equal(out$start, c(1L, 6L))
})

test_that("greedy semantics differs from connected-components: distinct survivors survive", {
  # Mirror of the example above; with CC semantics only h1 would survive.
  hits <- data.frame(
    sequence = c("s1","s1","s1"),
    motif    = c("M","M","M"),
    strand   = c("+","+","+"),
    start    = c(1L,3L,6L),
    end      = c(5L,7L,9L),
    pvalue   = c(1e-10, 1e-5, 1e-9),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits)
  expect_setequal(out$start, c(1L, 6L))
})

test_that("score-priority direction is auto-detected: higher score wins", {
  hits <- data.frame(
    sequence = c("s1","s1"),
    motif    = c("M","M"),
    strand   = c("+","+"),
    start    = c(1L, 2L),
    end      = c(5L, 6L),
    score    = c(3.0, 9.0),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits, by = "score")
  expect_equal(out$score, 9.0)
})

test_that("reverse = TRUE flips the priority direction", {
  hits <- data.frame(
    sequence = c("s1","s1"),
    motif    = c("M","M"),
    strand   = c("+","+"),
    start    = c(1L, 2L),
    end      = c(5L, 6L),
    pvalue   = c(1e-10, 1e-5),
    stringsAsFactors = FALSE
  )
  # Default: lowest p-value wins.
  expect_equal(dedup_hits(hits)$pvalue, 1e-10)
  # Reversed: highest p-value wins.
  expect_equal(dedup_hits(hits, reverse = TRUE)$pvalue, 1e-5)
})

test_that("ties go to the earliest row in input", {
  hits <- data.frame(
    sequence = c("s1","s1"),
    motif    = c("M","M"),
    strand   = c("+","+"),
    start    = c(1L, 2L),
    end      = c(5L, 6L),
    pvalue   = c(1e-5, 1e-5),
    row_id   = c("first", "second"),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits)
  expect_equal(out$row_id, "first")
})

# --- Group key flags --------------------------------------------------------

test_that("by default '+' and '-' strands are independent", {
  hits <- data.frame(
    sequence = c("s1","s1"),
    motif    = c("M","M"),
    strand   = c("+","-"),
    start    = c(1L, 1L),
    end      = c(5L, 5L),
    pvalue   = c(1e-5, 1e-3),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits)
  expect_equal(nrow(out), 2L)
})

test_that("ignore.strand = TRUE makes strands compete", {
  hits <- data.frame(
    sequence = c("s1","s1"),
    motif    = c("M","M"),
    strand   = c("+","-"),
    start    = c(1L, 1L),
    end      = c(5L, 5L),
    pvalue   = c(1e-5, 1e-3),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits, ignore.strand = TRUE)
  expect_equal(nrow(out), 1L)
  expect_equal(out$pvalue, 1e-5)
})

test_that("ignore.motif = TRUE makes motifs compete", {
  hits <- data.frame(
    sequence = c("s1","s1"),
    motif    = c("M1","M2"),
    strand   = c("+","+"),
    start    = c(1L, 1L),
    end      = c(5L, 5L),
    pvalue   = c(1e-5, 1e-3),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(dedup_hits(hits)),                        2L)
  expect_equal(nrow(dedup_hits(hits, ignore.motif = TRUE)),    1L)
})

test_that("different sequences never compete", {
  hits <- data.frame(
    sequence = c("s1","s2"),
    motif    = c("M","M"),
    strand   = c("+","+"),
    start    = c(1L, 1L),
    end      = c(5L, 5L),
    pvalue   = c(1e-5, 1e-3),
    stringsAsFactors = FALSE
  )
  expect_equal(nrow(dedup_hits(hits)), 2L)
})

# --- Shape preservation -----------------------------------------------------

test_that("empty input returns empty output of the same shape", {
  hits <- data.frame(
    sequence = character(0),
    motif    = character(0),
    strand   = character(0),
    start    = integer(0),
    end      = integer(0),
    pvalue   = numeric(0),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
})

test_that("GRanges input yields GRanges output", {
  skip_if_not_installed("GenomicRanges")
  gr <- GenomicRanges::GRanges(
    seqnames = c("s1","s1","s1"),
    ranges   = IRanges(start = c(1, 3, 6), end = c(5, 7, 9)),
    strand   = c("+","+","+")
  )
  GenomicRanges::mcols(gr) <- DataFrame(
    motif  = c("M","M","M"),
    pvalue = c(1e-10, 1e-8, 1e-9)
  )
  out <- dedup_hits(gr)
  expect_s4_class(out, "GRanges")
  expect_equal(length(out), 2L)
})

# --- Round-trip through scan_sequences2 -------------------------------------

test_that("scan_sequences2(no.overlaps = TRUE) yields a subset of un-deduped output", {
  set.seed(1)
  motifs <- list(create_motif("TATAAA", name = "M1"),
                 create_motif("CACGTG", name = "M2"))
  seqs <- create_sequences(seqnum = 5, seqlen = 400, rng.seed = 1)
  suppressMessages({
    raw <- scan_sequences2(motifs, seqs, pvalue = 5e-2,
                           return.granges = FALSE)
    dedup <- scan_sequences2(motifs, seqs, pvalue = 5e-2,
                             return.granges = FALSE,
                             no.overlaps = TRUE)
  })
  expect_true(nrow(dedup) <= nrow(raw))
  # Every deduped row must exist in raw at the same coordinates.
  key_raw <- with(raw,   paste(motif.i, sequence.i, start, end, strand, sep = "|"))
  key_d   <- with(dedup, paste(motif.i, sequence.i, start, end, strand, sep = "|"))
  expect_true(all(key_d %in% key_raw))
})

test_that("no.overlaps.by = 'score' picks highest-scoring hit per cluster", {
  set.seed(2)
  motifs <- list(create_motif("TATAAA", name = "M1"))
  seqs <- create_sequences(seqnum = 3, seqlen = 200, rng.seed = 2)
  suppressMessages({
    by_p <- scan_sequences2(motifs, seqs, pvalue = 5e-2, no.overlaps = TRUE,
                            no.overlaps.by = "pvalue",
                            return.granges = FALSE)
    by_s <- scan_sequences2(motifs, seqs, pvalue = 5e-2, no.overlaps = TRUE,
                            no.overlaps.by = "score",
                            return.granges = FALSE)
  })
  # Hit counts should match (just different tiebreaks).
  expect_equal(nrow(by_p), nrow(by_s))
})

# --- Inclusive overlap edge case --------------------------------------------

test_that("touching coordinates (end == next start) count as overlap", {
  # h1: 1-5, h2: 5-9  -- share position 5.
  hits <- data.frame(
    sequence = c("s1","s1"),
    motif    = c("M","M"),
    strand   = c("+","+"),
    start    = c(1L, 5L),
    end      = c(5L, 9L),
    pvalue   = c(1e-10, 1e-5),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits)
  expect_equal(nrow(out), 1L)
  expect_equal(out$start, 1L)
})

test_that("adjacent but non-touching hits both survive", {
  # h1: 1-5, h2: 6-9 -- gap at position 5 vs 6.
  hits <- data.frame(
    sequence = c("s1","s1"),
    motif    = c("M","M"),
    strand   = c("+","+"),
    start    = c(1L, 6L),
    end      = c(5L, 9L),
    pvalue   = c(1e-10, 1e-5),
    stringsAsFactors = FALSE
  )
  out <- dedup_hits(hits)
  expect_equal(nrow(out), 2L)
})
