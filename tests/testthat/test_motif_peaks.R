context("motif_peaks() / plot_motif_peaks()")

PEAK_COLS <- c("motif", "motif.i", "mode", "seq.length", "nhits",
               "best.window", "best.center", "hits.in", "hits.out",
               "expected.in", "enrichment", "log2.enrichment",
               "pvalue", "qvalue", "centers")

## Helper: plant a motif at a target centre position in n of the first n
## sequences of `seqs`. Position jitter +/- 4 bp.
plant_centred <- function(seqs, motif_str, n, target_center, rng.seed = 1) {
  set.seed(rng.seed)
  w  <- nchar(motif_str)
  ds <- as(motif_str, "DNAString")
  for (i in seq_len(n)) {
    pos <- target_center - w %/% 2L + sample(-4:4, 1)
    pos <- max(1L, min(length(seqs[[i]]) - w + 1L, pos))
    Biostrings::subseq(seqs[[i]], start = pos, width = w) <- ds
  }
  seqs
}

test_that("output shape and column types are stable", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 100, seqlen = 300, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 70, target_center = 150L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)

  r <- motif_peaks(hits, qvalue = 1)
  expect_true(is.data.frame(r))
  expect_equal(names(r), PEAK_COLS)
  expect_type(r$motif,        "character")
  expect_type(r$best.window,  "integer")
  expect_type(r$pvalue,       "double")
  expect_type(r$qvalue,       "double")
  expect_true(is.list(r$centers))
})

test_that("central mode finds a centrally-planted motif", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 200, seqlen = 500, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 150, target_center = 250L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)

  r <- motif_peaks(hits, qvalue = 1)
  expect_equal(nrow(r), 1L)
  expect_equal(r$best.center, 250)
  expect_lt(r$pvalue, 1e-20)
  expect_gt(r$enrichment, 5)
})

test_that("local mode finds an off-centre planted motif that central misses", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 200, seqlen = 500, rng.seed = 1)
  ## Plant off-centre at position 100
  seqs <- plant_centred(seqs, "TTGACATA", n = 150, target_center = 100L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)

  r_local   <- motif_peaks(hits, mode = "local",   qvalue = 1)
  r_central <- motif_peaks(hits, mode = "central", qvalue = 1)
  expect_equal(nrow(r_local), 1L)
  ## Local mode should find a centre near 100 with very low p.
  expect_lt(r_local$pvalue, 1e-20)
  expect_lt(abs(r_local$best.center - 100), 20)
  ## Central mode should produce a much larger p-value (the planted
  ## motif is nowhere near the centre).
  expect_gt(r_central$pvalue, r_local$pvalue * 1e10)
})

test_that("null fixture: random hits yield no enriched motifs", {
  set.seed(7)
  seqs <- create_sequences(seqnum = 50, seqlen = 300, rng.seed = 7)
  m <- create_motif("ACGT", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-2, return.granges = TRUE)
  ## Very strict qvalue should drop everything
  r <- motif_peaks(hits, qvalue = 1e-30)
  expect_equal(nrow(r), 0L)
  expect_equal(names(r), PEAK_COLS)
})

test_that("scan_sequences2 data.frame input is accepted (with explicit seq.length)", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 100, seqlen = 300, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 70, target_center = 150L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = FALSE)

  r <- motif_peaks(hits, seq.length = 300L, qvalue = 1)
  expect_equal(nrow(r), 1L)
})

test_that("scan_sequences2 GRanges path picks up seq.length automatically", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 100, seqlen = 300, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 70, target_center = 150L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)

  r <- motif_peaks(hits, qvalue = 1)
  expect_equal(r$seq.length, 300L)
})

test_that("scan_sequences (v1) GRanges path picks up seq.length automatically", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 100, seqlen = 300, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 70, target_center = 150L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- suppressWarnings(suppressMessages(
    scan_sequences(m, seqs, threshold = 0.001, threshold.type = "pvalue",
                   return.granges = TRUE, no.overlaps = FALSE)
  ))
  r <- motif_peaks(hits, qvalue = 1)
  expect_equal(nrow(r), 1L)
  expect_equal(r$seq.length, 300L)
})

test_that("scan_sequences (v1) data.frame input is accepted (uses `stop` column)", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 100, seqlen = 300, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 70, target_center = 150L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- suppressWarnings(suppressMessages(
    scan_sequences(m, seqs, threshold = 0.001, threshold.type = "pvalue",
                   return.granges = FALSE, no.overlaps = FALSE)
  ))
  ## v1 scan_sequences returns DataFrame; coerce to data.frame
  hits <- as.data.frame(hits)
  r <- motif_peaks(hits, seq.length = 300L, qvalue = 1)
  expect_equal(nrow(r), 1L)
})

test_that("data.frame input without seq.length errors informatively", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 30, seqlen = 200, rng.seed = 1)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-2, return.granges = FALSE)
  expect_error(motif_peaks(hits), regexp = "seq\\.length")
})

test_that("best-hit-per-sequence dedup: multiple hits per sequence count once", {
  set.seed(1)
  ## Same motif planted twice per sequence.
  seqs <- create_sequences(seqnum = 50, seqlen = 500, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 50, target_center = 240L,
                        rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 50, target_center = 260L,
                        rng.seed = 2)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)
  r <- motif_peaks(hits, qvalue = 1)
  ## nhits must be <= number of input sequences (50), not number of raw hits.
  expect_lte(r$nhits, 50L)
})

test_that("Bonferroni: wider window-search range yields larger p-value", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 200, seqlen = 500, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 150, target_center = 250L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)
  r_narrow <- motif_peaks(hits, qvalue = 1, window.step = 50L)
  r_wide   <- motif_peaks(hits, qvalue = 1, window.step = 5L)
  expect_gte(r_wide$pvalue, r_narrow$pvalue)
})

test_that("deterministic: same input twice -> identical output", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 100, seqlen = 300, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 60, target_center = 150L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)
  r1 <- motif_peaks(hits, qvalue = 1)
  r2 <- motif_peaks(hits, qvalue = 1)
  expect_equal(r1, r2)
})

test_that("plot_motif_peaks() returns a ggplot with the expected geoms", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 100, seqlen = 300, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 70, target_center = 150L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)
  r <- motif_peaks(hits, qvalue = 1)
  g <- plot_motif_peaks(r)
  expect_s3_class(g, "ggplot")
  geom_classes <- vapply(g$layers, function(l) class(l$geom)[1], character(1))
  expect_true(any(grepl("GeomRect", geom_classes)))
  expect_true(any(grepl("GeomBar|GeomHistogram", geom_classes)))
})

test_that("plot_motif_peaks() errors cleanly on a 0-row peaks input", {
  empty <- universalmotif:::empty_peaks_result()
  expect_error(plot_motif_peaks(empty), regexp = "0 rows")
})

## Plant a motif at a feature-frame-downstream position across mixed-strand
## features: '+' features get it at center + off, '-' features at the mirror
## (center - off) so that, after strand orientation, all align at center + off.
plant_oriented <- function(seqs, motif_str, n, strand_map, center, off,
                           rng.seed = 3) {
  set.seed(rng.seed)
  w  <- nchar(motif_str)
  ds <- as(motif_str, "DNAString")
  for (i in seq_len(n)) {
    pos <- if (strand_map[i] == "+") center + off else center - off
    pos <- pos + sample(-3:3, 1)
    pos <- max(1L, min(length(seqs[[i]]) - w + 1L, pos))
    Biostrings::subseq(seqs[[i]], start = pos, width = w) <- ds
  }
  seqs
}

test_that("seq.strand orients local mode (mirror peaks merge into one)", {
  set.seed(1)
  L <- 400L; n <- 200L
  seqs <- create_sequences(seqnum = n, seqlen = L, rng.seed = 1)
  smap <- rep(c("+", "-"), length.out = n)
  seqs <- plant_oriented(seqs, "TTGACATA", n = 150, strand_map = smap,
                         center = L %/% 2L, off = 60L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)

  sn <- GenomeInfoDb::seqlevels(hits)
  ss <- setNames(smap[seq_along(sn)], sn)

  r_no <- motif_peaks(hits, mode = "local", qvalue = 1)
  r_st <- motif_peaks(hits, mode = "local", seq.strand = ss, qvalue = 1)

  ## With orientation, the planted hits align downstream of centre and the
  ## peak is both more central-frame correct and far more significant.
  expect_lt(r_st$pvalue, r_no$pvalue)
  expect_gt(r_st$hits.in, r_no$hits.in)
  expect_lt(abs(r_st$best.center - (L %/% 2L + 60L)), 25)
})

test_that("seq.strand does not change central-mode p-values", {
  set.seed(1)
  L <- 400L; n <- 120L
  seqs <- create_sequences(seqnum = n, seqlen = L, rng.seed = 1)
  smap <- rep(c("+", "-"), length.out = n)
  seqs <- plant_oriented(seqs, "TTGACATA", n = 90, strand_map = smap,
                         center = L %/% 2L, off = 50L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)
  sn <- GenomeInfoDb::seqlevels(hits)
  ss <- setNames(smap[seq_along(sn)], sn)

  c0 <- motif_peaks(hits, mode = "central", qvalue = 1)
  c1 <- motif_peaks(hits, mode = "central", seq.strand = ss, qvalue = 1)
  expect_equal(c0$pvalue, c1$pvalue)
})

test_that("seq.strand errors on missing / invalid entries", {
  set.seed(1)
  seqs <- create_sequences(seqnum = 20, seqlen = 300, rng.seed = 1)
  seqs <- plant_centred(seqs, "TTGACATA", n = 15, target_center = 150L)
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)
  sn <- GenomeInfoDb::seqlevels(hits)
  ss <- setNames(rep("+", length(sn)), sn)

  expect_error(motif_peaks(hits, seq.strand = ss[-1]), regexp = "missing")
  bad <- ss; bad[1] <- "?"
  expect_error(motif_peaks(hits, seq.strand = bad), regexp = "\\+")
})
