context("scan_sequences()")

test_that("Results are accurate", {

  motif <- create_motif("AAAA", pseudocount = 1, nsites = 100)
  seq <- Biostrings::DNAStringSet("GGGAAAAGGGTTTTGGG")  # width 17
  # A four-base exact match has P = 1/256, so the default 1e-4
  # cutoff correctly yields no hits. Use an attainable cutoff here.
  res <- scan_sequences(motif, seq, RC = TRUE, threshold = 0.01, verbose = 0)

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
  r <- suppressWarnings(scan_sequences(motif, seq, threshold = 0.01,
                      threshold.type = "pvalue", calc.pvals = TRUE,
                      calc.qvals = TRUE, calc.qvals.method = "BH",
                      verbose = 0))
  expect_gt(nrow(r), 0L)
  if (nrow(r) > 0) {
    mLen <- ncol(motif@motif)
    mMax <- sum(Biostrings::width(seq) - mLen + 1)
    expected <- pmin(r$pvalue * mMax / rank(r$pvalue), 1)
    expect_equal(r$qvalue, expected, tolerance = 1e-12)
  }

})

test_that("scan_sequences() with multifreq motif on a short sequence gives no hits (regression: size_t underflow)", {

  m <- create_motif("ACGT", nsites = 100, pseudocount = 1)
  # Sequences the same length as the motif are used directly (no scan step)
  seqs_match <- Biostrings::DNAStringSet(rep("ACGT", 20))
  m <- add_multifreq(m, seqs_match, add.k = 2)
  short <- Biostrings::DNAStringSet("AC")
  res <- suppressWarnings(suppressMessages(scan_sequences(m, short, use.freq = 2,
      verbose = 0)))
  expect_equal(nrow(res), 0)

})

test_that("multifreq scan coordinates and matches cover the full k-let span", {

  train <- Biostrings::DNAStringSet(rep("ACGT", 20))
  m <- suppressMessages(create_motif(train, add.multifreq = 2,
                                     pseudocount = 1))
  seq <- Biostrings::DNAStringSet("TTACGTAA")
  hits <- suppressWarnings(scan_sequences(
    m, seq, use.freq = 2, RC = FALSE, threshold = -1e6,
    threshold.type = "logodds.abs", calc.pvals = FALSE, verbose = 0
  ))

  best <- which.max(hits$score)
  expect_equal(hits$start[best], 3L)
  expect_equal(hits$stop[best], 6L)
  expect_equal(hits$match[best], "ACGT")
  expect_true(all(hits$stop - hits$start + 1L == ncol(m@motif)))

})

test_that("all motifs must contain the requested multifreq matrix", {

  train <- Biostrings::DNAStringSet(rep("ACGTA", 20))
  m2 <- suppressMessages(create_motif(train, add.multifreq = 2,
                                      pseudocount = 1))
  m3 <- suppressMessages(create_motif(train, add.multifreq = 3,
                                      pseudocount = 1))
  expect_error(
    scan_sequences(list(m2, m3), Biostrings::DNAStringSet("ACGTACGTA"),
                   use.freq = 2, verbose = 0),
    "not all motifs have correct multifreqs", fixed = TRUE
  )

})

test_that("scan_sequences() suggests scan_sequences_lite() when arguments are compatible", {

  m <- create_motif("ACGTAC")
  seqs <- create_sequences(seqnum = 3, seqlen = 100, rng.seed = 1)

  old_opt <- getOption("universalmotif.suggest.scan_sequences_lite")
  on.exit(options(universalmotif.suggest.scan_sequences_lite = old_opt), add = TRUE)

  # Default call -- should fire.
  options(universalmotif.suggest.scan_sequences_lite = TRUE)
  expect_message(
    suppressWarnings(scan_sequences(m, seqs, threshold = 1e-2,
                                    threshold.type = "pvalue", RC = TRUE,
                                    calc.pvals = FALSE, calc.qvals = FALSE,
                                    verbose = 0)),
    "scan_sequences_lite()", fixed = TRUE
  )

  # Each disqualifying argument should suppress the hint.
  hint_silent <- function(...) {
    msgs <- capture_messages(
      suppressWarnings(scan_sequences(m, seqs, threshold = 1e-2,
                                      threshold.type = "pvalue", RC = TRUE,
                                      calc.pvals = FALSE, calc.qvals = FALSE,
                                      verbose = 0, ...))
    )
    expect_false(any(grepl("scan_sequences_lite", msgs, fixed = TRUE)))
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
  expect_false(any(grepl("scan_sequences_lite", msgs, fixed = TRUE)))

  # threshold.type other than pvalue: also suppress.
  msgs <- capture_messages(
    suppressWarnings(scan_sequences(m, seqs, threshold = 0.6,
                                    threshold.type = "logodds", RC = TRUE,
                                    calc.pvals = FALSE, calc.qvals = FALSE,
                                    verbose = 0))
  )
  expect_false(any(grepl("scan_sequences_lite", msgs, fixed = TRUE)))

  # The option silences the hint even when conditions are met.
  options(universalmotif.suggest.scan_sequences_lite = FALSE)
  msgs <- capture_messages(
    suppressWarnings(scan_sequences(m, seqs, threshold = 1e-2,
                                    threshold.type = "pvalue", RC = TRUE,
                                    calc.pvals = FALSE, calc.qvals = FALSE,
                                    verbose = 0))
  )
  expect_false(any(grepl("scan_sequences_lite", msgs, fixed = TRUE)))

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

test_that("gapped motifs scan on both strands, with zero-size gaps and after ungapped motifs", {
  seqs <- Biostrings::DNAStringSet(c(x = "TTACGTACGTTT", y = "TTACGTTTACGTTT"))
  g <- add_gap(create_motif("ACGTACGT", name = "g"), gaploc = 4, mingap = 1, maxgap = 2)
  res <- scan_sequences(g, seqs, threshold = 1e-3, RC = TRUE)
  expect_setequal(res$strand, c("+", "-"))
  expect_true(all(res$match == "ACGT..ACGT"))
  g0 <- add_gap(create_motif("ACGTACGT", name = "g0"), gaploc = 4, mingap = 0, maxgap = 1)
  res <- scan_sequences(g0, seqs, threshold = 1e-3)
  expect_equal(res$match, "ACGTACGT")
  res <- scan_sequences(list(create_motif("GGGG", name = "u"), g), seqs, threshold = 1e-3)
  expect_equal(unique(res$motif.i), 2L)
})

test_that("hits of every motif get their own P-values next to gapped motifs", {
  seqs <- Biostrings::DNAStringSet(c(x = "TTACGTACGTTTGATAAGTT", y = "TTACGTTTACGTTTAGATAAT"))
  g <- add_gap(create_motif("ACGTACGT", name = "g"), gaploc = 4, mingap = 1, maxgap = 2)
  u <- create_motif("GATAA", name = "u")
  both <- scan_sequences(list(g, u), seqs, threshold = 1e-2, RC = FALSE)
  alone_u <- scan_sequences(u, seqs, threshold = 1e-2, RC = FALSE)
  alone_g <- scan_sequences(g, seqs, threshold = 1e-2, RC = FALSE)
  expect_equal(both$pvalue[both$motif.i == 2], alone_u$pvalue)
  expect_equal(both$pvalue[both$motif.i == 1], alone_g$pvalue)
})

## gapped_fixture() and read_hits() are in helper-gapped.R.

scan_full_fixture <- function(...) {
  f <- gapped_fixture()
  suppressMessages(scan_sequences(f$motifs, f$seqs, RC = TRUE, warn.NA = FALSE,
                                  ...))
}

## A hit's identity, with '-' strand hits given by their leftmost position as
## scan_sequences_lite() gives them.
hit_key <- function(x) {
  x <- as.data.frame(x)
  stop <- if (is.null(x$stop)) x$end else x$stop
  paste(x$motif, x$sequence, pmin(x$start, stop), pmax(x$start, stop),
        x$strand, x$match)
}

test_that("gapped scans give the same hits as scanning every combination of gap sizes", {
  f <- gapped_fixture()
  ## A site of g2 in s5, and the same site with an unknown letter inside its
  ## gap in s4, which rules it out.
  seqs <- c(f$seqs, Biostrings::DNAStringSet(c(s4 = "CCAGGTNTCAGGACC",
                                               s5 = "CCAGGTATCAGGACC")))
  for (m in f$motifs[2:3]) {
    m <- convert_type(m, "PWM")
    got <- suppressWarnings(suppressMessages(scan_sequences(m, seqs,
      threshold = 2, threshold.type = "logodds.abs", RC = TRUE,
      calc.pvals = FALSE, correct.gaps = FALSE)))
    want <- lapply(ungap_single(m), function(x) {
      gaps <- get_gaplocs(x)
      x@gapinfo@isgapped <- FALSE
      x@name <- m@name
      res <- suppressWarnings(suppressMessages(scan_sequences(x, seqs,
        threshold = 2, threshold.type = "logodds.abs", RC = TRUE,
        calc.pvals = FALSE)))
      if (nrow(res) && length(gaps)) {
        res$match <- vapply(res$match, function(s) {
          s <- strsplit(s, "")[[1]]
          s[gaps] <- "."
          paste(s, collapse = "")
        }, character(1), USE.NAMES = FALSE)
      }
      as.data.frame(res)
    })
    want <- do.call(rbind, want)
    expect_gt(nrow(got), 0)
    expect_setequal(paste(hit_key(got), got$score),
                    paste(hit_key(want), want$score))
    if (m@name == "g2") {
      site <- function(s) any(got$sequence == s & got$start == 3 &
                              got$stop == 13 & got$strand == "+")
      expect_true(site("s5"))
      expect_false(site("s4"))
    }
  }
})

test_that("correct.gaps corrects gapped P-values and P-value thresholds", {
  combos <- c(u = 1, g = 3, g2 = 4)
  raw <- scan_full_fixture(threshold = 0.01, correct.gaps = FALSE)
  cor <- scan_full_fixture(threshold = 0.01)
  expect_lt(nrow(cor), nrow(raw))
  m <- match(hit_key(cor), hit_key(raw))
  expect_false(anyNA(m))
  expect_equal(cor$score, raw$score[m])
  expect_equal(cor$pvalue, unname(pmin(raw$pvalue[m] * combos[cor$motif], 1)))
  keep <- raw$pvalue * combos[raw$motif] <= 0.01
  expect_setequal(hit_key(cor), hit_key(raw)[keep])
  ## The score thresholds are those of the corrected cutoffs.
  f <- gapped_fixture()
  expect_equal(unique(cor$thresh.score[cor$motif == "g2"]),
               motif_pvalue(f$motifs[[3]], pvalue = 0.01 / 4))
  ## Motifs without gaps are not affected.
  expect_identical(as.data.frame(cor[cor$motif == "u", ]),
                   as.data.frame(raw[raw$motif == "u", ]))
})

test_that("correct.gaps leaves logodds thresholds alone but corrects P-values", {
  raw <- scan_full_fixture(threshold = 0.6, threshold.type = "logodds",
                           correct.gaps = FALSE)
  cor <- scan_full_fixture(threshold = 0.6, threshold.type = "logodds")
  expect_identical(hit_key(cor), hit_key(raw))
  combos <- c(u = 1, g = 3, g2 = 4)
  expect_equal(cor$pvalue, unname(pmin(raw$pvalue * combos[raw$motif], 1)))
})

test_that("scan_sequences() and scan_sequences_lite() agree on gapped motifs", {
  f <- gapped_fixture()
  for (correct in c(TRUE, FALSE)) {
    full <- scan_full_fixture(threshold = 0.01, correct.gaps = correct)
    lite <- suppressMessages(scan_sequences_lite(f$motifs, f$seqs,
      pvalue = 0.01, correct.gaps = correct, return.granges = FALSE))
    expect_setequal(paste(hit_key(full), full$score, signif(full$pvalue, 10)),
                    paste(hit_key(lite), lite$score, signif(lite$pvalue, 10)))
  }
})

test_that("gapped motifs are scanned on the strands respect.strand asks for", {
  f <- gapped_fixture()
  g <- f$motifs[[2]]
  g@strand <- "-"
  res <- suppressMessages(scan_sequences(list(f$motifs[[1]], g), f$seqs,
    threshold = 0.003, respect.strand = TRUE, warn.NA = FALSE))
  expect_identical(unique(res$strand[res$motif == "g"]), "-")
  both <- scan_full_fixture(threshold = 0.003)
  expect_setequal(hit_key(res[res$motif == "g", ]),
                  hit_key(both[both$motif == "g" & both$strand == "-", ]))
})

test_that("sequences shorter than a motif give no hits", {
  m <- create_motif("ACGTACGT", name = "m")
  s <- Biostrings::DNAStringSet(c(a = "ACGT", b = "TTACGTACGTTT"))
  expect_equal(unique(as.character(GenomicRanges::seqnames(scan_sequences_lite(m, s, pvalue = 1e-3)))), "b")
  expect_equal(unique(scan_sequences(m, s, threshold = 1e-3)$sequence), "b")
})
