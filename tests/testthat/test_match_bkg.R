context("match_bkg() / plot_match_bkg()")

## Helper: synthesize DNA sequences with a target per-sequence GC fraction
## and length. GC is approximate -- sampled iid per position.
make_seqs <- function(n, len, gc, seed = 1) {
  set.seed(seed)
  alph <- c("A", "C", "G", "T")
  prob <- c((1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2)
  ## When len is a single number, build n equal-length seqs; when a vector,
  ## per-sequence lengths.
  if (length(len) == 1L) len <- rep(len, n)
  out <- vapply(seq_len(n), function(i)
    paste0(sample(alph, len[i], replace = TRUE, prob = prob),
           collapse = ""),
    character(1))
  Biostrings::DNAStringSet(out)
}

## Helper: a per-sequence covariate vector, normally distributed around
## `mean`. Seeded so tests are deterministic.
make_cov <- function(n, mean, sd = 0.02, seed = 1) {
  set.seed(seed)
  stats::rnorm(n, mean, sd)
}

test_that("output shape, class, and length are correct", {
  set.seed(1)
  target   <- make_seqs(20, 200, 0.5, seed = 1)
  universe <- make_seqs(200, 200, 0.5, seed = 2)
  bkg <- match_bkg(target, universe, n.per.target = 1L)
  expect_s4_class(bkg, "DNAStringSet")
  expect_equal(length(bkg), length(target) * 1L)
})

test_that("RNAStringSet round-trips as RNAStringSet", {
  set.seed(1)
  target   <- Biostrings::RNAStringSet(c("ACGU", "ACGUACGU", "AAAACCGG"))
  universe <- Biostrings::RNAStringSet(c("UUUU", "AACCCCGG", "GGGGCCCC",
                                          "ACGUACGUAC", "AUAUAUAU",
                                          "GGCCGGCCGGCC", "AACCGGUU",
                                          "UCAGUCAG", "CAGUACGU",
                                          "UAUAUAUAUA", "GCGCGCGC",
                                          "ACACACAC"))
  bkg <- suppressWarnings(
    match_bkg(target, universe, n.per.target = 1L)
  )
  expect_s4_class(bkg, "RNAStringSet")
})

test_that("GC match quality: target GC ~ matched-bkg GC", {
  set.seed(1)
  ## Target: 50 seqs at GC ~0.7, length 200
  target <- make_seqs(50, 200, 0.7, seed = 11)
  ## Universe: 1000 seqs, GC uniform [0.1, 0.9], length 200
  universe <- do.call(c, lapply(seq(0.1, 0.9, by = 0.05), function(g)
    make_seqs(50, 200, g, seed = round(g * 100))))
  bkg <- match_bkg(target, universe, n.per.target = 1L)
  gc.t <- universalmotif:::compute_gc(target)
  gc.b <- universalmotif:::compute_gc(bkg)
  expect_lt(abs(mean(gc.t) - mean(gc.b)), 0.05)
})

test_that("length match quality: target length ~ matched-bkg length", {
  set.seed(1)
  target_lens <- rep(200, 40)
  target <- make_seqs(40, target_lens, 0.5, seed = 21)
  universe_lens <- sample(seq(50, 500, by = 10), 400, replace = TRUE)
  universe <- make_seqs(400, universe_lens, 0.5, seed = 22)
  bkg <- suppressWarnings(
    match_bkg(target, universe, n.per.target = 1L)
  )
  ## Most matches should fall within +/-50bp of the target length (one
  ## quantile bin worth on this fixture).
  diffs <- abs(width(bkg) - 200)
  expect_lt(mean(diffs), 100)
})

test_that("determinism: caller's set.seed -> identical output", {
  target   <- make_seqs(15, 200, 0.5, seed = 31)
  universe <- make_seqs(150, 200, 0.5, seed = 32)
  set.seed(42)
  b1 <- match_bkg(target, universe)
  set.seed(42)
  b2 <- match_bkg(target, universe)
  expect_equal(as.character(b1), as.character(b2))
})

test_that("unique = TRUE: no duplicate universe indices in indices output", {
  target   <- make_seqs(10, 200, 0.5, seed = 41)
  universe <- make_seqs(200, 200, 0.5, seed = 42)
  df <- match_bkg(target, universe, n.per.target = 3L, unique = TRUE,
                  return.indices = TRUE)
  expect_equal(anyDuplicated(df$universe.i), 0L)
})

test_that("return.indices = TRUE has the documented columns", {
  target   <- make_seqs(5, 200, 0.5, seed = 51)
  universe <- make_seqs(50, 200, 0.5, seed = 52)
  df <- match_bkg(target, universe, n.per.target = 2L,
                  return.indices = TRUE)
  expect_s3_class(df, "data.frame")
  expect_equal(names(df), c("target.i", "target.gc", "target.length",
                            "universe.i", "universe.gc", "universe.length"))
  expect_equal(nrow(df), 5L * 2L)
})

test_that("n.per.target = 5 returns 5N sequences", {
  target   <- make_seqs(4, 200, 0.5, seed = 61)
  universe <- make_seqs(100, 200, 0.5, seed = 62)
  bkg <- match_bkg(target, universe, n.per.target = 5L)
  expect_equal(length(bkg), 4L * 5L)
})

test_that("AA / non-DNA-RNA input is rejected", {
  aa <- Biostrings::AAStringSet(c("ACDEFG", "ACDEFGHI", "GHIKLM"))
  expect_error(match_bkg(aa, aa), regexp = "DNA/RNA")
})

test_that("DNA target + RNA universe is rejected", {
  d <- Biostrings::DNAStringSet(c("ACGT", "ACGTACGT"))
  r <- Biostrings::RNAStringSet(c("ACGU", "ACGUACGU", "AACC"))
  expect_error(match_bkg(d, r), regexp = "alphabet")
})

test_that("universe too small (unique=TRUE) warns and still returns", {
  target   <- make_seqs(10, 200, 0.5, seed = 71)
  universe <- make_seqs(5, 200, 0.5, seed = 72)
  expect_warning(
    bkg <- match_bkg(target, universe, n.per.target = 2L,
                     unique = TRUE),
    regexp = "smaller than"
  )
  expect_equal(length(bkg), 20L)
})

test_that("empty-bin fallback works (sparse universe + extreme target)", {
  ## Universe: 20 seqs with GC ~ 0.3, target: 5 seqs with GC ~ 0.8.
  ## Target's home GC bin will be empty in the universe.
  target   <- make_seqs(5,  200, 0.8, seed = 81)
  universe <- make_seqs(40, 200, 0.3, seed = 82)
  expect_warning(
    bkg <- match_bkg(target, universe),
    regexp = "ring-expansion"
  )
  expect_equal(length(bkg), 5L)
})

test_that("plot_match_bkg() returns a ggplot with density layers", {
  target   <- make_seqs(20, 200, 0.5, seed = 91)
  universe <- make_seqs(200, 200, 0.5, seed = 92)
  bkg <- match_bkg(target, universe)
  g <- plot_match_bkg(target, bkg)
  expect_s3_class(g, "ggplot")
  geom_classes <- vapply(g$layers, function(l) class(l$geom)[1], character(1))
  expect_true(any(grepl("GeomDensity", geom_classes)))
})

test_that("plot_match_bkg(by = 'gc') panels only the GC axis", {
  target   <- make_seqs(10, 200, 0.5, seed = 101)
  universe <- make_seqs(50, 200, 0.5, seed = 102)
  bkg <- match_bkg(target, universe)
  g <- plot_match_bkg(target, bkg, by = "gc")
  expect_s3_class(g, "ggplot")
})

## --- genome / exclude arg validation -----------------------------------------

test_that("supplying both `universe` and `genome` errors", {
  target   <- make_seqs(5, 100, 0.5, seed = 1)
  universe <- make_seqs(50, 100, 0.5, seed = 2)
  expect_error(
    match_bkg(target, universe = universe, genome = "fake"),
    regexp = "mutually exclusive"
  )
})

test_that("supplying neither `universe` nor `genome` errors", {
  target <- make_seqs(5, 100, 0.5, seed = 1)
  expect_error(match_bkg(target), regexp = "exactly one")
})

test_that("non-BSgenome `genome` argument errors", {
  target <- make_seqs(5, 100, 0.5, seed = 1)
  expect_error(match_bkg(target, genome = list(foo = 1)),
               regexp = "BSgenome")
})

test_that("genome shortcut works when BSgenome is available", {
  skip_if_not_installed("BSgenome.Athaliana.TAIR.TAIR9")
  suppressPackageStartupMessages({
    requireNamespace("BSgenome.Athaliana.TAIR.TAIR9")
    requireNamespace("GenomicRanges")
    requireNamespace("Biostrings")
  })
  ath <- BSgenome.Athaliana.TAIR.TAIR9::Athaliana
  ## 20 random 200-bp Arabidopsis windows as the "target"
  set.seed(2026)
  chr.lens <- GenomeInfoDb::seqlengths(ath)[paste0("Chr", 1:5)]
  chrs   <- sample(names(chr.lens), 20, replace = TRUE,
                   prob = chr.lens / sum(chr.lens))
  starts <- vapply(chrs,
                   function(ch) sample.int(chr.lens[ch] - 200, 1),
                   integer(1))
  target.gr <- GenomicRanges::GRanges(
    seqnames = chrs,
    ranges   = IRanges::IRanges(start = starts, width = 200)
  )
  target.seqs <- Biostrings::getSeq(ath, target.gr)

  set.seed(2026)
  bkg <- match_bkg(target.seqs, genome = ath, exclude = target.gr,
                   n.candidates = 5000L)
  expect_s4_class(bkg, "DNAStringSet")
  expect_equal(length(bkg), length(target.seqs))
  ## All widths should equal a target width (matched on length).
  expect_true(all(width(bkg) %in% width(target.seqs)))
  ## GC distributions should be close (within 0.1 abs mean diff).
  gc.t <- mean(Biostrings::letterFrequency(target.seqs, "GC", as.prob = TRUE))
  gc.b <- mean(Biostrings::letterFrequency(bkg,         "GC", as.prob = TRUE))
  expect_lt(abs(gc.t - gc.b), 0.1)
})

## --- covariate matching ------------------------------------------------------

test_that("covariates = NULL matches the no-covariate path exactly", {
  target   <- make_seqs(15, 200, 0.5, seed = 201)
  universe <- make_seqs(150, 200, 0.5, seed = 202)
  set.seed(7)
  a <- match_bkg(target, universe)
  set.seed(7)
  b <- match_bkg(target, universe, covariates = NULL,
                 universe.covariates = NULL)
  expect_identical(as.character(a), as.character(b))
})

test_that("single covariate: matched-bkg covariate ~ target covariate", {
  target   <- make_seqs(50,   200, 0.5, seed = 211)
  universe <- make_seqs(1000, 200, 0.5, seed = 212)
  sig.t <- make_cov(50, 0.7, sd = 0.02, seed = 213)
  set.seed(214); sig.u <- runif(1000)
  set.seed(215)
  df <- match_bkg(target, universe,
                  covariates          = data.frame(signal = sig.t),
                  universe.covariates = data.frame(signal = sig.u),
                  return.indices      = TRUE)
  expect_lt(abs(mean(df$target.signal) - mean(df$universe.signal)), 0.06)
})

test_that("two covariates: both axes matched within tolerance", {
  target   <- make_seqs(40,   200, 0.5, seed = 221)
  universe <- make_seqs(2000, 200, 0.5, seed = 222)
  a.t <- make_cov(40, 0.3, sd = 0.02, seed = 2231)
  b.t <- make_cov(40, 0.8, sd = 0.02, seed = 2232)
  set.seed(224); a.u <- runif(2000); b.u <- runif(2000)
  set.seed(225)
  df <- suppressWarnings(
    match_bkg(target, universe,
              covariates          = data.frame(aa = a.t, bb = b.t),
              universe.covariates = data.frame(aa = a.u, bb = b.u),
              n.bins.covariates   = 5L,
              return.indices      = TRUE)
  )
  expect_lt(abs(mean(df$target.aa) - mean(df$universe.aa)), 0.1)
  expect_lt(abs(mean(df$target.bb) - mean(df$universe.bb)), 0.1)
  expect_equal(nrow(df), 40L)
})

test_that("encode_bins reduces to the original two-axis index", {
  g <- c(1L, 3L, 5L, 2L); l <- c(2L, 1L, 4L, 4L)
  ng <- 5L; nl <- 4L
  expect_equal(universalmotif:::encode_bins(cbind(g, l), c(ng, nl)),
               (g - 1L) * nl + l)
  ## A single coordinate vector encodes to a length-1 key.
  expect_equal(universalmotif:::encode_bins(c(3L, 4L), c(ng, nl)),
               (3L - 1L) * nl + 4L)
})

test_that("ring_offsets: magnitudes sum to ring; 2-D matches Manhattan ring", {
  for (r in 0:3) {
    offs <- universalmotif:::ring_offsets(r, 2L)
    sums <- vapply(offs, function(d) sum(abs(d)), integer(1))
    expect_true(all(sums == r))
    expect_equal(length(offs), if (r == 0L) 1L else 4L * r)
    keys <- vapply(offs, paste, character(1), collapse = ",")
    expect_equal(anyDuplicated(keys), 0L)            # no duplicate offsets
  }
  ## 3-D ring sizes are 4r^2 + 2 for r >= 1 (1, 6, 18 for r = 0, 1, 2).
  expect_equal(length(universalmotif:::ring_offsets(0L, 3L)), 1L)
  expect_equal(length(universalmotif:::ring_offsets(1L, 3L)), 6L)
  expect_equal(length(universalmotif:::ring_offsets(2L, 3L)), 18L)
})

test_that("return.indices gains target.<cov>/universe.<cov> columns", {
  target   <- make_seqs(5,  200, 0.5, seed = 231)
  universe <- make_seqs(60, 200, 0.5, seed = 232)
  set.seed(233); s.t <- runif(5); s.u <- runif(60)
  df <- match_bkg(target, universe, n.per.target = 2L,
                  covariates          = data.frame(signal = s.t),
                  universe.covariates = data.frame(signal = s.u),
                  return.indices      = TRUE)
  expect_equal(names(df), c("target.i", "target.gc", "target.length",
                            "universe.i", "universe.gc", "universe.length",
                            "target.signal", "universe.signal"))
  expect_equal(nrow(df), 5L * 2L)
  expect_true(all(df$target.signal    %in% s.t))
  expect_true(all(df$universe.signal  %in% s.u))
})

test_that("covariate validation errors", {
  target   <- make_seqs(5,  200, 0.5, seed = 241)
  universe <- make_seqs(40, 200, 0.5, seed = 242)
  set.seed(243)
  ct <- data.frame(signal = runif(5))
  cu <- data.frame(signal = runif(40))

  expect_error(match_bkg(target, genome = "x", covariates = ct),
               regexp = "cannot be combined with `genome`")
  expect_error(match_bkg(target, universe, universe.covariates = cu),
               regexp = "without `covariates`")
  expect_error(match_bkg(target, universe, covariates = ct),
               regexp = "requires a matching")
  expect_error(match_bkg(target, universe,
                         covariates          = data.frame(signal = runif(4)),
                         universe.covariates = cu),
               regexp = "row")
  expect_error(match_bkg(target, universe, covariates = ct,
                         universe.covariates = data.frame(signal = runif(39))),
               regexp = "row")
  expect_error(match_bkg(target, universe, covariates = ct,
                         universe.covariates = data.frame(other = runif(40))),
               regexp = "identical column")
  expect_error(match_bkg(target, universe,
                         covariates          = data.frame(signal = letters[1:5]),
                         universe.covariates = cu),
               regexp = "must be numeric")
  ct.na <- ct; ct.na$signal[1] <- NA
  expect_error(match_bkg(target, universe, covariates = ct.na,
                         universe.covariates = cu),
               regexp = "NA/NaN/Inf")
  expect_error(match_bkg(target, universe,
                         covariates          = data.frame(gc = runif(5)),
                         universe.covariates = data.frame(gc = runif(40))),
               regexp = "reserved")
  expect_error(match_bkg(target, universe, covariates = ct,
                         universe.covariates = cu, n.bins.covariates = 1L),
               regexp = ">= 2")
})

test_that("constant covariate emits a message and still returns", {
  target   <- make_seqs(10,  200, 0.5, seed = 251)
  universe <- make_seqs(100, 200, 0.5, seed = 252)
  ct <- data.frame(signal = rep(0.5, 10))
  cu <- data.frame(signal = rep(0.5, 100))
  expect_message(
    bkg <- match_bkg(target, universe, covariates = ct,
                     universe.covariates = cu),
    regexp = "constant"
  )
  expect_equal(length(bkg), 10L)
})

test_that("bin.type.covariates equalwidth differs from quantile on skewed data", {
  target   <- make_seqs(20,  200, 0.5, seed = 261)
  universe <- make_seqs(400, 200, 0.5, seed = 262)
  set.seed(263); s.t <- rexp(20); s.u <- rexp(400)
  ct <- data.frame(signal = s.t); cu <- data.frame(signal = s.u)
  set.seed(264)
  q <- suppressWarnings(
    match_bkg(target, universe, covariates = ct, universe.covariates = cu,
              bin.type.covariates = "quantile", return.indices = TRUE))
  set.seed(264)
  e <- suppressWarnings(
    match_bkg(target, universe, covariates = ct, universe.covariates = cu,
              bin.type.covariates = "equalwidth", return.indices = TRUE))
  expect_false(identical(q$universe.i, e$universe.i))
})

test_that("determinism with covariates: same seed -> identical output", {
  target   <- make_seqs(12,  200, 0.5, seed = 271)
  universe <- make_seqs(120, 200, 0.5, seed = 272)
  set.seed(273); s.t <- runif(12); s.u <- runif(120)
  ct <- data.frame(signal = s.t); cu <- data.frame(signal = s.u)
  set.seed(9)
  b1 <- suppressWarnings(
    match_bkg(target, universe, covariates = ct, universe.covariates = cu))
  set.seed(9)
  b2 <- suppressWarnings(
    match_bkg(target, universe, covariates = ct, universe.covariates = cu))
  expect_equal(as.character(b1), as.character(b2))
})

test_that("3 covariates: sparse universe triggers ring-expansion fallback", {
  target   <- make_seqs(8,  200, 0.5, seed = 281)
  universe <- make_seqs(40, 200, 0.5, seed = 282)
  set.seed(283)
  ct <- data.frame(a = runif(8),  b = runif(8),  c = runif(8))
  cu <- data.frame(a = runif(40), b = runif(40), c = runif(40))
  expect_warning(
    bkg <- match_bkg(target, universe, covariates = ct,
                     universe.covariates = cu, n.bins.covariates = 15L),
    regexp = "ring-expansion"
  )
  expect_equal(length(bkg), 8L)
})

test_that("plot_match_bkg(indices=) facets gc, length and covariates", {
  target   <- make_seqs(20,  200, 0.5, seed = 291)
  universe <- make_seqs(200, 200, 0.5, seed = 292)
  set.seed(293); s.t <- runif(20); s.u <- runif(200)
  df <- match_bkg(target, universe,
                  covariates          = data.frame(signal = s.t),
                  universe.covariates = data.frame(signal = s.u),
                  return.indices      = TRUE)
  g <- plot_match_bkg(indices = df)
  expect_s3_class(g, "ggplot")
  expect_setequal(levels(g$data$axis),
                  c("GC fraction", "Sequence length", "signal"))
  ## `by` can subset to a single covariate panel.
  g2 <- plot_match_bkg(indices = df, by = "signal")
  expect_setequal(levels(g2$data$axis), "signal")
})
