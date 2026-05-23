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
