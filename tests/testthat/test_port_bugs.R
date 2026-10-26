# Bugs found while porting the package to Rust and Python.

test_that("get_matches() rounds scores up and handles single columns", {
  m <- create_motif(matrix(c(1, 0.75, -1, -2, 0, -1, -1, -1), 4,
    dimnames = list(DNA_BASES, NULL)), type = "PWM")
  expect_identical(get_matches(m, 0.7505), "AA")
  expect_identical(sort(get_matches(m, 0.75)), c("AA", "CA"))
  one <- suppressMessages(create_motif("A", nsites = 10, pseudocount = 1))
  expect_true("A" %in% get_matches(one, 0))
})

test_that("score_match() gives the scores scanning reports", {
  m <- create_motif("TGACGTCA", nsites = 20, pseudocount = 1)
  seqs <- Biostrings::DNAStringSet(c("ATGACGTCAGGTGACGTAATTGACCTCAT",
    "CCTGACGACATTGAAGTCACCATGACGTCA"))
  h <- scan_sequences(m, seqs, threshold = 0.3, threshold.type = "logodds",
    RC = TRUE, verbose = 0, calc.pvals = FALSE)
  expect_gt(nrow(h), 0)
  expect_identical(score_match(m, h$match), h$score)
})

test_that("normalize() gives a PCM without a pseudocount one", {
  pcm <- create_motif(matrix(c(10, 0, 0, 0, 0, 10, 0, 0), 4,
    dimnames = list(DNA_BASES, NULL)), type = "PCM", pseudocount = 0)
  n <- suppressMessages(normalize(pcm))
  expect_equal(n@pseudocount, 1)
  expect_identical(n@motif, pcm@motif)
})

test_that("an IC column of zeros converts back to a uniform column", {
  expect_equal(icm_to_ppm(c(0, 0, 0, 0)), rep(0.25, 4))
  u <- create_motif(matrix(c(0.25, 0.25, 0.25, 0.25, 0.7, 0.1, 0.1, 0.1), 4,
    dimnames = list(DNA_BASES, NULL)), type = "PPM")
  back <- convert_type(convert_type(u, "ICM"), "PPM")
  expect_equal(unname(back@motif[, 1]), rep(0.25, 4))
})

test_that("create_motif(bkg = NULL) uses the default background", {
  expect_s4_class(create_motif("ACGT", bkg = NULL), "universalmotif")
  expect_s4_class(create_motif(10, bkg = NULL), "universalmotif")
  expect_s4_class(create_motif(matrix(c(1, 0, 0, 0, 0, 1, 0, 0), 4,
    dimnames = list(DNA_BASES, NULL)), bkg = NULL), "universalmotif")
})

test_that("read_transfac() keeps counts from a single site", {
  f <- tempfile()
  writeLines(c("ID one", "P0\tA\tC\tG\tT", "01\t1\t0\t0\t0", "02\t0\t1\t0\t0",
    "XX", "//"), f)
  m <- read_transfac(f)
  expect_equal(unname(colSums(m@motif)), c(1, 1))
  expect_equal(m@nsites, 1)
})

test_that("negative letters mirror the positive ones", {
  f <- as.data.frame(fontDFroboto)[, c("x", "y", "order", "group")]
  f$group <- as.character(f$group)
  d <- make_position_polygon_data(c(A = 0.5, C = 0.2, G = -0.5, T = -0.2), f,
    sort.positions = TRUE, position.id = 1, flip.neg = TRUE)
  y <- split(d$y, d$group)
  # The largest letter is furthest from zero on both sides.
  expect_gt(max(y$A), max(y$C))
  expect_lt(min(y$G), min(y$T))
})
