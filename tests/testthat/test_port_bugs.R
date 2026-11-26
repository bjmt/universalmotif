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

test_that("read_motifs() skips unknown fields and write_motifs() keeps full precision", {
  m1 <- create_motif(matrix(c(12345.6, 1, 1, 1, 0.1, 2, 3, 12341.6), 4,
    dimnames = list(DNA_BASES, NULL)), type = "PCM")
  m2 <- convert_type(create_motif("ACGTRYN"), "PWM")
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  write_motifs(list(m1, m2), f)
  r <- read_motifs(f)
  expect_identical(r[[1]]@motif, m1@motif)
  expect_identical(r[[2]]@motif, m2@motif)
  lines <- readLines(f)
  lines <- append(lines, c("newfield: x", "another: 2"),
    after = grep("^name", lines)[1])
  writeLines(lines, f)
  expect_warning(r <- read_motifs(f), "newfield, another")
  expect_length(r, 2)
  expect_null(attr(r[[1]], "unknown"))
  expect_identical(r[[1]]@motif, m1@motif)
})

test_that("MEME, JASPAR and TRANSFAC files keep gaps in the alternative name", {
  m <- add_gap(create_motif("AGGTCAAGGTCA", nsites = 10, name = "nr"),
    gaploc = c(3, 6), mingap = c(0, 1), maxgap = c(1, 5))
  expect_identical(gap_notation(m), "AGG[0-1]TCA[1-5]AGGTCA")
  for (fmt in c("meme", "jaspar", "transfac")) {
    f <- tempfile()
    on.exit(unlink(f), add = TRUE)
    expect_silent(get(paste0("write_", fmt))(m, f))
    r <- suppressMessages(get(paste0("read_", fmt))(f))
    expect_identical(r@altname, "AGG[0-1]TCA[1-5]AGGTCA")
    expect_true(r@gapinfo@isgapped)
    expect_equal(r@gapinfo@gaploc, c(3, 6))
    expect_equal(r@gapinfo@mingap, c(0, 1))
    expect_equal(r@gapinfo@maxgap, c(1, 5))
  }
  # An alternative name is replaced with a warning.
  m@altname <- "orig"
  expect_warning(write_meme(m, tempfile()), "replaced by their gaps")
  # Notation that does not fit the motif is left as a name.
  expect_null(parse_gap_notation("AGG[1-2]TCA", 7))
  expect_null(parse_gap_notation("[1-2]AGG", 3))
  u <- create_motif("AGGTCA", altname = "AGG[1-2]TC")
  expect_false(with_noted_gaps(u)@gapinfo@isgapped)
  # A motif whose gaps were removed does not write a stale notation.
  u@altname <- "AGG[1-2]TCA"
  expect_length(file_altname(u), 0)
  # Custom alphabets have no consensus, so their gaps are still lost.
  cu <- add_gap(create_motif("ABCABC", alphabet = "ABC", nsites = 5),
    gaploc = 3, mingap = 1, maxgap = 2)
  expect_warning(write_meme(cu, tempfile()), "cannot preserve gap definitions")
})

test_that("lite merging and logos align motifs narrower than min.overlap", {
  a <- create_motif("TTGACGTCAAGG", nsites = 20, pseudocount = 0, name = "a")
  b <- create_motif("GACA", nsites = 20, name = "b")
  m <- merge_motifs_lite(c(a, b))
  expect_gt(max(abs(m@motif - convert_type(a, "PPM")@motif)), 0.4)
  b <- create_motif("GACGA", nsites = 20, name = "b")
  p <- view_motifs_lite(c(a, b), return.raw = TRUE)
  expect_equal(unname(which(colSums(p[[2]]) > 0)), 3:7)
})

test_that("compare_motifs_lite() leaves self-comparisons out of the q-values", {
  set.seed(2)
  q <- create_motif("TGACGTCAGG", nsites = 20, name = "q")
  hits <- lapply(1:4, function(i) create_motif(c("TGACGTCAGG", "TGACTTCAGG",
    "TGAGGTCAGC")[sample(3, 1)], nsites = 20, name = paste0("hit", i)))
  rand <- lapply(1:20, function(i) create_motif(10, nsites = 20,
    name = paste0("r", i)))
  r <- compare_motifs_lite(c(list(q), hits, rand), compare.to = 1, qvalue = 1)
  self <- r$target == "q"
  expect_identical(r$qvalue[self], r$pvalue[self])
  expect_equal(r$qvalue[!self], p.adjust(r$pvalue[!self], "BH"))
  m <- compare_motifs_lite(c(list(q), hits, rand), matrix.out = "qvalue")
  p <- compare_motifs_lite(c(list(q), hits, rand), matrix.out = "pvalue")
  expect_equal(diag(m), diag(p))
  expect_equal(unname(m[1, -1]), unname(p.adjust(p[1, -1], "BH")))
})
