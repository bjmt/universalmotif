context("compare_motifs2()")

suppressPackageStartupMessages({
  library(Biostrings)
})

# Reusable small fixture: six DNA motifs with deterministic seeds.
make_motifs <- function(seed = 1, n = 6, width = 8) {
  set.seed(seed)
  lapply(seq_len(n), function(i)
    create_motif(paste(sample(c("A","C","G","T"), width, replace = TRUE),
                       collapse = ""), name = paste0("M", i)))
}

suppressMessages({  ## quiet "Added a pseudocount" notes inside the loops

# ---- Hand-checked properties ----------------------------------------------

test_that("matrix mode is symmetric and has 1.0 on the diagonal", {
  motifs <- make_motifs()
  m <- compare_motifs2(motifs)
  expect_true(isSymmetric(m))
  expect_equal(unname(diag(m)), rep(1, length(motifs)),
               tolerance = 1e-9)
})

test_that("matrix-mode score is in [-1, 1]", {
  motifs <- make_motifs()
  m <- compare_motifs2(motifs)
  expect_true(all(m <= 1 + 1e-9))
  expect_true(all(m >= -1 - 1e-9))
})

test_that("self-comparison long-format hit has overlap == width, strand '+', score 1", {
  motifs <- make_motifs(n = 3)
  hits <- compare_motifs2(motifs, compare.to = 1, qvalue = 1)
  self <- hits[hits$subject == hits$target, ]
  expect_equal(nrow(self), 1L)
  expect_equal(self$overlap, ncol(motifs[[1]]@motif))
  expect_equal(self$strand, "+")
  expect_equal(self$score, 1, tolerance = 1e-9)
  expect_equal(self$offset, 0L)
})

test_that("long-format columns and order match the documented contract", {
  motifs <- make_motifs()
  hits <- compare_motifs2(motifs, compare.to = 1, qvalue = 1)
  expect_identical(
    colnames(hits),
    c("subject", "subject.i", "target", "target.i",
      "offset", "strand", "overlap",
      "score", "subject.consensus", "target.consensus",
      "pvalue", "qvalue")
  )
  expect_true(all(diff(hits$qvalue) >= -1e-12))  # sorted by qvalue asc
})

# ---- Matrix.out variants --------------------------------------------------

test_that("matrix.out variants return distinct matrices with same shape", {
  motifs <- make_motifs()
  ms <- compare_motifs2(motifs, matrix.out = "score")
  mp <- compare_motifs2(motifs, matrix.out = "pvalue")
  mq <- compare_motifs2(motifs, matrix.out = "qvalue")
  expect_equal(dim(ms), dim(mp))
  expect_equal(dim(mp), dim(mq))
  expect_equal(dimnames(ms), dimnames(mp))
  expect_equal(dimnames(mp), dimnames(mq))
  expect_false(identical(unname(ms), unname(mp)))
  expect_false(identical(unname(mp), unname(mq)))
  expect_true(all(mp >= 0 - 1e-12 & mp <= 1 + 1e-12))
  expect_true(all(mq >= 0 - 1e-12 & mq <= 1 + 1e-12))
})

test_that("matrix.out = qvalue diagonal is the self-comparison q-value", {
  motifs <- make_motifs()
  mq <- compare_motifs2(motifs, matrix.out = "qvalue")
  # Self-comparison should be the most significant entry of each row.
  expect_true(all(diag(mq) <= apply(mq, 1, min) + 1e-12))
})

# ---- Both null modes ------------------------------------------------------

test_that("score/offset/overlap/strand are identical between null modes", {
  motifs <- make_motifs()
  a_emp <- compare_motifs2(motifs, compare.to = 1, qvalue = 1, null = "empirical")
  a_par <- compare_motifs2(motifs, compare.to = 1, qvalue = 1, null = "parametric")
  # The two outputs may differ in row order due to q-value-based sorting,
  # so reorder both by (subject.i, target.i) before comparing.
  k <- function(x) order(x$subject.i, x$target.i)
  a_emp <- a_emp[k(a_emp), , drop = FALSE]
  a_par <- a_par[k(a_par), , drop = FALSE]
  for (col in c("score", "offset", "overlap", "strand",
                "subject.consensus", "target.consensus")) {
    expect_equal(a_emp[[col]], a_par[[col]], info = col)
  }
})

test_that("parametric mode is independent of which motif is the query (symmetric pvalue)", {
  motifs <- make_motifs()
  mp <- compare_motifs2(motifs, null = "parametric", matrix.out = "pvalue")
  expect_true(isSymmetric(mp, tol = 1e-9))
})

test_that("custom bkg changes parametric p-values but not empirical p-values", {
  motifs <- make_motifs()
  uniform <- c(.25, .25, .25, .25)
  skewed  <- c(A = .4, C = .1, G = .1, T = .4)
  ## p-values for the same alignment in matrix form
  mp_uniform <- compare_motifs2(motifs, null = "parametric",
                                bkg = uniform, matrix.out = "pvalue")
  mp_skewed  <- compare_motifs2(motifs, null = "parametric",
                                bkg = skewed,  matrix.out = "pvalue")
  expect_false(isTRUE(all.equal(unname(mp_uniform), unname(mp_skewed))))
  ## empirical doesn't use bkg for the null, so should be unchanged
  me_uniform <- compare_motifs2(motifs, null = "empirical",
                                bkg = uniform, matrix.out = "pvalue")
  me_skewed  <- compare_motifs2(motifs, null = "empirical",
                                bkg = skewed,  matrix.out = "pvalue")
  expect_equal(me_uniform, me_skewed)
})

# ---- min.overlap ----------------------------------------------------------

test_that("increasing min.overlap can only reduce score-magnitude (best alignment is among a subset)", {
  motifs <- make_motifs(n = 4, width = 10)
  m3 <- compare_motifs2(motifs, min.overlap = 3, matrix.out = "score")
  m7 <- compare_motifs2(motifs, min.overlap = 7, matrix.out = "score")
  # diag is always 1 in both
  expect_true(all(abs(diag(m3) - 1) < 1e-9))
  expect_true(all(abs(diag(m7) - 1) < 1e-9))
  # off-diagonal cells: best from the larger set is at least as good
  expect_true(all(abs(m3) + 1e-9 >= abs(m7)))
})

# ---- qvalue filtering -----------------------------------------------------

test_that("qvalue = 0.01 is a strict subset of qvalue = 0.5", {
  motifs <- make_motifs(n = 8)
  loose <- compare_motifs2(motifs, compare.to = 1, qvalue = 0.5)
  tight <- compare_motifs2(motifs, compare.to = 1, qvalue = 0.01)
  expect_true(nrow(tight) <= nrow(loose))
  k_loose <- with(loose, paste(subject.i, target.i, sep = "|"))
  k_tight <- with(tight, paste(subject.i, target.i, sep = "|"))
  expect_true(all(k_tight %in% k_loose))
})

test_that("very stringent qvalue can return 0 rows with the right shape", {
  motifs <- make_motifs(n = 8)
  hits <- compare_motifs2(motifs, compare.to = 1, qvalue = 1e-30)
  expect_s3_class(hits, "data.frame")
  expect_equal(nrow(hits), 0L)
  expect_identical(
    colnames(hits),
    c("subject", "subject.i", "target", "target.i",
      "offset", "strand", "overlap",
      "score", "subject.consensus", "target.consensus",
      "pvalue", "qvalue")
  )
})

# ---- Determinism ----------------------------------------------------------

test_that("compare_motifs2() is deterministic", {
  motifs <- make_motifs()
  expect_equal(compare_motifs2(motifs), compare_motifs2(motifs))
  expect_equal(
    compare_motifs2(motifs, compare.to = 1, qvalue = 0.5),
    compare_motifs2(motifs, compare.to = 1, qvalue = 0.5)
  )
})

# ---- Error handling -------------------------------------------------------

test_that("AA motifs are rejected", {
  m <- create_motif("LLNN")
  expect_error(compare_motifs2(list(m, m)), "DNA/RNA")
})

test_that("invalid arguments raise informative errors", {
  motifs <- make_motifs()
  expect_error(compare_motifs2(motifs, qvalue = 0),     "in \\(0, 1\\]")
  expect_error(compare_motifs2(motifs, qvalue = 2),     "in \\(0, 1\\]")
  expect_error(compare_motifs2(motifs, min.overlap = 0),"positive integer")
  expect_error(compare_motifs2(motifs, RC = "yes"),     "single logical")
  expect_error(compare_motifs2(motifs, bkg = c(0.5, 0.5)), "length-4")
  expect_error(compare_motifs2(motifs, bkg = c(.3, .3, .3, .3)), "sum to 1")
})

test_that("compare.to as character resolves by motif name", {
  motifs <- make_motifs()
  h1 <- compare_motifs2(motifs, compare.to = "M2", qvalue = 1)
  h2 <- compare_motifs2(motifs, compare.to = 2,    qvalue = 1)
  expect_equal(h1, h2)
})

})  # end suppressMessages
