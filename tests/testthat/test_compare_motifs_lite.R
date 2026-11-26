context("compare_motifs_lite()")

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
  m <- compare_motifs_lite(motifs)
  expect_true(isSymmetric(m))
  expect_equal(unname(diag(m)), rep(1, length(motifs)),
               tolerance = 1e-9)
})

test_that("matrix-mode score is in [-1, 1]", {
  motifs <- make_motifs()
  m <- compare_motifs_lite(motifs)
  expect_true(all(m <= 1 + 1e-9))
  expect_true(all(m >= -1 - 1e-9))
})

test_that("self-comparison long-format hit has overlap == width, strand '+', score 1", {
  motifs <- make_motifs(n = 3)
  hits <- compare_motifs_lite(motifs, compare.to = 1, qvalue = 1)
  self <- hits[hits$subject == hits$target, ]
  expect_equal(nrow(self), 1L)
  expect_equal(self$overlap, ncol(motifs[[1]]@motif))
  expect_equal(self$strand, "+")
  expect_equal(self$score, 1, tolerance = 1e-9)
  expect_equal(self$offset, 0L)
})

test_that("long-format columns and order match the documented contract", {
  motifs <- make_motifs()
  hits <- compare_motifs_lite(motifs, compare.to = 1, qvalue = 1)
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
  ms <- compare_motifs_lite(motifs, matrix.out = "score")
  mp <- compare_motifs_lite(motifs, matrix.out = "pvalue")
  mq <- compare_motifs_lite(motifs, matrix.out = "qvalue")
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
  mq <- compare_motifs_lite(motifs, matrix.out = "qvalue")
  # Self-comparison should be the most significant entry of each row.
  expect_true(all(diag(mq) <= apply(mq, 1, min) + 1e-12))
})

# ---- Both null modes ------------------------------------------------------

test_that("score/offset/overlap/strand are identical between null modes", {
  motifs <- make_motifs()
  a_emp <- compare_motifs_lite(motifs, compare.to = 1, qvalue = 1, null = "empirical")
  a_par <- compare_motifs_lite(motifs, compare.to = 1, qvalue = 1, null = "parametric")
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
  mp <- compare_motifs_lite(motifs, null = "parametric", matrix.out = "pvalue")
  expect_true(isSymmetric(mp, tol = 1e-9))
})

test_that("custom bkg changes parametric p-values but not empirical p-values", {
  motifs <- make_motifs()
  uniform <- c(.25, .25, .25, .25)
  skewed  <- c(A = .4, C = .1, G = .1, T = .4)
  ## p-values for the same alignment in matrix form
  mp_uniform <- compare_motifs_lite(motifs, null = "parametric",
                                bkg = uniform, matrix.out = "pvalue")
  mp_skewed  <- compare_motifs_lite(motifs, null = "parametric",
                                bkg = skewed,  matrix.out = "pvalue")
  expect_false(isTRUE(all.equal(unname(mp_uniform), unname(mp_skewed))))
  ## empirical doesn't use bkg for the null, so should be unchanged
  me_uniform <- compare_motifs_lite(motifs, null = "empirical",
                                bkg = uniform, matrix.out = "pvalue")
  me_skewed  <- compare_motifs_lite(motifs, null = "empirical",
                                bkg = skewed,  matrix.out = "pvalue")
  expect_equal(me_uniform, me_skewed)
})

# ---- min.overlap ----------------------------------------------------------

test_that("increasing min.overlap can only reduce score-magnitude (best alignment is among a subset)", {
  motifs <- make_motifs(n = 4, width = 10)
  m3 <- compare_motifs_lite(motifs, min.overlap = 3, matrix.out = "score")
  m7 <- compare_motifs_lite(motifs, min.overlap = 7, matrix.out = "score")
  # diag is always 1 in both
  expect_true(all(abs(diag(m3) - 1) < 1e-9))
  expect_true(all(abs(diag(m7) - 1) < 1e-9))
  # off-diagonal cells: best from the larger set is at least as good
  expect_true(all(abs(m3) + 1e-9 >= abs(m7)))
})

# ---- qvalue filtering -----------------------------------------------------

test_that("qvalue = 0.01 is a strict subset of qvalue = 0.5", {
  motifs <- make_motifs(n = 8)
  loose <- compare_motifs_lite(motifs, compare.to = 1, qvalue = 0.5)
  tight <- compare_motifs_lite(motifs, compare.to = 1, qvalue = 0.01)
  expect_true(nrow(tight) <= nrow(loose))
  k_loose <- with(loose, paste(subject.i, target.i, sep = "|"))
  k_tight <- with(tight, paste(subject.i, target.i, sep = "|"))
  expect_true(all(k_tight %in% k_loose))
})

test_that("very stringent qvalue can return 0 rows with the right shape", {
  motifs <- make_motifs(n = 8)
  hits <- compare_motifs_lite(motifs, compare.to = 1, qvalue = 1e-30)
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

test_that("compare_motifs_lite() is deterministic", {
  motifs <- make_motifs()
  expect_equal(compare_motifs_lite(motifs), compare_motifs_lite(motifs))
  expect_equal(
    compare_motifs_lite(motifs, compare.to = 1, qvalue = 0.5),
    compare_motifs_lite(motifs, compare.to = 1, qvalue = 0.5)
  )
})

# ---- Error handling -------------------------------------------------------

test_that("AA motifs are rejected", {
  m <- create_motif("LLNN")
  expect_error(compare_motifs_lite(list(m, m)), "DNA/RNA")
})

test_that("invalid arguments raise informative errors", {
  motifs <- make_motifs()
  expect_error(compare_motifs_lite(motifs, qvalue = 0),     "in \\(0, 1\\]")
  expect_error(compare_motifs_lite(motifs, qvalue = 2),     "in \\(0, 1\\]")
  expect_error(compare_motifs_lite(motifs, min.overlap = 0),"positive integer")
  expect_error(compare_motifs_lite(motifs, RC = "yes"),     "single logical")
  expect_error(compare_motifs_lite(motifs, bkg = c(0.5, 0.5)), "length-4")
  expect_error(compare_motifs_lite(motifs, bkg = c(.3, .3, .3, .3)), "sum to 1")
})

test_that("compare.to as character resolves by motif name", {
  motifs <- make_motifs()
  h1 <- compare_motifs_lite(motifs, compare.to = "M2", qvalue = 1)
  h2 <- compare_motifs_lite(motifs, compare.to = 2,    qvalue = 1)
  expect_equal(h1, h2)
})

})  # end suppressMessages

test_that("compare_motifs_lite runs on RNA motifs", {
  m1 <- create_motif("ACGU", alphabet = "RNA", name = "rna1")
  m2 <- create_motif("ACGU", alphabet = "RNA", name = "rna2")
  res <- compare_motifs_lite(list(m1, m2))
  expect_true(isSymmetric(res))
  expect_equal(unname(diag(res)), c(1, 1), tolerance = 1e-6)
})

test_that("compare_motifs and compare_motifs_lite both run on a shared DNA fixture (v1/v2 parity smoke)", {
  ## v1 and v2 use different scoring stacks, so numerical parity is not
  ## expected; this is a smoke check that both produce a symmetric square
  ## matrix of the right size for the same input.
  motifs <- list(create_motif("TTGACATA", name = "a"),
                 create_motif("CTTGACAT", name = "b"),
                 create_motif("GGGCCCCC", name = "c"))
  m1 <- compare_motifs(motifs)
  m2 <- compare_motifs_lite(motifs)
  expect_equal(dim(m1), c(3L, 3L))
  expect_equal(dim(m2), c(3L, 3L))
  expect_true(isSymmetric(m1))
  expect_true(isSymmetric(m2))
})

test_that("compare_motifs_lite gives the same similarity matrix at nthreads = 1 and 2", {
  motifs <- list(create_motif("TTGACATA", name = "a"),
                 create_motif("CTTGACAT", name = "b"),
                 create_motif("GGGCCCCC", name = "c"))
  a <- compare_motifs_lite(motifs, nthreads = 1)
  b <- compare_motifs_lite(motifs, nthreads = 2)
  expect_equal(a, b)
})

test_that("consensus strings use U for RNA and follow the motif consensus rule", {
  r1 <- create_motif("ACGUUA", alphabet = "RNA", name = "r1")
  r2 <- create_motif("ACGUUAG", alphabet = "RNA", name = "r2")
  res <- compare_motifs_lite(c(r1, r2), compare.to = 1, qvalue = 1)
  expect_false(any(grepl("T", c(res$subject.consensus, res$target.consensus))))
  for (col in list(c(.4, .3, .2, .1), c(.6, .2, .1, .1))) {
    m <- create_motif(matrix(c(col, diag(4)), 4, dimnames = list(DNA_BASES, NULL)), name = "a")
    res <- compare_motifs_lite(c(m, m), compare.to = 1, qvalue = 1)
    expect_equal(res$subject.consensus[1], m@consensus)
  }
})

test_that("the parametric null gives near-perfect alignments a P-value above 0", {
  mat <- matrix(c(0, .015, .295, .69, .983, 0, .011, .005, .003, .883, .001, .113,
                  .181, .481, .069, .269, .021, .001, .737, .24, .008, .005, .542, .445,
                  .007, .021, .964, .008, .004, .001, .811, .183), 4,
                dimnames = list(DNA_BASES, NULL))
  m <- create_motif(sweep(mat, 2, colSums(mat), "/"), name = "a")
  p <- compare_motifs_lite(c(m, create_motif("ACGTACGT", name = "b")),
                           null = "parametric", matrix.out = "pvalue")
  expect_gt(p[1, 1], 0)
})

# ---- Complete scores, Euclidean distances and overhangs ----------------------
#
# The reference values below were worked out independently for the fixture
# motifs of helper-lite-fixture.R and are listed target by target.

expect_lite_rows <- function(res, offset, strand, overlap, score, pvalue,
                             tol = 1e-5) {
  res <- by_target(res)
  expect_equal(res$offset, offset)
  expect_equal(res$strand, strand)
  expect_equal(res$overlap, overlap)
  expect_equal(res$score, score, tolerance = tol)
  expect_equal(res$pvalue, pvalue, tolerance = tol)
}

test_that("the default comparison of the fixture agrees with the reference values", {
  ms <- lite_fixture()
  expect_lite_rows(compare_motifs_lite(ms, compare.to = 1, qvalue = 1),
    offset = c(0, -1, 1, 0, -3, 0, 1),
    strand = c("+", "+", "+", "-", "+", "-", "+"),
    overlap = c(9, 8, 8, 9, 5, 9, 8),
    score = c(1, 0.9367918, 0.7810601, 0.06043549, 0.6910652, 0.4062361, 0.66503),
    pvalue = c(4.43299e-12, 0.000225354, 0.0054934, 1, 0.0969893, 0.427816, 0.0213852))
  expect_equal(unname(compare_motifs_lite(ms)[6, ]),
               c(0.4062361, 0.4719982, 0.6103467, 0.2069896, 0.697356, 1, 0.7857823),
               tolerance = 1e-6)
})

test_that("the new arguments at their defaults leave results unchanged", {
  ms <- lite_fixture()
  for (out in c("score", "pvalue", "qvalue"))
    expect_identical(
      compare_motifs_lite(ms, matrix.out = out),
      compare_motifs_lite(ms, matrix.out = out, complete = FALSE, method = "PCC",
                          overhang = "ignore", use.gaps = TRUE))
  expect_identical(
    compare_motifs_lite(ms, compare.to = 2:3, qvalue = 1),
    compare_motifs_lite(ms, compare.to = 2:3, qvalue = 1, complete = FALSE,
                        method = "PCC", overhang = "ignore", use.gaps = TRUE))
})

test_that("complete scores agree with the reference values", {
  ms <- lite_fixture()
  expect_lite_rows(compare_motifs_lite(ms, compare.to = 1, qvalue = 1, complete = TRUE),
    offset = c(0, -1, 1, 0, 1, 0, 1),
    strand = c("+", "+", "+", "-", "-", "-", "+"),
    overlap = c(9, 8, 8, 9, 8, 9, 8),
    score = c(1, 0.9367918, 0.7810601, 0.06043549, 0.4038858, 0.4062361, 0.66503),
    pvalue = c(2.46277e-13, 0.000155432, 0.00176302, 0.995201, 0.240369, 0.15286,
               0.00599699))
  expect_lite_rows(compare_motifs_lite(ms, compare.to = 6, qvalue = 1, complete = TRUE,
                                       null = "parametric"),
    offset = c(0, -1, 1, 0, -3, 0, 1),
    strand = c("-", "-", "-", "-", "-", "+", "-"),
    overlap = c(9, 8, 8, 9, 5, 9, 8),
    score = c(0.4062361, 0.4719982, 0.6103467, 0.1292942, 0.697356, 1, 0.7857823),
    pvalue = c(0.113621, 0.122021, 0.0169857, 0.95813, 0.271192, 8.86235e-15,
               0.000226626))
  expect_lite_rows(compare_motifs_lite(ms, compare.to = 4, qvalue = 1, complete = TRUE,
                                       RC = FALSE),
    offset = c(0, 2, 1, 0, -3, -2, 1),
    strand = rep("+", 7),
    overlap = c(9, 7, 8, 9, 5, 7, 8),
    score = c(0.02237042, 0.141224, 0.05874193, 1, 0.459737, 0.2069896, 0.3007845),
    pvalue = c(0.961207, 0.930408, 0.963884, 6.56739e-14, 0.532366, 0.761735,
               0.305558))
})

test_that("Euclidean distances agree with the reference values", {
  ms <- lite_fixture()
  expect_lite_rows(compare_motifs_lite(ms, compare.to = 1, qvalue = 1, method = "EUCL"),
    offset = c(0, -1, 1, -3, -3, -2, 1),
    strand = rep("+", 7),
    overlap = c(9, 8, 8, 6, 5, 7, 8),
    score = c(0, 0.1953143, 0.2819487, 0.6635426, 0.3015772, 0.5168532, 0.3987325),
    pvalue = c(1.53923e-15, 0.000188241, 0.00469218, 1, 0.0481933, 0.772742,
               0.0634991))
  expect_lite_rows(compare_motifs_lite(ms, compare.to = 1, qvalue = 1, method = "EUCL",
                                       complete = TRUE, null = "parametric"),
    offset = c(0, -1, 1, -3, -3, 0, 1),
    strand = c("+", "+", "+", "+", "+", "-", "+"),
    overlap = c(9, 8, 8, 6, 5, 9, 8),
    score = c(0, 0.1953143, 0.2819487, 0.6635426, 0.3015772, 0.5624529, 0.3987325),
    pvalue = c(1.84632e-16, 2.42108e-07, 3.43428e-05, 0.999592, 0.0877776, 0.279622,
               0.00273695))
  d <- compare_motifs_lite(ms, method = "EUCL")
  expect_true(isSymmetric(d))
  expect_equal(unname(diag(d)), rep(0, 7))
  expect_true(all(d >= 0 & d <= sqrt(2)))
})

test_that("background overhangs agree with the reference values", {
  ms <- lite_fixture()
  expect_lite_rows(compare_motifs_lite(ms, compare.to = 1, qvalue = 1,
                                       overhang = "background"),
    offset = c(0, -1, 1, 0, -3, 0, 1),
    strand = c("+", "+", "+", "-", "+", "-", "+"),
    overlap = c(9, 8, 8, 9, 5, 9, 8),
    score = c(1, 0.8327038, 0.6942756, 0.06043549, 0.4065089, 0.4062361, 0.6259106),
    pvalue = c(4.43299e-12, 0.000217059, 0.0054934, 1, 0.0969893, 0.427816,
               0.0213852))
  expect_lite_rows(compare_motifs_lite(ms, compare.to = 1, qvalue = 1,
                                       overhang = "background", method = "EUCL",
                                       bkg = c(.3, .2, .2, .3)),
    offset = c(0, -1, 1, 4, -3, -2, 1),
    strand = c("+", "+", "+", "-", "+", "+", "+"),
    overlap = c(9, 8, 8, 5, 5, 7, 8),
    score = c(0, 0.2409087, 0.3086445, 0.6571249, 0.4423452, 0.5301697, 0.401942),
    pvalue = c(1.53923e-15, 0.000186881, 0.0046728, 1, 0.0760153, 0.642945,
               0.0634991))
})

test_that("background overhangs with correlations only lower the scores of partial overlaps", {
  ms <- lite_fixture()
  a <- by_target(compare_motifs_lite(ms, compare.to = 3, qvalue = 1))
  b <- by_target(compare_motifs_lite(ms, compare.to = 3, qvalue = 1,
                                     overhang = "background"))
  expect_equal(a$offset, b$offset)
  expect_equal(a$strand, b$strand)
  full <- a$overlap == 9 & vapply(ms, function(m) ncol(m@motif), integer(1)) == 9
  expect_equal(a$score[full], b$score[full])
  expect_true(all(b$score[!full] < a$score[!full]))
})

test_that("comparison options are checked", {
  ms <- lite_fixture(c("M1", "M2"))
  expect_error(compare_motifs_lite(ms, complete = TRUE, overhang = "background"),
               "cannot be combined")
  expect_error(compare_motifs_lite(ms, method = "ALLR"), "should be one of")
  expect_error(compare_motifs_lite(ms, overhang = "both"), "should be one of")
  expect_error(compare_motifs_lite(ms, complete = NA), "single logical")
  expect_error(compare_motifs_lite(ms, use.gaps = "yes"), "single logical")
})

test_that("complete scores depend on which motif is the query", {
  ms <- lite_fixture()
  p <- compare_motifs_lite(ms, complete = TRUE, matrix.out = "pvalue")
  expect_false(isSymmetric(p))
  expect_true(all(p >= 0 & p <= 1))
})

# ---- Gapped motifs ---------------------------------------------------------------

gapped_fixture <- function() {
  ms <- lite_fixture()
  ms[[6]] <- add_gap(ms[[6]], gaploc = 4, mingap = 0, maxgap = 3)
  ms
}

test_that("gapped comparisons agree with the reference values", {
  x <- gapped_fixture()
  res <- compare_motifs_lite(x, compare.to = 6, qvalue = 1)
  expect_lite_rows(res,
    offset = c(4, -1, 1, 1, -3, 0, 1),
    strand = c("+", "-", "-", "-", "-", "+", "-"),
    overlap = c(7, 8, 9, 9, 5, 9, 8),
    score = c(0.7354258, 0.5917948, 0.8273505, 0.3572281, 0.697356, 1, 0.7857823),
    pvalue = c(0.420072, 0.491505, 0.0167411, 1, 0.602345, 2.60559e-08, 0.0129109))
  res <- by_target(res)
  expect_equal(res$subject.gaps, c("2", "1", "3", "2", "0", "0", "0"))
  expect_equal(res$target.gaps, c("", "", "", "", "", "0", ""))
  expect_equal(res$subject.consensus[3], "GANNNGTCA")

  res <- compare_motifs_lite(x, compare.to = 1, qvalue = 1)
  expect_lite_rows(res,
    offset = c(0, -1, 1, 0, -3, -1, 1),
    strand = c("+", "+", "+", "-", "+", "-", "+"),
    overlap = c(9, 8, 8, 9, 5, 9, 8),
    score = c(1, 0.9367918, 0.7810601, 0.06043549, 0.6910652, 0.542895, 0.66503),
    pvalue = c(5.79021e-13, 0.000206615, 0.00456879, 1, 0.0929189, 0.732373,
               0.0218785))
  expect_equal(by_target(res)$target.gaps, c("", "", "", "", "", "1", ""))

  expect_lite_rows(compare_motifs_lite(x, compare.to = 6, qvalue = 1, complete = TRUE,
                                       method = "EUCL"),
    offset = c(0, -1, 1, 0, -3, 0, 1),
    strand = c("-", "-", "-", "-", "-", "+", "-"),
    overlap = c(9, 8, 9, 9, 5, 12, 8),
    score = c(0.5055545, 0.4427527, 0.2680973, 0.6160448, 0.420438, 0, 0.2942647),
    pvalue = c(0.631995, 0.625379, 0.0799677, 1, 1, 1.90362e-15, 0.00843349))

  y <- c(list(add_gap(lite_fixture("M1")[[1]], gaploc = c(3, 6), mingap = c(1, 0),
                      maxgap = c(2, 1)), x[[6]]), lite_fixture(c("M2", "M3", "M4")))
  res <- compare_motifs_lite(y, compare.to = 1, qvalue = 1)
  expect_lite_rows(res,
    offset = c(0, 0, -1, 3, 3),
    strand = c("+", "+", "+", "+", "-"),
    overlap = c(10, 9, 8, 8, 9),
    score = c(1, 0.6654105, 0.7382328, 0.7317124, 0.437684),
    pvalue = c(8.99908e-10, 0.69032, 0.138201, 0.233113, 1))
  expect_equal(by_target(res)$subject.gaps, c("1,0", "2,0", "2,0", "2,0", "2,1"))
  expect_equal(by_target(res)$target.gaps, c("1,0", "0", "", "", ""))
})

test_that("gapped comparison matrices agree with the reference values", {
  x <- gapped_fixture()
  s <- compare_motifs_lite(x)
  expect_equal(unname(s[6, ]), c(0.7354258, 0.5917948, 0.8273505, 0.3572281,
                                 0.697356, 1, 0.7857823), tolerance = 1e-6)
  expect_equal(unname(s[, 6]), c(0.542895, 0.5917948, 0.6103467, 0.3572281,
                                 0.697356, 1, 0.7857823), tolerance = 1e-6)
  q <- compare_motifs_lite(x, matrix.out = "qvalue")
  expect_equal(unname(q[6, ]), c(0.72281424, 0.72281424, 0.050223364, 1,
                                 0.72281424, 2.6055943e-08, 0.050223364),
               tolerance = 1e-5)
  expect_equal(unname(q[1, ]), c(5.7902096e-13, 1.2396920e-03, 1.3706371e-02, 1,
                                 1.3937840e-01, 8.7884706e-01, 4.3757085e-02),
               tolerance = 1e-5)
})

test_that("use.gaps = FALSE compares gapped motifs as if their parts followed one another", {
  x <- gapped_fixture()
  ms <- lite_fixture()
  a <- compare_motifs_lite(x, compare.to = 6, qvalue = 1, use.gaps = FALSE)
  b <- compare_motifs_lite(ms, compare.to = 6, qvalue = 1)
  expect_identical(a, b)
  expect_identical(compare_motifs_lite(x, use.gaps = FALSE), compare_motifs_lite(ms))
})

test_that("motifs without gaps give the same results either way", {
  ms <- lite_fixture()
  expect_identical(compare_motifs_lite(ms, compare.to = 1:2, qvalue = 1),
                   compare_motifs_lite(ms, compare.to = 1:2, qvalue = 1,
                                       use.gaps = FALSE))
  expect_false("subject.gaps" %in% colnames(
    compare_motifs_lite(ms, compare.to = 1, qvalue = 1)))
})

test_that("gapped motifs with too many combinations of gap sizes are refused", {
  ms <- lite_fixture()
  big <- add_gap(ms[[1]], gaploc = c(3, 6), mingap = c(0, 0), maxgap = c(5, 5))
  expect_error(compare_motifs_lite(c(list(big), ms[2:3])),
               "36 combinations of gap sizes, more than the 32")
  expect_equal(dim(compare_motifs_lite(c(list(big), ms[2:3]), use.gaps = FALSE)),
               c(3L, 3L))
  ok <- add_gap(ms[[1]], gaploc = c(3, 6), mingap = c(0, 0), maxgap = c(3, 7))
  expect_equal(dim(compare_motifs_lite(c(list(ok), ms[2:3]))), c(3L, 3L))
})
