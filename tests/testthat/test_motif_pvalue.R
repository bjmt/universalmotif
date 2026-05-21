context("motif_pvalue()")

test_that("p-values from scores are ok", {

  m <- create_motif("SGDGNTGGAY", pseudocount = 1, nsites = 88)
  res <- motif_pvalue(m, 1)
  expect_equal(round(res, 7), 0.0009766)

})

test_that("scores from p-values are ok", {

  m <- create_motif("SGDGNTGGAY", pseudocount = 1, nsites = 88)
  res <- motif_pvalue(m, pvalue = 0.001, k = 12)
  expect_equal(round(res, 3), -0.037)

})

test_that("motif_pvalue() allow.nonfinite=TRUE doesn't error on score path (regression: sanitize_input rejected -Inf scores)", {

  m <- create_motif("ATCGTACGTG")
  s <- motif_score(m, 0, allow.nonfinite = TRUE)
  expect_true(is.numeric(s))
  pval <- motif_pvalue(m, score = s, allow.nonfinite = TRUE, method = "exhaustive")
  expect_true(is.numeric(pval))
  expect_true(pval >= 0 && pval <= 1)

})

test_that("motif_pvalue(method='exhaustive') agrees with 'dynamic' for motifs wider than 2*k (regression: wrong inner-loop index)", {

  # k=4, motif width 12 -> nsplit=3 -> triggers the mot_split.size() > 2 branch
  set.seed(1)
  m <- create_motif(create_sequences(seqlen = 12, seqnum = 50), pseudocount = 1)
  s <- 1.0
  p_dyn <- motif_pvalue(m, score = s, method = "dynamic")
  p_exh <- motif_pvalue(m, score = s, k = 4, method = "exhaustive")
  expect_equal(as.numeric(p_dyn), as.numeric(p_exh), tolerance = 0.05)

})

test_that("motif_{pvalue,score}_dynamic_batch_cpp matches the per-motif path bit-for-bit", {
  # The batched C++ entry point parallelises over motifs but must produce
  # identical numbers to the legacy `mapply(motif_pvalue_dynamic_single_cpp, ...)`
  # path. nthreads is also exercised.
  set.seed(1)
  mk <- function() {
    consensus <- paste(sample(c("A","C","G","T"), 8, replace = TRUE), collapse = "")
    convert_type(suppressMessages(normalize(create_motif(consensus))), "PWM")
  }
  motifs <- replicate(6, mk(), simplify = FALSE)
  mats   <- lapply(motifs, function(m) m@motif)
  bkgs   <- lapply(motifs, function(m) m@bkg[seq_len(nrow(m@motif))])
  scores <- lapply(seq_along(motifs), function(i) sort(runif(7, -3, 5)))
  pvals  <- lapply(seq_along(motifs), function(i)
                   sort(c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5)))

  old_p <- mapply(motif_pvalue_dynamic_single_cpp, mats, bkgs, scores,
                  SIMPLIFY = FALSE)
  new_p1 <- motif_pvalue_dynamic_batch_cpp(mats, bkgs, scores, nthreads = 1)
  new_p4 <- motif_pvalue_dynamic_batch_cpp(mats, bkgs, scores, nthreads = 4)
  expect_identical(old_p, new_p1)
  expect_identical(new_p1, new_p4)

  old_s <- mapply(motif_score_dynamic_single_cpp, mats, bkgs, pvals,
                  SIMPLIFY = FALSE)
  new_s1 <- motif_score_dynamic_batch_cpp(mats, bkgs, pvals, nthreads = 1)
  new_s4 <- motif_score_dynamic_batch_cpp(mats, bkgs, pvals, nthreads = 4)
  expect_identical(old_s, new_s1)
  expect_identical(new_s1, new_s4)
})

test_that("motif_pvalue() dynamic batched path agrees with per-motif loop on a real fixture", {
  # End-to-end smoke: the same answer should come back whether we call
  # motif_pvalue() with a list of motifs + list of score-vectors (batched
  # under the hood) or loop in R.
  set.seed(2)
  motifs <- lapply(seq_len(4), function(i)
    suppressMessages(normalize(create_motif(
      paste(sample(c("A","C","G","T"), 7, replace = TRUE), collapse = "")))))
  scores <- lapply(seq_along(motifs), function(i) sort(runif(5, -3, 6)))

  one_shot <- motif_pvalue(motifs, score = scores, method = "dynamic")
  loop     <- lapply(seq_along(motifs),
                     function(i) motif_pvalue(motifs[[i]], score = scores[[i]],
                                              method = "dynamic"))
  expect_equal(one_shot, loop)
})
