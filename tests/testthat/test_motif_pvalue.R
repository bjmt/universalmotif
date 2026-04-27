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
