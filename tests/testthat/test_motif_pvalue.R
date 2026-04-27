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
