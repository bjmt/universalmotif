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
