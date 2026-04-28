context("motif_score()")

test_that("motif_score() returns a finite result with allow.nonfinite=TRUE and threshold.type='total' (regression: was stop()ing)", {

  m <- create_motif("ATCGTACGTG")
  s <- motif_score(m, 0, allow.nonfinite = TRUE)
  expect_true(is.numeric(s))
  expect_false(any(is.nan(s)))

})

test_that("motif_score() threshold.type='fromzero' still works with allow.nonfinite=TRUE", {

  m <- create_motif("ATCGTACGTG")
  s <- motif_score(m, c(0, 0.5, 1), allow.nonfinite = TRUE, threshold.type = "fromzero")
  expect_true(is.numeric(s))
  expect_equal(length(s), 3)
  expect_true(s[1] <= s[2] && s[2] <= s[3])

})
