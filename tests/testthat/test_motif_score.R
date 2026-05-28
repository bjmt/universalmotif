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

test_that("motif_score() rejects thresholds outside [0, 1]", {
  m <- create_motif("ACGT")
  expect_error(motif_score(m, threshold = 1.5), regexp = "between 0 and 1")
  expect_error(motif_score(m, threshold = -0.1), regexp = "between 0 and 1")
})

test_that("motif_score() at threshold = 0 equals the min possible score (with allow.nonfinite)", {
  m <- create_motif("ACGT")
  ## threshold = 0 means 0% of the [min, max] range -> the min score itself.
  s <- motif_score(m, threshold = 0, allow.nonfinite = TRUE, threshold.type = "total")
  expect_true(is.finite(s))
})

test_that("motif_score() at threshold = 1 equals the max possible score (with allow.nonfinite)", {
  m <- create_motif("ACGT")
  s <- motif_score(m, threshold = 1, allow.nonfinite = TRUE, threshold.type = "total")
  expect_true(is.finite(s))
})

test_that("motif_score() with allow.nonfinite = FALSE normalises -Inf PWM via pseudocount", {
  ## A one-hot PCM with pseudocount = 0 has -Inf entries in PWM space.
  ## With allow.nonfinite = FALSE, motif_score() silently normalises via
  ## pseudocount rather than propagating -Inf; the score must stay finite.
  m <- create_motif("ACGT", pseudocount = 0)
  s <- suppressMessages(motif_score(m, threshold = 0.5,
                                     allow.nonfinite = FALSE))
  expect_true(is.finite(s))
})
