context("utils-motif()")

test_that("get_matches() works when passed a length-1 list (regression: undefined `i`)", {

  m <- create_motif("ACGT", pseudocount = 1, nsites = 100)
  expect_error(get_matches(list(m), 0), regexp = NA)

})

test_that("get_matches() error message references the correct min score (regression: score.range[2] instead of [1])", {

  m <- create_motif("ACGT", pseudocount = 1, nsites = 100)
  rng <- motif_score(m, c(0, 1))
  min_score <- as.character(round(rng[[1]], 3))
  expect_error(get_matches(m, rng[[1]] - 100), regexp = min_score, fixed = TRUE)

})

test_that("score_match() returns -Inf, not NA, with allow.nonfinite = TRUE (regression: as.integer(-Inf) -> NA)", {

  m <- create_motif("ACGT", pseudocount = 0, nsites = 100)
  s <- score_match(m, "TTTT", allow.nonfinite = TRUE)
  expect_true(is.infinite(s) && s < 0)
  expect_false(is.na(s))

})
