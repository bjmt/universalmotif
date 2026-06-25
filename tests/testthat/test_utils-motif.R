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

## Per-utility positive/negative cases. Most of these have only been exercised
## indirectly via higher-level functions; the cases below pin down the leaves.

test_that("ppm_to_pwm produces finite values for a uniform background", {
  ppm <- c(0.25, 0.25, 0.25, 0.25)
  out <- ppm_to_pwm(ppm, bkg = c(0.25, 0.25, 0.25, 0.25),
                     pseudocount = 1, nsites = 100)
  expect_equal(length(out), 4L)
  expect_true(all(is.finite(out)))
})

test_that("ppm_to_pwm with pseudocount = 0 produces -Inf for zero probabilities", {
  ppm <- c(1, 0, 0, 0)
  out <- ppm_to_pwm(ppm, bkg = c(0.25, 0.25, 0.25, 0.25),
                     pseudocount = 0, nsites = 100)
  expect_true(any(is.infinite(out)))
})

test_that("ppm_to_pcm scales probabilities to integer-like counts", {
  out <- ppm_to_pcm(c(0.5, 0.25, 0.125, 0.125), nsites = 100)
  expect_equal(sum(out), 100, tolerance = 1e-6)
  expect_true(all(out >= 0))
})

test_that("pcm_to_ppm normalises a count vector", {
  out <- pcm_to_ppm(c(50, 25, 12, 13))
  expect_equal(sum(out), 1, tolerance = 1e-9)
  expect_equal(out[1], 0.5, tolerance = 1e-9)
})

test_that("position_icscore returns a non-negative score for a uniform column", {
  out <- position_icscore(c(0.25, 0.25, 0.25, 0.25),
                          bkg = c(0.25, 0.25, 0.25, 0.25),
                          type = "PPM")
  expect_equal(out, 0, tolerance = 1e-6)
})

test_that("position_icscore approaches 2 bits for a one-hot DNA column", {
  out <- position_icscore(c(1, 0, 0, 0),
                          bkg = c(0.25, 0.25, 0.25, 0.25),
                          type = "PPM")
  ## The internal pseudocount adjustment keeps the value just under 2.
  expect_true(out > 1.5 && out <= 2)
})

test_that("consensus_to_ppm round-trips a single letter", {
  out <- consensus_to_ppm("A")
  expect_equal(length(out), 4L)
  expect_equal(out[1], 0.997, tolerance = 0.01)
})

test_that("consensus_to_ppmAA produces a 20-element AA vector", {
  out <- consensus_to_ppmAA("A")
  expect_equal(length(out), 20L)
  expect_equal(sum(out), 1, tolerance = 1e-6)
})

test_that("get_consensus picks the dominant letter from a PPM column", {
  out <- get_consensus(c(0.9, 0.05, 0.025, 0.025), alphabet = "DNA",
                       type = "PPM")
  expect_equal(out, "A")
})

test_that("get_consensusAA returns an AA letter for an AA PPM column", {
  position <- rep(0, 20)
  position[1] <- 1
  out <- get_consensusAA(position, type = "PPM")
  expect_true(nchar(out) == 1L)
  expect_true(out %in% AA_STANDARD2)
})

test_that("get_consensus()/get_consensusAA() error cleanly on bad input (regression: crash on a universalmotif object)", {
  ## A whole motif object (or any non-numeric value) used to abort the R
  ## session during the S4 -> double conversion in C++; a too-short vector
  ## triggered out-of-bounds indexing. Both must now be catchable errors.
  m <- create_motif("CCNNAA")
  expect_error(get_consensus(m), "numeric vector of length 4")
  expect_error(get_consensus(m, type = "CWM"), "numeric vector of length 4")
  expect_error(get_consensus(c(0.5, 0.5)), "numeric vector of length 4")
  expect_error(get_consensusAA(m), "numeric vector of length 20")
  expect_error(get_consensusAA(c(0.5, 0.5)), "numeric vector of length 20")
})

test_that("add_gap toggles the gapinfo slot; ungap clears it", {
  m <- create_motif("ACGTACGT", nsites = 100)
  gapped <- add_gap(m, gaploc = 4L, mingap = 2L, maxgap = 2L)
  expect_true(gapped@gapinfo@isgapped)
  expect_equal(gapped@gapinfo@gaploc, 4L)
  expect_equal(gapped@gapinfo@mingap, 2L)
  expect_equal(gapped@gapinfo@maxgap, 2L)
  back <- ungap(gapped)
  expect_false(back@gapinfo@isgapped)
})

test_that("prob_match returns a finite probability for a matching string", {
  m <- create_motif("ACGT", pseudocount = 1, nsites = 100)
  p <- prob_match(m, "ACGT")
  expect_true(is.numeric(p))
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0))
})
