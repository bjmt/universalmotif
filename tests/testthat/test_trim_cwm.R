context("trim_cwm()")

flanked_cwm <- function() {
  weak  <- matrix(0.02, nrow = 4, ncol = 3,
                  dimnames = list(c("A","C","G","T"), NULL))
  strong <- matrix(c(0.9, -0.3, -0.3, -0.3,
                    -0.3,  0.9, -0.3, -0.3,
                    -0.3, -0.3,  0.9, -0.3,
                    -0.3, -0.3, -0.3,  0.9),
                  nrow = 4, byrow = FALSE,
                  dimnames = list(c("A","C","G","T"), NULL))
  mat <- cbind(weak, strong, weak)
  create_motif(mat, type = "CWM", name = "flanked")
}

test_that("default trim drops weak edge columns and preserves CWM type", {
  cwm <- flanked_cwm()
  before_w <- ncol(cwm["motif"])
  trimmed <- trim_cwm(cwm)
  expect_equal(trimmed@type, "CWM")
  expect_lt(ncol(trimmed@motif), before_w)
  ## The strong core (4 columns) should survive.
  expect_gte(ncol(trimmed@motif), 4L)
})

test_that("abs.threshold takes precedence and uses an absolute cutoff", {
  cwm <- flanked_cwm()
  ## Set abs.threshold above the weak flanks (sum abs = 0.08) but
  ## below the strong core (sum abs = 1.8). Only the strong core
  ## should remain.
  trimmed <- trim_cwm(cwm, abs.threshold = 0.5)
  expect_equal(trimmed@type, "CWM")
  expect_equal(ncol(trimmed@motif), 4L)
})

test_that("trim.from = 'left' only trims the left edge", {
  cwm <- flanked_cwm()
  trimmed <- trim_cwm(cwm, trim.from = "left")
  ## Weak right flank stays; weak left flank goes.
  expect_lt(ncol(trimmed@motif), ncol(cwm@motif))
  expect_gt(ncol(trimmed@motif), 4L)
})

test_that("trim_cwm errors on non-CWM input", {
  cwm <- flanked_cwm()
  ppm <- convert_type(cwm, "PPM")
  expect_error(trim_cwm(ppm), "requires CWM motifs")
})

test_that("trim_cwm on a uniformly-strong CWM returns the input unchanged", {
  strong <- matrix(c(0.9, -0.3, -0.3, -0.3,
                    -0.3,  0.9, -0.3, -0.3,
                    -0.3, -0.3,  0.9, -0.3,
                    -0.3, -0.3, -0.3,  0.9),
                  nrow = 4, byrow = FALSE,
                  dimnames = list(c("A","C","G","T"), NULL))
  cwm <- create_motif(strong, type = "CWM", name = "all-strong")
  trimmed <- trim_cwm(cwm)
  expect_equal(ncol(trimmed@motif), ncol(cwm@motif))
})
