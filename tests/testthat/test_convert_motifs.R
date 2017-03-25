library(universalmotif)
context("Test motif conversion tools")

test_motif <- create_motif("AVWWKTTG")

test_that("TFBSTools conversion works", {
  if (!requireNamespace("TFBSTools", quietly = TRUE)) skip("missing pkg")
  test_pfmatrix <- convert_motifs(test_motif, "TFBSTools-PFMatrix")
  test_pwmatrix <- convert_motifs(test_motif, "TFBSTools-PWMatrix")
  expect_identical(convert_motifs(test_pfmatrix)@consensus, test_motif@consensus)
  expect_equal(convert_motifs(test_pfmatrix)@icscore, test_motif@icscore)
  expect_identical(convert_motifs(test_pwmatrix)@consensus, test_motif@consensus)
  expect_equal(convert_motifs(test_pwmatrix)@icscore, test_motif@icscore)
})

test_that("seqLogo conversion works", {
  if (!requireNamespace("seqLogo", quietly = TRUE)) skip("missing pkg")
  test_pwm <- convert_motifs(test_motif, "seqLogo-pwm")
  expect_identical(convert_motifs(test_pwm)@consensus, test_motif@consensus)
  expect_equal(convert_motifs(test_pwm)@icscore, test_motif@icscore)
})

test_that("motifStack conversion works", {
  if (!requireNamespace("motifStack", quietly = TRUE)) skip("missing pkg")
  test_pcm <- convert_motifs(test_motif, "motifStack-pcm")
  test_pfm <- convert_motifs(test_motif, "motifStack-pfm")
  expect_identical(convert_motifs(test_pcm)@consensus, test_motif@consensus)
  expect_equal(convert_motifs(test_pcm)@icscore, test_motif@icscore)
  expect_identical(convert_motifs(test_pfm)@consensus, test_motif@consensus)
  expect_equal(convert_motifs(test_pfm)@icscore, test_motif@icscore)
})

test_that("PWMEnrich conversion works", {
  if (!requireNamespace("PWMEnrich", quietly = TRUE)) skip("missing pkg")
  test_pwm2 <- convert_motifs(test_motif, "PWMEnrich-PWM")
  expect_identical(convert_motifs(test_pwm2)@consensus, test_motif@consensus)
  expect_equal(convert_motifs(test_pwm2)@icscore, test_motif@icscore,
               tolerance = 0.1)
})

test_that("rGADEM conversion works", {
  if (!requireNamespace("rGADEM", quietly = TRUE)) skip("missing pkg")
  test_motif2 <- convert_motifs(test_motif, "rGADEM-motif")
  expect_identical(convert_motifs(test_motif2)@consensus, test_motif@consensus)
  expect_equal(convert_motifs(test_motif2)@icscore, test_motif@icscore)
})
