context("convert_motifs()")

test_that("motif conversion works", {

  m <- matrix(rep(0.25, 16), nrow = 4)
  rownames(m) <- c("A", "C", "G", "T")
  mm <- convert_motifs(m)

  expect_equal(mm@alphabet, "DNA")

})

test_that("list conversion to TFBSTools *MatrixList classes works", {

  skip_if_not_installed("TFBSTools")

  m <- create_motif(name = "m1")
  n <- create_motif(name = "m2")

  pfml <- convert_motifs(list(m, n), class = "TFBSTools-PFMatrixList")
  expect_s4_class(pfml, "PFMatrixList")
  expect_equal(length(pfml), 2)
  expect_equal(unname(vapply(pfml, function(x) x@name, character(1))),
               c("m1", "m2"))

  pwml <- convert_motifs(list(m, n), class = "TFBSTools-PWMatrixList")
  expect_s4_class(pwml, "PWMatrixList")
  expect_equal(length(pwml), 2)

  icml <- convert_motifs(list(m, n), class = "TFBSTools-ICMatrixList")
  expect_s4_class(icml, "ICMatrixList")
  expect_equal(length(icml), 2)

  pwml2 <- convert_motifs(as.list(pwml), class = "TFBSTools-PWMatrixList")
  expect_s4_class(pwml2, "PWMatrixList")
  expect_equal(length(pwml2), 2)

})
