context("Test motif alphabet switching")

test_that("motif alphabet switching works", {

  m1 <- create_motif()
  m2 <- switch_alph(m1)

  expect_s4_class(m2, "universalmotif")
  expect_equal(m2@alphabet, "RNA")

})
