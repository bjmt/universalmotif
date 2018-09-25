context("Test motif alphabet switching")

test_that("motif alphabet switching works", {

  m1 <- create_motif()
  m2 <- shuffle_motifs(m)

  expect_s4_class(m2, "universalmotif")
  expect_true(m1@consensus != m2@consensus)

})
