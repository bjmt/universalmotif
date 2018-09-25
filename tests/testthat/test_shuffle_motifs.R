context("Test motif shuffling")

test_that("motif shuffling works", {

  m1 <- create_motif()
  expect_s4_class(shuffle_motifs(m1), "universalmotif")

})
