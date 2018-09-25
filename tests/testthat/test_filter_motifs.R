context("Test motif filtering")

test_that("motif filtering works", {

  m1 <- create_motif("AAAAAAA", name = "motif1")
  m2 <- create_motif("NNNNN", name = "motif2")
  m <- list(m1, m2)

  expect_equal(length(filter_motifs(m, name = "motif1")), 1)
  expect_equal(length(filter_motifs(m, width = 6)), 1)
  expect_equal(length(filter_motifs(m, icscore = 12)), 1)

})
