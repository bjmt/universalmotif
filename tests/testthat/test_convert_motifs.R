context("Test motif conversion")
library(universalmotif)

test_that("motif conversion works", {

  m <- matrix(rep(0.25, 16), nrow = 4)
  rownames(m) <- c("A", "C", "G", "T")
  mm <- convert_motifs(m)

  expect_equal(mm@alphabet, "DNA")

})
