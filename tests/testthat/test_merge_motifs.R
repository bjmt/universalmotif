context("Test motif merging")
library(universalmotif)

test_that("motif merging works", {

  m1 <- create_motif("AAAAAAA", nsites = 10)
  m2 <- create_motif("TTTTTTT", nsites = 10)

  m3 <- merge_motifs(list(m1, m2), tryRC = FALSE, min.overlap = 10)

  expect_equal(m3@consensus, "WWWWWWW")

})
