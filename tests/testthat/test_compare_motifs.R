context("Test motif comparison")
library(universalmotif)

test_that("basic comparison works", {

  motif1 <- create_motif("TTTWWW")
  motif2 <- create_motif("GGGWWW")

  res <- compare_motifs(list(motif1, motif2), progress = FALSE)

  expect_true(res[1, 1] == 1 && res[2, 2] == 1)
  expect_true(res[1, 2] == res[2, 1])
  expect_true(res[1, 2] > 0.3 && res[1, 2] < 0.4)

})
