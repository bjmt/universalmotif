context("trim_motifs()")

test_that("motif trimming works", {

  m1 <- create_motif("NNCCCNN", nsites = 100)
  m2 <- trim_motifs(m1)

  expect_equal(ncol(m2), 3)

})
