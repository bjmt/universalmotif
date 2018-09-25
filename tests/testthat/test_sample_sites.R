context("Test site sampling")
library(universalmotif)

test_that("site sampling works", {

  m <- create_motif("AAAAA", nsites = 100)
  s <- sample_sites(m)

  expect_true(all(as.character(s) == "AAAAA"))

})
