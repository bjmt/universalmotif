context("Test site sampling")

test_that("site sampling works", {

  m <- create_motif("AAAAA", nsites = 100)
  s <- sample_sites(m)

  expect_true(all(as.character(s) == "AAAAA"))

})
