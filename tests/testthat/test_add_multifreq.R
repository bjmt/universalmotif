context("add_multifreq()")

test_that("filling multifreq works", {

  m1 <- create_motif("AAAAA", nsites=10)
  seqs <- Biostrings::DNAStringSet(rep(c("AAAAA", "ATAAA"), 3))

  m2 <- add_multifreq(m1, seqs, add.k = 2)

  m2.multi <- m2["multifreq"]$`2`

  expect_equal(unname(m2.multi[1, ]), c(0.5, 0.5, 1, 1))

})
