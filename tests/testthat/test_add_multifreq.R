context("add_multifreq()")

test_that("filling multifreq works", {

  m1 <- create_motif("AAAAA", nsites=10)

  seqs <- Biostrings::DNAStringSet(rep(c("AAAAA", "ATAAA"), 3))
  seqs2 <- Biostrings::DNAStringSet(rep(c("AAAAAC", "ATAAAC"), 3))

  m2 <- add_multifreq(m1, seqs, add.k = 2)
  m2.2 <- add_multifreq(m1, seqs2, add.k = 2)

  m2.multi <- m2["multifreq"]$`2`

  expect_equal(unname(m2.multi[1, ]), c(0.5, 0.5, 1, 1))
  expect_equal(m2@multifreq[[1]], m2.2@multifreq[[1]])

  expect_warning(add_multifreq(m2, seqs, add.k = 2))

  m3 <- create_motif("QQQQQ", alphabet = "QWERTY")
  seqs3 <- Biostrings::BStringSet(rep(c("QQQQQ", "QWQQQ"), 3))

  m4 <- add_multifreq(m3, seqs3, add.k = 2)

  expect_equal(unname(m4@multifreq[[1]][8, ]), c(0.5, 0.5, 1, 1))

  s <- create_sequences(seqlen = 10)
  expect_equal(add_multi(s, 2), add_multi_ANY(s, 2, Biostrings::DNA_BASES))

})
