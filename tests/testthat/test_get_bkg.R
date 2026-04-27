context("get_bkg()")
library(Biostrings)

test_that("get_bkg() == oligonucleotideFrequency()", {

  seqs.DNA <- create_sequences()

  bkg.DNA <- get_bkg(seqs.DNA, k = 3)[["count"]]

  bkg.DNA2 <- oligonucleotideFrequency(seqs.DNA, 3, 1, as.prob = FALSE)
  bkg.DNA2 <- unname(colSums(bkg.DNA2))

  expect_equal(bkg.DNA, bkg.DNA2)

})

test_that("get_bkg() rejects fractional window.size that resolves to 0 (regression)", {

  short_seq <- Biostrings::DNAStringSet("ACGT")
  expect_error(get_bkg(short_seq, k = 1, window = TRUE, window.size = 0.001))

})
