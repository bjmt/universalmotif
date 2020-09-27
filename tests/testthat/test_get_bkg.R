context("get_bkg()")
library(Biostrings)

test_that("get_bkg() == oligonucleotideFrequency()", {

  seqs.DNA <- create_sequences()

  bkg.DNA <- get_bkg(seqs.DNA, k = 3, as.prob = FALSE)[["count"]]

  bkg.DNA2 <- oligonucleotideFrequency(seqs.DNA, 3, 1, as.prob = FALSE)
  bkg.DNA2 <- unname(colSums(bkg.DNA2))

  expect_equal(bkg.DNA, bkg.DNA2)

})
