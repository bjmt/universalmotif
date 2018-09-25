context("Test sequence shuffling")

test_that("sequence shuffling works", {

  seqs <- create_sequences()
  seqs <- shuffle_sequences(seqs)

  expect_s4_class(seqs, "DNAStringSet")

})
