context("shuffle_sequences()")

test_that("sequence shuffling works", {

  seqs <- create_sequences()
  s <- shuffle_sequences(seqs)
  l <- shuffle_sequences(seqs, method = "linear", k = 2)
  m <- shuffle_sequences(seqs, method = "markov", k = 2)

  expect_s4_class(seqs, "DNAStringSet")
  expect_true(any(s != l))
  expect_true(any(s != m))
  expect_true(any(l != m))

})
