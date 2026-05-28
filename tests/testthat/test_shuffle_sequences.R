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

test_that("shuffle_sequences() rejects window.size < k (regression)", {

  seqs <- Biostrings::DNAStringSet("ACGTACGTACGT")
  expect_error(shuffle_sequences(seqs, window = TRUE, window.size = 2, k = 3))

})

test_that("k = 1 routes to the dedicated 1-let path for every method", {
  set.seed(1)
  seqs <- Biostrings::DNAStringSet(paste(sample(c("A","C","G","T"), 200,
                                                  replace = TRUE),
                                          collapse = ""))
  for (m in c("euler", "linear", "markov")) {
    s <- shuffle_sequences(seqs, method = m, k = 1)
    expect_s4_class(s, "DNAStringSet")
    expect_equal(width(s), width(seqs))
    ## k = 1 preserves the mononucleotide composition exactly.
    orig <- table(strsplit(as.character(seqs)[1], "")[[1]])
    new  <- table(strsplit(as.character(s)[1], "")[[1]])
    expect_equal(sort(names(orig)), sort(names(new)))
    expect_equal(as.vector(orig[names(new)]), as.vector(new[names(new)]))
  }
})

test_that("k larger than sequence length is rejected", {
  seqs <- Biostrings::DNAStringSet("ACGT")
  expect_error(shuffle_sequences(seqs, method = "euler", k = 10),
               regexp = "shortest sequence length")
  expect_error(shuffle_sequences(seqs, method = "markov", k = 10),
               regexp = "shortest sequence length")
  expect_error(shuffle_sequences(seqs, method = "linear", k = 10),
               regexp = "shortest sequence length")
})

test_that("seeded shuffle is reproducible at nthreads = 2", {
  set.seed(1)
  seqs <- Biostrings::DNAStringSet(vapply(seq_len(4), function(i) {
    paste(sample(c("A","C","G","T"), 200, replace = TRUE), collapse = "")
  }, character(1)))
  a <- shuffle_sequences(seqs, method = "euler", k = 2, nthreads = 2,
                          rng.seed = 42)
  b <- shuffle_sequences(seqs, method = "euler", k = 2, nthreads = 2,
                          rng.seed = 42)
  expect_equal(as.character(a), as.character(b))
})
