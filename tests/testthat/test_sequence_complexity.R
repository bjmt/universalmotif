context("sequence_complexity()")

library(Biostrings)

methods <- c("WoottonFederhen", "WoottonFederhenFast",
             "Trifonov", "TrifonovFast", "DUST")

test_that("each method runs and returns a DataFrame with one row per window", {
  seqs <- DNAStringSet(c("ACGTACGTACGTACGTACGT",
                         "AAAAAAAAAAAAAAAAAAAA"))
  for (m in methods) {
    res <- sequence_complexity(seqs, window.size = 10, window.overlap = 5,
                                method = m)
    expect_s4_class(res, "DataFrame")
    expect_true(all(c("sequence", "start", "stop", "complexity") %in%
                    colnames(res)))
    expect_equal(nrow(res), 2L * 3L)
    expect_true(all(is.finite(res$complexity)))
  }
})

test_that("homopolymer scores at the low-complexity end for each method", {
  seqs <- DNAStringSet(c("AAAAAAAAAAAAAAAAAAAA",
                         "ACGTACGTACGTACGTACGT"))
  for (m in methods) {
    res <- sequence_complexity(seqs, window.size = 20, method = m)
    homo <- res$complexity[res$sequence == "1"]
    mixed <- res$complexity[res$sequence == "2"]
    if (m == "DUST") {
      ## DUST is reversed: higher = less complex.
      expect_true(homo > mixed)
    } else {
      expect_true(homo < mixed)
    }
  }
})

test_that("AA sequences work with WoottonFederhen", {
  aaseq <- AAStringSet(c("ACDEFGHIKLMNPQRSTVWY",
                         "AAAAAAAAAAAAAAAAAAAA"))
  res <- sequence_complexity(aaseq, window.size = 20,
                              method = "WoottonFederhen")
  expect_s4_class(res, "DataFrame")
  expect_equal(nrow(res), 2L)
  expect_true(all(is.finite(res$complexity)))
  expect_true(res$complexity[res$sequence == "1"] >
              res$complexity[res$sequence == "2"])
})

test_that("nthreads parity (1 vs 2) on a single method", {
  seqs <- DNAStringSet(vapply(seq_len(8), function(i) {
    set.seed(i)
    paste0(sample(c("A","C","G","T"), 200, replace = TRUE), collapse = "")
  }, character(1)))
  a <- sequence_complexity(seqs, window.size = 20, window.overlap = 10,
                            method = "DUST", nthreads = 1)
  b <- sequence_complexity(seqs, window.size = 20, window.overlap = 10,
                            method = "DUST", nthreads = 2)
  expect_equal(a$complexity, b$complexity)
})
