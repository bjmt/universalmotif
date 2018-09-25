context("Test motif scanning")

test_that("Results are accurate", {

  motif <- create_motif("AAAA", pseudocount = 1, nsites = 100)
  seq <- Biostrings::DNAStringSet("GGGAAAAGGGTTTTGGG")  # width 17
  res <- scan_sequences(motif, seq, RC = TRUE, verbose = 0, progress = FALSE)

  expect_equal(res$start[1], 4)
  expect_equal(res$motif[1], "motif")
  expect_equal(res$sequence[1], "1")
  expect_equal(res$stop[1], 7)
  expect_equal(res$match[1], "AAAA")
  expect_equal(res$strand[1], "+")
  expect_equal(res$score[1], 7.956, tolerance = 0.001)
  expect_equal(res$`max.score`[1], 7.957, tolerance = 0.001)
  expect_equal(res$`score.pct`[1], 99.99, tolerance = 0.01) 
  
  expect_equal(res$score[1], res$score[2])
  expect_equal(res$`max.score`[1], res$`max.score`[2])
  expect_equal(res$strand[2], "-")
  expect_equal(res$start[2], 14)
  expect_equal(res$stop[2], 11)

})
