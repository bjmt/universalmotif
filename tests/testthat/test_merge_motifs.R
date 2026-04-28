context("merge_motifs()")

test_that("motif merging works", {

  m1 <- create_motif("AAAAAAA", nsites = 10)
  m2 <- create_motif("TTTTTTT", nsites = 10)

  m3 <- merge_motifs(list(m1, m2), tryRC = FALSE, min.overlap = 10)

  expect_equal(m3@consensus, "WWWWWWW")

})

test_that("merge_motifs() preserves bkgsites when exactly one motif has it (regression: > 1 should be >= 1)", {

  m1 <- create_motif("ACGT", pseudocount = 1, nsites = 100)
  m2 <- create_motif("ACGT", pseudocount = 1, nsites = 100)
  m1["bkgsites"] <- 500
  merged <- merge_motifs(c(m1, m2))
  expect_equal(merged@bkgsites, 500)

})

test_that("merge_motifs() leaves pval/qval/eval empty when no input has them (regression: ifelse(., numeric()))", {

  m1 <- create_motif("ACGT", pseudocount = 1, nsites = 100)
  m2 <- create_motif("ACGT", pseudocount = 1, nsites = 100)
  merged <- merge_motifs(c(m1, m2))
  expect_equal(length(merged@pval), 0L)
  expect_equal(length(merged@qval), 0L)
  expect_equal(length(merged@eval), 0L)
  expect_false(any(is.na(c(merged@pval, merged@qval, merged@eval))))

})
