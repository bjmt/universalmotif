context("trim_motifs()")

test_that("motif trimming works", {

  m1 <- create_motif("NNCCCNN", nsites = 100)
  m2 <- trim_motifs(m1)

  expect_equal(ncol(m2), 3)

})

test_that("trim_motifs() preserves successfully-trimmed motifs when one is fully trimmed (regression: return(NULL) on partial trim)", {

  m1 <- create_motif("ACGTACGT", pseudocount = 1, nsites = 100)
  m2 <- create_motif(matrix(rep(0.25, 16), nrow = 4,
                             dimnames = list(c("A", "C", "G", "T"), NULL)),
                     alphabet = "DNA", type = "PPM")
  out <- suppressMessages(trim_motifs(c(m1, m2), min.ic = 0.5))
  out_list <- if (is.list(out)) out else list(out)
  expect_true(length(out_list) >= 1)
  expect_equal(out_list[[1]]@name, m1@name)

})
