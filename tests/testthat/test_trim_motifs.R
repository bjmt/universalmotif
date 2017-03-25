library(universalmotif)
context("Testing trim_motifs")

test_motif <- create_motif("NNNAAANTTTTNNN")
test_motif2 <- create_motif("NNNNANNNN")
test_motif3 <- create_motif("NNN")

test_that("trim_motifs works", {
  expect_identical(trim_motifs(test_motif)@consensus, "AAANTTTT")
  expect_identical(trim_motifs(test_motif2)@consensus, "A")
  expect_identical(trim_motifs(test_motif3), NULL)
})
