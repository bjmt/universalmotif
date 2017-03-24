library(universalmotif)
context("Testing trim_motifs")

test_motif <- create_motif("NNNAAANTTTTNNN")

test_that("trim_motifs leaves inner Ns alone", {
  expect_identical(trim_motifs(test_motif)@consensus, "AAANTTTT")
})
