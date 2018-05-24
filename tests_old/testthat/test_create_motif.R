library(universalmotif)
context("Test create_motif")

test_consensus <- create_motif(consensus = "YYDWARTT", name = "test")

test_that("create_motif from a consensus string works", {
  expect_equal(motif_slots(test_consensus, "icscore"),
               c("icscore" = 10.06), tolerance = 0.009)
  expect_identical(motif_slots(test_consensus, "consensus"),
                   c("consensus" = "YYDWARTT"))
})
