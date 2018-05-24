library(universalmotif)
context("Test filter_motifs")

test_motif1 <- create_motif("AKVAATGCV", nsites = 56)
motif_slots(test_motif1, "strand") <- "-"
motif_slots(test_motif1, "extrachar") <- c("hello" = "goodbye")

test_motif2 <- create_motif("UAUAUAWW", nsites = 40, alphabet = "RNA",
                            out_type = "PCM")
motif_slots(test_motif2, "extranum") <- c("thenum" = 5)

test_list <- list(test_motif1, test_motif2)

test_that("filtering works", {
  expect_equal(length(filter_motifs(test_list, alphabet = "DNA")), 1)
  expect_equal(length(filter_motifs(test_list, nsites = 50)), 1)
  expect_equal(length(filter_motifs(test_list, min_width = 9)), 1)
  expect_equal(length(filter_motifs(test_list, type = "PPM")), 1)
  expect_equal(length(filter_motifs(test_list, icscore = 13.06)), 1)
  expect_equal(length(filter_motifs(test_list, strand = "-")), 2)
  expect_equal(length(filter_motifs(test_list, strand = "+")), 1)
  expect_equal(length(filter_motifs(test_list,
                                    extrachar = c("hello" = "goodbye"))), 1)
  expect_equal(length(filter_motifs(test_list,
                                    extrachar = c("goodbye" = "goodbye"))), 0)
})                            
