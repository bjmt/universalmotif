context("enrich_motifs()")

test_that("motif enrichment works", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  s1 <- Biostrings::DNAStringSet(rep("TTTAAA", 100))
  s2 <- Biostrings::DNAStringSet(c(rep("CCCGGG", 100), "TTTAAA"))
  r <- suppressWarnings(
         enrich_motifs(m, s1, s2, verbose = 0, progress = FALSE, threshold = 0.8,
                       threshold.type = "logodds", max.p = 0.1, max.q = 0.1,
                       max.e = 0.1)
        )

  r2 <- suppressWarnings(
          enrich_motifs(m, s1, s2, verbose = 0, progress = FALSE,
                        search.mode = "both", return.scan.results = TRUE)
        )

  expect_true(is.data.frame(r))
  expect_equal(nrow(r), 1)

  expect_true(is.list(r2))
  expect_equal(names(r2), c("enrichment.report", "input.scan", "bkg.scan"))

})
