context("enrich_motifs()")

test_that("motif enrichment works", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  s1 <- Biostrings::DNAStringSet(rep("TTTAAA", 100))
  s2 <- Biostrings::DNAStringSet(c(rep("CCCGGG", 100), "TTTAAA"))
  r <- suppressWarnings(
         enrich_motifs(m, s1, s2, verbose = 0, threshold = 0.8,
                       threshold.type = "logodds", max.p = 0.1, max.q = 0.1,
                       max.e = 0.1)
        )

  r2 <- suppressWarnings(
          enrich_motifs(m, s1, s2, verbose = 0,
                        return.scan.results = TRUE)
        )

  expect_true(is(r, "DataFrame"))
  expect_equal(nrow(r), 1)

  expect_true(is(r2, "DataFrame"))
  expect_equal(names(S4Vectors::metadata(r2)),
               c("scan.target", "scan.bkg", "args"))

})

test_that("enrich_motifs() gives a clear error for sequences shorter than motif (regression: was crashing in fisher.test)", {

  m <- create_motif("TTTAAA", pseudocount = 1)
  too_short <- Biostrings::DNAStringSet("A")
  expect_error(enrich_motifs(m, too_short, verbose = 0), regexp = "widest motif")

})

test_that("enrich_motifs runs on RNA motifs and sequences", {
  suppressMessages(suppressWarnings({
    m <- create_motif("ACGU", alphabet = "RNA", nsites = 10)
    set.seed(1)
    seqs <- Biostrings::RNAStringSet(vapply(seq_len(20), function(i) {
      paste0(sample(c("A","C","G","U"), 80, replace = TRUE), collapse = "")
    }, character(1)))
    bkg <- shuffle_sequences(seqs, k = 1)
    res <- enrich_motifs(m, seqs, bkg.sequences = bkg, verbose = 0,
                         threshold = 0.5, threshold.type = "logodds.abs")
    expect_true(is.data.frame(res) || isS4(res))
  }))
})
