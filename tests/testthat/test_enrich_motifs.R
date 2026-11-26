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

test_that("enrich_motifs() corrects the threshold of gapped motifs as scan_sequences() does", {
  op <- options(universalmotif.suggest.enrich_motifs_lite = FALSE)
  on.exit(options(op), add = TRUE)
  j <- read_jaspar(system.file("extdata", "jaspar.txt", package = "universalmotif"))
  m <- j[[which.max(vapply(j, function(x) ncol(x@motif), integer(1)))]]
  m <- add_gap(m, gaploc = 5, mingap = 0, maxgap = 5)
  s <- create_sequences(seqnum = 50, seqlen = 1000, rng.seed = 1)
  b <- create_sequences(seqnum = 50, seqlen = 1000, rng.seed = 2)
  for (cg in c(FALSE, TRUE)) {
    r <- suppressMessages(enrich_motifs(m, s, b, max.p = 1, max.q = 1, max.e = Inf,
      threshold = 1e-3, correct.gaps = cg, RC = FALSE, no.overlaps = FALSE))
    x <- scan_sequences(m, s, threshold = 1e-3, threshold.type = "pvalue",
      RC = FALSE, correct.gaps = cg, verbose = 0)
    expect_equal(r$target.hits, nrow(x))
  }
  r0 <- suppressMessages(enrich_motifs(m, s, b, max.p = 1, max.q = 1, max.e = Inf,
    threshold = 1e-3, correct.gaps = FALSE, RC = FALSE, no.overlaps = FALSE))
  r1 <- suppressMessages(enrich_motifs(m, s, b, max.p = 1, max.q = 1, max.e = Inf,
    threshold = 1e-3, RC = FALSE, no.overlaps = FALSE))
  expect_lt(r1$target.hits, r0$target.hits)
})
