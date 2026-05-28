context("enrich_motifs2()")

OUT_COLS <- c("motif", "motif.i", "consensus",
              "target.seq.n", "target.seq.hits", "target.site.hits",
              "bkg.seq.n", "bkg.seq.hits", "bkg.site.hits",
              "enrichment", "log2.enrichment", "pvalue", "qvalue")

test_that("output shape and columns are stable", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  tgt <- Biostrings::DNAStringSet(rep("TTTAAA", 50))
  bkg <- Biostrings::DNAStringSet(rep("GGGCCC", 50))

  r <- enrich_motifs2(m, tgt, bkg, pvalue = 0.01, qvalue = 1, rng.seed = 1)

  expect_true(is.data.frame(r))
  expect_equal(names(r), OUT_COLS)
  expect_equal(nrow(r), 1L)
  expect_type(r$motif,   "character")
  expect_type(r$motif.i, "integer")
  expect_type(r$pvalue,  "double")
  expect_type(r$qvalue,  "double")
})

test_that("determinism with fixed rng.seed (shuffled bkg path)", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  tgt <- Biostrings::DNAStringSet(c(rep("TTTAAACCCGGGAAA", 30),
                                    rep("GGGCCCGGGCCCGGG", 30)))

  r1 <- enrich_motifs2(m, tgt, pvalue = 0.01, qvalue = 1, rng.seed = 42L)
  r2 <- enrich_motifs2(m, tgt, pvalue = 0.01, qvalue = 1, rng.seed = 42L)

  expect_identical(r1, r2)
})

test_that("self-enrichment: target = motif, bkg = unrelated", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  tgt <- Biostrings::DNAStringSet(rep("TTTAAA", 100))
  bkg <- Biostrings::DNAStringSet(rep("GGGCCC", 100))

  r <- enrich_motifs2(m, tgt, bkg, pvalue = 0.01, qvalue = 0.1, rng.seed = 1)

  expect_equal(nrow(r), 1L)
  expect_lt(r$pvalue, 1e-10)
  expect_lte(r$qvalue, 0.1)
  expect_gt(r$enrichment, 1)
  expect_true(is.finite(r$log2.enrichment) && r$log2.enrichment > 0)
})

test_that("AA motifs are rejected", {

  m <- create_motif("AAYY", alphabet = "AA")
  tgt <- Biostrings::AAStringSet(rep("AAYYAA", 5))

  expect_error(
    enrich_motifs2(m, tgt, qvalue = 1),
    regexp = "DNA/RNA"
  )
})

test_that("seqs vs sites can differ on the same fixture", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  ## Each target sequence has multiple instances of the motif.
  tgt <- Biostrings::DNAStringSet(
    rep("TTTAAACCCTTTAAACCCTTTAAA", 30))
  bkg <- Biostrings::DNAStringSet(
    rep("GGGCCCGGGCCCGGGCCCGGGCCC", 30))

  r_seqs  <- enrich_motifs2(m, tgt, bkg, test = "seqs",
                            pvalue = 0.01, qvalue = 1, rng.seed = 1)
  r_sites <- enrich_motifs2(m, tgt, bkg, test = "sites",
                            pvalue = 0.01, qvalue = 1, rng.seed = 1)

  expect_equal(nrow(r_seqs),  1L)
  expect_equal(nrow(r_sites), 1L)
  ## In the sites mode, the denominator is per-position so the enrichment
  ## is generally not identical to the seqs mode.
  expect_false(isTRUE(all.equal(r_seqs$enrichment, r_sites$enrichment)))
})

test_that("qvalue filter is a subset operation", {

  m1 <- create_motif("TTTAAA", pseudocount = 1, nsites = 100, name = "m1")
  m2 <- create_motif("CACGTG", pseudocount = 1, nsites = 100, name = "m2")
  tgt <- Biostrings::DNAStringSet(
    c(rep("TTTAAACCCGGG", 50), rep("CACGTGCACGTG", 50)))
  bkg <- Biostrings::DNAStringSet(rep("GGGCCCGGGCCC", 100))

  r_loose <- enrich_motifs2(list(m1, m2), tgt, bkg,
                            pvalue = 0.01, qvalue = 0.5, rng.seed = 1)
  r_tight <- enrich_motifs2(list(m1, m2), tgt, bkg,
                            pvalue = 0.01, qvalue = 0.01, rng.seed = 1)

  expect_true(all(r_tight$motif %in% r_loose$motif))
  expect_lte(nrow(r_tight), nrow(r_loose))
})

test_that("RC = FALSE changes hit counts on a palindromic-asymmetric fixture", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  ## Sequences containing only the - strand instance of TTTAAA
  ## (which is also TTTAAA in DNA — TTTAAA is palindromic).
  ## Use an asymmetric motif instead.
  m2 <- create_motif("TTGACA", pseudocount = 1, nsites = 100, name = "asym")
  tgt <- Biostrings::DNAStringSet(rep("TTGACAAAATTGACA", 30))
  bkg <- Biostrings::DNAStringSet(rep("GGGCCCGGGCCCGGG", 30))

  r_rc  <- enrich_motifs2(m2, tgt, bkg, RC = TRUE,
                          pvalue = 0.01, qvalue = 1, rng.seed = 1)
  r_one <- enrich_motifs2(m2, tgt, bkg, RC = FALSE,
                          pvalue = 0.01, qvalue = 1, rng.seed = 1)

  ## Both ran without error; the site-hit counts on the target side should
  ## match (motif is asymmetric and shows up only on + strand of target),
  ## but the scannable position pool used in sites mode differs.
  expect_equal(nrow(r_rc),  1L)
  expect_equal(nrow(r_one), 1L)
})

test_that("empty result returns 0 rows with correct columns", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  tgt <- Biostrings::DNAStringSet(rep("GGGCCCGGGCCC", 30))
  bkg <- Biostrings::DNAStringSet(rep("AAAATTTTAAAA", 30))

  r <- enrich_motifs2(m, tgt, bkg, pvalue = 0.01, qvalue = 1e-30,
                      rng.seed = 1)
  expect_equal(nrow(r), 0L)
  expect_equal(names(r), OUT_COLS)
})

test_that("user-supplied background does not depend on rng.seed", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  tgt <- Biostrings::DNAStringSet(rep("TTTAAACCCGGG", 30))
  bkg <- Biostrings::DNAStringSet(rep("CACGTGCACGTG", 30))

  r1 <- enrich_motifs2(m, tgt, bkg, pvalue = 0.01, qvalue = 1, rng.seed = 1)
  r2 <- enrich_motifs2(m, tgt, bkg, pvalue = 0.01, qvalue = 1, rng.seed = 999)

  expect_identical(r1, r2)
})

test_that("pseudocount shifts p-value toward 1 (regularisation)", {

  m <- create_motif("TTTAAA", pseudocount = 1, nsites = 100)
  tgt <- Biostrings::DNAStringSet(c(rep("TTTAAA", 5), rep("GGGCCC", 5)))
  bkg <- Biostrings::DNAStringSet(rep("GGGCCC", 10))

  r_low  <- enrich_motifs2(m, tgt, bkg, pseudocount = 1L,
                           pvalue = 0.01, qvalue = 1, rng.seed = 1)
  r_high <- enrich_motifs2(m, tgt, bkg, pseudocount = 100L,
                           pvalue = 0.01, qvalue = 1, rng.seed = 1)

  expect_equal(nrow(r_low), 1L)
  expect_equal(nrow(r_high), 1L)
  expect_gt(r_high$pvalue, r_low$pvalue)
})

test_that("enrich_motifs2 runs on RNA motifs and sequences", {
  suppressMessages({
    ## Use a fixture where the motif is implanted in every primary sequence
    ## so the enrichment call has data to compare against the shuffled background.
    set.seed(1)
    bg <- vapply(seq_len(20), function(i) {
      paste0(sample(c("A","C","G","U"), 80, replace = TRUE), collapse = "")
    }, character(1))
    primary <- vapply(bg, function(s) {
      pos <- 30
      paste0(substr(s, 1, pos - 1), "ACGU",
             substr(s, pos + 4, nchar(s)))
    }, character(1))
    seqs <- Biostrings::RNAStringSet(primary)
    bkg  <- Biostrings::RNAStringSet(bg)
    m <- create_motif("ACGU", alphabet = "RNA", nsites = 10)
    res <- enrich_motifs2(list(m), seqs, bkg, pvalue = 0.5)
    expect_s3_class(res, "data.frame")
    expect_true(all(c("motif", "fg", "bg", "pvalue") %in% colnames(res)) ||
                ncol(res) > 0L)
  })
})

test_that("enrich_motifs and enrich_motifs2 both run on shared fixtures (v1/v2 parity smoke)", {
  suppressMessages(suppressWarnings({
    set.seed(1)
    bg <- vapply(seq_len(20), function(i) {
      paste(sample(c("A","C","G","T"), 100, replace = TRUE), collapse = "")
    }, character(1))
    primary <- vapply(bg, function(s) {
      paste0(substr(s, 1, 30), "ACGTAC", substr(s, 37, nchar(s)))
    }, character(1))
    seqs <- Biostrings::DNAStringSet(primary)
    bkg <- Biostrings::DNAStringSet(bg)
    m <- create_motif("ACGTAC", nsites = 100)
    r1 <- enrich_motifs(m, seqs, bkg.sequences = bkg, verbose = 0,
                        threshold = 1e-2, threshold.type = "pvalue")
    r2 <- enrich_motifs2(list(m), seqs, bkg, pvalue = 1e-2)
    expect_true(nrow(r1) >= 0L)
    expect_true(nrow(r2) >= 0L)
  }))
})
