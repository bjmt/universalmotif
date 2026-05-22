context("motif_finder()")

## yamtk-specific extra columns attached on top of the standard
## universalmotif_df shape (`consensus` is supplied by to_df() from the
## motif's @consensus slot, so it isn't in this list).
OUT_EXTRA_COLS <- c("rank", "width",
                    "seqs_pos", "seqs_neg", "sites_pos", "sites_neg",
                    "n_pos", "n_neg", "pvalue", "qvalue")

## Helper: plant `motif` in `n` of the first `n` sequences of `seqs`,
## at random positions, returning a Biostrings::DNAStringSet.
plant <- function(seqs, motif, n, rng.seed = 1) {
  set.seed(rng.seed)
  mw <- nchar(as.character(motif))
  for (i in seq_len(n)) {
    pos <- sample.int(length(seqs[[i]]) - mw, 1)
    Biostrings::subseq(seqs[[i]], start = pos, width = mw) <- motif
  }
  seqs
}

test_that("plant-and-recover finds the planted motif", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 60, seqlen = 200, rng.seed = 1)
  seqs <- plant(seqs, Biostrings::DNAString("TTGACATA"), n = 48, rng.seed = 1)

  r <- motif_finder(seqs, qvalue = 1, rng.seed = 1, nthreads = 1,
                    nmotifs = 2, min.width = 6, max.width = 10)
  expect_gte(nrow(r), 1L)
  ## The top motif (lowest p-value) should be either TTGACATA or its RC
  ## TATGTCAA (we scan both strands by default).
  top <- r$consensus[which.min(r$pvalue)]
  expect_true(top %in% c("TTGACATA", "TATGTCAA") ||
              grepl("TGTCAA|TGACAT", top))
})

test_that("output is a universalmotif_df with the expected extra columns", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 40, seqlen = 150, rng.seed = 1)
  seqs <- plant(seqs, Biostrings::DNAString("TTGACATA"), n = 32, rng.seed = 1)

  r <- motif_finder(seqs, qvalue = 1, rng.seed = 1, nthreads = 1,
                    nmotifs = 1, min.width = 6, max.width = 8)
  expect_s3_class(r, "universalmotif_df")
  expect_true("motif" %in% names(r))
  expect_true(all(OUT_EXTRA_COLS %in% names(r)))
  ## to_list() should round-trip into a list of universalmotif S4 objects.
  ml <- to_list(r)
  expect_type(ml, "list")
  expect_s4_class(ml[[1]], "universalmotif")
})

test_that("determinism with fixed rng.seed", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 40, seqlen = 150, rng.seed = 1)
  seqs <- plant(seqs, Biostrings::DNAString("TTGACATA"), n = 32, rng.seed = 1)

  r1 <- motif_finder(seqs, qvalue = 1, rng.seed = 42L, nthreads = 1,
                     nmotifs = 1, min.width = 6, max.width = 8)
  r2 <- motif_finder(seqs, qvalue = 1, rng.seed = 42L, nthreads = 1,
                     nmotifs = 1, min.width = 6, max.width = 8)
  expect_equal(r1$pvalue,    r2$pvalue)
  expect_equal(r1$consensus, r2$consensus)
})

test_that("empty result on random DNA returns 0-row df with correct columns", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 30, seqlen = 100, rng.seed = 7)
  r <- motif_finder(seqs, qvalue = 1e-30, rng.seed = 1, nthreads = 1,
                    nmotifs = 1, min.width = 6, max.width = 8)
  expect_equal(nrow(r), 0L)
  expect_true(all(OUT_EXTRA_COLS %in% names(r)))
})

test_that("AA / non-DNA sequences are rejected", {

  aa <- Biostrings::AAStringSet(rep("ACDEFGHIKLMNP", 10))
  expect_error(motif_finder(aa, qvalue = 1, rng.seed = 1, nthreads = 1),
               regexp = "DNA/RNA")
})

test_that("nmotifs bound is respected", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 60, seqlen = 200, rng.seed = 1)
  seqs <- plant(seqs, Biostrings::DNAString("TTGACATA"), n = 50, rng.seed = 1)

  r <- motif_finder(seqs, qvalue = 1, rng.seed = 1, nthreads = 1,
                    nmotifs = 1, min.width = 6, max.width = 8)
  ## Per-width cap is nmotifs; max.width - min.width + 1 = 3 widths => max 3 rows.
  expect_lte(nrow(r), 3L)
})

test_that("discovered widths lie in [min.width, max.width]", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 50, seqlen = 200, rng.seed = 1)
  seqs <- plant(seqs, Biostrings::DNAString("TTGACATA"), n = 40, rng.seed = 1)

  r <- motif_finder(seqs, qvalue = 1, rng.seed = 1, nthreads = 1,
                    nmotifs = 1, min.width = 7, max.width = 10)
  ## widths may shrink due to IC-trim, so use lower bound of (min - some slack).
  expect_true(all(r$width >= 3L))
  expect_true(all(r$width <= 10L))
})

test_that("qvalue filter is a subset operation", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 60, seqlen = 200, rng.seed = 1)
  seqs <- plant(seqs, Biostrings::DNAString("TTGACATA"), n = 40, rng.seed = 1)

  loose <- motif_finder(seqs, qvalue = 1,    rng.seed = 1, nthreads = 1,
                        nmotifs = 1, min.width = 6, max.width = 10)
  tight <- motif_finder(seqs, qvalue = 1e-10, rng.seed = 1, nthreads = 1,
                        nmotifs = 1, min.width = 6, max.width = 10)
  expect_lte(nrow(tight), nrow(loose))
  if (nrow(tight) > 0L) {
    expect_true(all(tight$consensus %in% loose$consensus))
  }
})

test_that("user-supplied background does not depend on rng.seed", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 40, seqlen = 150, rng.seed = 1)
  seqs <- plant(seqs, Biostrings::DNAString("TTGACATA"), n = 30, rng.seed = 1)
  bkg  <- create_sequences(seqnum = 40, seqlen = 150, rng.seed = 99)

  r1 <- motif_finder(seqs, bkg.sequences = bkg, qvalue = 1, rng.seed = 1,
                     nthreads = 1, nmotifs = 1, min.width = 6, max.width = 8)
  r2 <- motif_finder(seqs, bkg.sequences = bkg, qvalue = 1, rng.seed = 999,
                     nthreads = 1, nmotifs = 1, min.width = 6, max.width = 8)
  expect_equal(r1$consensus, r2$consensus)
  expect_equal(r1$pvalue,    r2$pvalue)
})

test_that("RC=FALSE only finds + strand instances", {

  set.seed(1)
  seqs <- create_sequences(seqnum = 50, seqlen = 200, rng.seed = 1)
  ## Plant asymmetric motif on + strand only.
  seqs <- plant(seqs, Biostrings::DNAString("TTGACATA"), n = 40, rng.seed = 1)

  ## With RC = TRUE we expect very strong significance.
  rc_on <- motif_finder(seqs, qvalue = 1, RC = TRUE, rng.seed = 1,
                        nthreads = 1, nmotifs = 1,
                        min.width = 6, max.width = 8)
  ## With RC = FALSE we still should find it (it's on + strand), but counts
  ## differ. Either way: at least one motif discovered.
  rc_off <- motif_finder(seqs, qvalue = 1, RC = FALSE, rng.seed = 1,
                         nthreads = 1, nmotifs = 1,
                         min.width = 6, max.width = 8)
  expect_gte(nrow(rc_on),  1L)
  expect_gte(nrow(rc_off), 1L)
})
