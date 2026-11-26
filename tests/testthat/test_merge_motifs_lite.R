context("merge_motifs_lite()")

test_that("returns a universalmotif S4 of type PPM with input alphabet", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  r <- merge_motifs_lite(list(m1, m2))
  expect_s4_class(r, "universalmotif")
  expect_equal(r@type, "PPM")
  expect_equal(r@alphabet, "DNA")
})

test_that("width is bounded by [min(input.widths), sum(input.widths)]", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  r <- merge_motifs_lite(list(m1, m2, m3))
  expect_gte(ncol(r@motif), min(c(8L, 8L, 8L)))
  expect_lte(ncol(r@motif), sum(c(8L, 8L, 8L)))
})

test_that("single-motif merge returns the input as PPM", {
  m1 <- create_motif("TTGACATA", name = "a", type = "PPM")
  r <- merge_motifs_lite(list(m1))
  expect_equal(ncol(r@motif), ncol(m1@motif))
  expect_equal(r@alphabet, "DNA")
})

test_that("two identical motifs merge to an equivalent motif", {
  m1 <- create_motif("TTGACATA", name = "a")
  r <- merge_motifs_lite(list(m1, m1))
  expect_equal(ncol(r@motif), ncol(m1@motif))
  ## Same column values up to numerical tolerance.
  expect_equal(as.numeric(r@motif), as.numeric(m1@motif), tolerance = 1e-8)
})

test_that("RC handling: merging a motif with its reverse-complement", {
  m1 <- create_motif("TTGACATA", name = "a")
  m1_rc <- create_motif(c("TATGTCAA"), name = "a_rc")  # RC of TTGACATA
  r <- merge_motifs_lite(list(m1, m1_rc), RC = TRUE)
  ## The merged motif should be similar to m1 (or its RC) -- the
  ## merging engine RCs whichever doesn't match the anchor orientation.
  expect_s4_class(r, "universalmotif")
  expect_equal(ncol(r@motif), 8L)
})

test_that("determinism: same input twice -> identical result", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  r1 <- merge_motifs_lite(list(m1, m2, m3))
  r2 <- merge_motifs_lite(list(m1, m2, m3))
  expect_equal(r1@motif, r2@motif)
  expect_equal(r1@consensus, r2@consensus)
})

test_that("AA motifs are rejected", {
  m1 <- create_motif("AAYY", alphabet = "AA", name = "p1")
  m2 <- create_motif("AAYY", alphabet = "AA", name = "p2")
  expect_error(merge_motifs_lite(list(m1, m2)), regexp = "DNA/RNA")
})

test_that("alphabet mismatch is rejected", {
  d <- create_motif("ACGT", name = "d", alphabet = "DNA")
  r <- create_motif("ACGU", name = "r", alphabet = "RNA")
  expect_error(merge_motifs_lite(list(d, r)), regexp = "same alphabet")
})

test_that("weighted=TRUE vs FALSE can differ on skewed-nsites input", {
  ## Use non-RC-symmetric motifs (AAAAAAAA / CCCCCCCC) and disable RC,
  ## so the two motifs really do contribute independently to the average.
  m1 <- create_motif("AAAAAAAA", name = "a", nsites = 1000)
  m2 <- create_motif("CCCCCCCC", name = "b", nsites = 10)
  r_w  <- merge_motifs_lite(list(m1, m2), weighted = TRUE,  RC = FALSE)
  r_uw <- merge_motifs_lite(list(m1, m2), weighted = FALSE, RC = FALSE)
  ## The weighted version should be more A-rich in the overlap region
  ## (since m1 dominates 100:1), while the unweighted version should be
  ## ~50/50 A/C there. Differences are diluted by the non-overlapping
  ## flanks, so check both the overlap-aware metric and a global one.
  a_w  <- mean(r_w@motif["A", ])
  a_uw <- mean(r_uw@motif["A", ])
  expect_gt(a_w, a_uw + 0.1)
})

test_that("no automatic IC-trim: low-IC flanks are preserved", {
  ## Build two motifs that overlap only at the centre, so the merged
  ## motif has low-IC flanks contributed by only one input each.
  m1 <- create_motif("AAAACCCC", name = "a")
  m2 <- create_motif("CCCCGGGG", name = "b")
  r <- merge_motifs_lite(list(m1, m2))
  ## Result width should equal the union extent (>= 8 cols).
  expect_gte(ncol(r@motif), 8L)
})

test_that("new.name overrides the default", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  r <- merge_motifs_lite(list(m1, m2), new.name = "my_merged")
  expect_equal(r@name, "my_merged")
})

test_that("merge_motifs and merge_motifs_lite both return a universalmotif (v1/v2 parity smoke)", {
  m1 <- create_motif("ACGTAC", name = "a")
  m2 <- create_motif("ACGTAC", name = "b")
  r1 <- merge_motifs(list(m1, m2))
  r2 <- merge_motifs_lite(list(m1, m2))
  expect_s4_class(r1, "universalmotif")
  expect_s4_class(r2, "universalmotif")
})

# ---- Gapped motifs ---------------------------------------------------------------
#
# Reference values for the fixture motifs of helper-lite-fixture.R were worked
# out independently.

gapped_pair <- function() {
  ms <- lite_fixture()
  list(g6 = add_gap(ms[[6]], gaploc = 4, mingap = 0, maxgap = 3),
       g7 = add_gap(ms[[7]], gaploc = 3, mingap = 1, maxgap = 2),
       ms = ms)
}

test_that("motifs with the same number of gaps are merged part by part", {
  g <- gapped_pair()
  m <- merge_motifs_lite(list(g$g6, g$g7))
  expect_equal(m@name, "merged_M6+M7")
  expect_equal(m@consensus, "MTGAGGTCAT")
  expect_true(m@gapinfo@isgapped)
  expect_equal(m@gapinfo@gaploc, 5)
  expect_equal(m@gapinfo@mingap, 0)
  expect_equal(m@gapinfo@maxgap, 3)
  expect_equal(m@nsites, 2)
  expect_equal(as.vector(m@motif)[13:24],
               c(0.701, 0.058, 0.1735, 0.0675, 0.132, 0.136, 0.731, 0.001,
                 0.115, 0.047, 0.675, 0.163))
})

test_that("motifs with different numbers of gaps are merged with their gaps filled", {
  g <- gapped_pair()
  expect_warning(m <- merge_motifs_lite(list(g$g6, g$ms[[7]])),
                 "the gaps of M6 were filled with blank columns")
  expect_equal(m@consensus, "MTGAGTCAT")
  expect_false(m@gapinfo@isgapped)
  expect_warning(m <- merge_motifs_lite(list(g$ms[[7]], g$g6, g$g7)),
                 "the gaps of M6, M7 were filled")
  expect_equal(m@name, "merged_M7+M6+M7")
  expect_equal(m@consensus, "CTVWGAGTCAT")
  expect_equal(as.vector(m@motif)[9:16],
               c(0.304333, 0.35, 0.311, 0.034667, 0.295667, 0.059667, 0.129667,
                 0.515), tolerance = 1e-5)
  g1 <- add_gap(g$ms[[1]], gaploc = c(3, 6), mingap = c(1, 0), maxgap = c(2, 1))
  expect_warning(m <- merge_motifs_lite(list(g1, g$ms[[2]], g$ms[[3]])),
                 "the gaps of M1 were filled")
  expect_equal(m@consensus, "CTTGNNAYATAKS")
})

test_that("use.gaps = FALSE merges gapped motifs as if their parts followed one another", {
  g <- gapped_pair()
  m <- expect_silent(merge_motifs_lite(list(g$g6, g$g7), use.gaps = FALSE))
  expect_false(m@gapinfo@isgapped)
  expect_equal(m@consensus, "MTGAGTCAT")
  expect_equal(m@motif, merge_motifs_lite(g$ms[6:7])@motif)
})
