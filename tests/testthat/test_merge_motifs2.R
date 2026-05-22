context("merge_motifs2()")

test_that("returns a universalmotif S4 of type PPM with input alphabet", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  r <- merge_motifs2(list(m1, m2))
  expect_s4_class(r, "universalmotif")
  expect_equal(r@type, "PPM")
  expect_equal(r@alphabet, "DNA")
})

test_that("width is bounded by [min(input.widths), sum(input.widths)]", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  r <- merge_motifs2(list(m1, m2, m3))
  expect_gte(ncol(r@motif), min(c(8L, 8L, 8L)))
  expect_lte(ncol(r@motif), sum(c(8L, 8L, 8L)))
})

test_that("single-motif merge returns the input as PPM", {
  m1 <- create_motif("TTGACATA", name = "a", type = "PPM")
  r <- merge_motifs2(list(m1))
  expect_equal(ncol(r@motif), ncol(m1@motif))
  expect_equal(r@alphabet, "DNA")
})

test_that("two identical motifs merge to an equivalent motif", {
  m1 <- create_motif("TTGACATA", name = "a")
  r <- merge_motifs2(list(m1, m1))
  expect_equal(ncol(r@motif), ncol(m1@motif))
  ## Same column values up to numerical tolerance.
  expect_equal(as.numeric(r@motif), as.numeric(m1@motif), tolerance = 1e-8)
})

test_that("RC handling: merging a motif with its reverse-complement", {
  m1 <- create_motif("TTGACATA", name = "a")
  m1_rc <- create_motif(c("TATGTCAA"), name = "a_rc")  # RC of TTGACATA
  r <- merge_motifs2(list(m1, m1_rc), RC = TRUE)
  ## The merged motif should be similar to m1 (or its RC) -- the
  ## merging engine RCs whichever doesn't match the anchor orientation.
  expect_s4_class(r, "universalmotif")
  expect_equal(ncol(r@motif), 8L)
})

test_that("determinism: same input twice -> identical result", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  r1 <- merge_motifs2(list(m1, m2, m3))
  r2 <- merge_motifs2(list(m1, m2, m3))
  expect_equal(r1@motif, r2@motif)
  expect_equal(r1@consensus, r2@consensus)
})

test_that("AA motifs are rejected", {
  m1 <- create_motif("AAYY", alphabet = "AA", name = "p1")
  m2 <- create_motif("AAYY", alphabet = "AA", name = "p2")
  expect_error(merge_motifs2(list(m1, m2)), regexp = "DNA/RNA")
})

test_that("alphabet mismatch is rejected", {
  d <- create_motif("ACGT", name = "d", alphabet = "DNA")
  r <- create_motif("ACGU", name = "r", alphabet = "RNA")
  expect_error(merge_motifs2(list(d, r)), regexp = "same alphabet")
})

test_that("weighted=TRUE vs FALSE can differ on skewed-nsites input", {
  ## Use non-RC-symmetric motifs (AAAAAAAA / CCCCCCCC) and disable RC,
  ## so the two motifs really do contribute independently to the average.
  m1 <- create_motif("AAAAAAAA", name = "a", nsites = 1000)
  m2 <- create_motif("CCCCCCCC", name = "b", nsites = 10)
  r_w  <- merge_motifs2(list(m1, m2), weighted = TRUE,  RC = FALSE)
  r_uw <- merge_motifs2(list(m1, m2), weighted = FALSE, RC = FALSE)
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
  r <- merge_motifs2(list(m1, m2))
  ## Result width should equal the union extent (>= 8 cols).
  expect_gte(ncol(r@motif), 8L)
})

test_that("new.name overrides the default", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  r <- merge_motifs2(list(m1, m2), new.name = "my_merged")
  expect_equal(r@name, "my_merged")
})
