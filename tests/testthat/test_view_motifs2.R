context("view_motifs2()")

test_that("view_motifs2 returns a ggplot for a single motif", {
  m <- create_motif("CACGTG", name = "ebox")
  p <- view_motifs2(m)
  expect_s3_class(p, "ggplot")
})

test_that("view_motifs2 returns a ggplot for a multi-motif list", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  p <- view_motifs2(list(m1, m2, m3))
  expect_s3_class(p, "ggplot")
})

test_that("view_motifs2 detects reverse-complement and tags name with RC.text", {
  ## TTGACATA is non-palindromic, so its RC is genuinely different
  ## and the aligner should flag it.
  m1 <- create_motif("TTGACATA", name = "creb")
  m2 <- motif_rc(m1); m2["name"] <- "creb_rc"
  raw <- view_motifs2(list(m1, m2), return.raw = TRUE)
  expect_equal(length(raw), 2L)
  expect_true(any(grepl("\\[RC\\]", names(raw))))
})

test_that("view_motifs2 with return.raw produces matrices of equal width", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  raw <- view_motifs2(list(m1, m2, m3), return.raw = TRUE)
  expect_equal(length(raw), 3L)
  expect_equal(length(unique(vapply(raw, ncol, integer(1)))), 1L)
})

test_that("view_motifs2 sort.by accepts none / ic / similarity", {
  motifs <- lapply(c("CACGTG", "TGACGT", "AAATTT", "GCGCGC"),
                   function(s) create_motif(s, name = s))
  expect_s3_class(view_motifs2(motifs, sort.by = "none"), "ggplot")
  expect_s3_class(view_motifs2(motifs, sort.by = "ic"), "ggplot")
  expect_s3_class(view_motifs2(motifs, sort.by = "similarity"), "ggplot")
})

test_that("view_motifs2 names.pos = 'right' returns a ggplot", {
  motifs <- list(create_motif("CACGTG", name = "ebox"),
                 create_motif("TTGACATA", name = "creb"))
  expect_s3_class(view_motifs2(motifs, names.pos = "right"), "ggplot")
})

test_that("view_motifs2 rejects non-DNA / non-RNA alphabets", {
  aa <- create_motif("ACDEFGHIKL", alphabet = "AA", name = "aa1")
  expect_error(view_motifs2(aa), "DNA/RNA")
})

test_that("view_motifs2 rejects mixed alphabets", {
  dna <- create_motif("ACGT", name = "dna1")
  rna <- create_motif("ACGU", alphabet = "RNA", name = "rna1")
  expect_error(view_motifs2(list(dna, rna)), "alphabet")
})

test_that("view_motifs2 dedup.names default rewrites duplicates", {
  m1 <- create_motif("CACGTG", name = "same")
  m2 <- create_motif("TGACGT", name = "same")
  raw <- view_motifs2(list(m1, m2), return.raw = TRUE)
  expect_equal(length(unique(names(raw))), 2L)
})
