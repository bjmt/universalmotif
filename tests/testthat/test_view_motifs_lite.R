context("view_motifs_lite()")

test_that("view_motifs_lite returns a ggplot for a single motif", {
  m <- create_motif("CACGTG", name = "ebox")
  p <- view_motifs_lite(m)
  expect_s3_class(p, "ggplot")
})

test_that("view_motifs_lite returns a ggplot for a multi-motif list", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  p <- view_motifs_lite(list(m1, m2, m3))
  expect_s3_class(p, "ggplot")
})

test_that("view_motifs_lite detects reverse-complement and tags name with RC.text", {
  ## TTGACATA is non-palindromic, so its RC is genuinely different
  ## and the aligner should flag it.
  m1 <- create_motif("TTGACATA", name = "creb")
  m2 <- motif_rc(m1); m2["name"] <- "creb_rc"
  raw <- view_motifs_lite(list(m1, m2), return.raw = TRUE)
  expect_equal(length(raw), 2L)
  expect_true(any(grepl("\\[RC\\]", names(raw))))
})

test_that("view_motifs_lite with return.raw produces matrices of equal width", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  raw <- view_motifs_lite(list(m1, m2, m3), return.raw = TRUE)
  expect_equal(length(raw), 3L)
  expect_equal(length(unique(vapply(raw, ncol, integer(1)))), 1L)
})

test_that("view_motifs_lite sort.by accepts none / ic / similarity", {
  motifs <- lapply(c("CACGTG", "TGACGT", "AAATTT", "GCGCGC"),
                   function(s) create_motif(s, name = s))
  expect_s3_class(view_motifs_lite(motifs, sort.by = "none"), "ggplot")
  expect_s3_class(view_motifs_lite(motifs, sort.by = "ic"), "ggplot")
  expect_s3_class(view_motifs_lite(motifs, sort.by = "similarity"), "ggplot")
})

test_that("view_motifs_lite names.pos = 'right' returns a ggplot", {
  motifs <- list(create_motif("CACGTG", name = "ebox"),
                 create_motif("TTGACATA", name = "creb"))
  expect_s3_class(view_motifs_lite(motifs, names.pos = "right"), "ggplot")
})

test_that("view_motifs_lite rejects non-DNA / non-RNA alphabets", {
  aa <- create_motif("ACDEFGHIKL", alphabet = "AA", name = "aa1")
  expect_error(view_motifs_lite(aa), "DNA/RNA")
})

test_that("view_motifs_lite rejects mixed alphabets", {
  dna <- create_motif("ACGT", name = "dna1")
  rna <- create_motif("ACGU", alphabet = "RNA", name = "rna1")
  expect_error(view_motifs_lite(list(dna, rna)), "alphabet")
})

test_that("all motifs must contain the requested multifreq matrix", {
  train <- Biostrings::DNAStringSet(rep("ACGTA", 20))
  m2 <- suppressMessages(create_motif(train, add.multifreq = 2,
                                      pseudocount = 1))
  m3 <- suppressMessages(create_motif(train, add.multifreq = 3,
                                      pseudocount = 1))
  expect_error(view_motifs_lite(list(m2, m3), use.freq = 2),
               "not all motifs have corresponding multifreq matrix",
               fixed = TRUE)
})

test_that("view_motifs_lite dedup.names default rewrites duplicates", {
  m1 <- create_motif("CACGTG", name = "same")
  m2 <- create_motif("TGACGT", name = "same")
  raw <- view_motifs_lite(list(m1, m2), return.raw = TRUE)
  expect_equal(length(unique(names(raw))), 2L)
})

test_that("view_motifs_lite aligns gapped motifs with their gaps filled", {
  ms <- lite_fixture()
  g6 <- add_gap(ms[[6]], gaploc = 4, mingap = 0, maxgap = 3)
  g7 <- add_gap(ms[[7]], gaploc = 3, mingap = 1, maxgap = 2)
  ## Which columns of each aligned motif hold letters, for reference layouts
  ## worked out independently.
  filled <- function(raw) vapply(raw, function(m)
    paste(as.integer(colSums(m) > 0), collapse = ""), character(1))
  raw <- view_motifs_lite(list(g6, ms[[7]], ms[[3]]), return.raw = TRUE,
                          sort.by = "none")
  expect_equal(names(raw), c("M6", "M7 [RC]", "M3 [RC]"))
  expect_equal(unname(filled(raw)), c("0111111111", "0111111110", "1111111110"))
  raw <- view_motifs_lite(list(ms[[1]], g6, g7), return.raw = TRUE,
                          sort.by = "none")
  expect_equal(names(raw), c("M1 [RC]", "M6", "M7 [RC]"))
  expect_equal(unname(filled(raw)), c("0111111111", "1111011111", "1111101110"))
  ## Without their gaps, the motifs are laid out as the ungapped ones are.
  expect_identical(
    view_motifs_lite(list(ms[[1]], g6, g7), return.raw = TRUE, use.gaps = FALSE),
    view_motifs_lite(list(ms[[1]], ms[[6]], ms[[7]]), return.raw = TRUE))
  expect_s3_class(view_motifs_lite(list(ms[[1]], g6, g7)), "ggplot")
})
