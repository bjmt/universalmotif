context("motif_tree2()")

test_that("motif_tree2 returns a ggtree on a motif list", {
  skip_if_not_installed("ape")
  skip_if_not_installed("ggtree")
  motifs <- lapply(c("CACGTG", "TGACGT", "AAATTT", "GCGCGC", "TTGACATA",
                     "CTTGACAT"),
                   function(s) create_motif(s, name = s))
  tr <- suppressWarnings(motif_tree2(motifs, linecol = "none"))
  expect_s3_class(tr, "ggplot")
  expect_s3_class(tr, "ggtree")
})

test_that("motif_tree2 draws tip logos by default on a rectangular tree", {
  skip_if_not_installed("ape")
  skip_if_not_installed("ggtree")
  motifs <- lapply(c("CACGTG", "TGACGT", "AAATTT", "GCGCGC", "TTGACATA",
                     "CTTGACAT"),
                   function(s) create_motif(s, name = s))
  with_logo <- suppressWarnings(motif_tree2(motifs, linecol = "none"))
  no_logo <- suppressWarnings(motif_tree2(motifs, linecol = "none",
                                          tiplogo = FALSE))
  expect_s3_class(with_logo, "ggtree")
  ## The tip-logo layer is the one extra layer.
  expect_equal(length(with_logo$layers), length(no_logo$layers) + 1L)
  d <- suppressWarnings(
    ggplot2::layer_data(with_logo, length(with_logo$layers)))
  expect_gt(nrow(d), 0L)
})

test_that("motif_tree2 warns and skips tip logos on non-linear layouts", {
  skip_if_not_installed("ape")
  skip_if_not_installed("ggtree")
  motifs <- lapply(c("CACGTG", "TGACGT", "AAATTT", "GCGCGC"),
                   function(s) create_motif(s, name = s))
  expect_warning(
    p <- motif_tree2(motifs, linecol = "none", layout = "circular"),
    "rectangular"
  )
  expect_s3_class(p, "ggtree")
})

test_that("motif_tree2 short-circuits on dist input", {
  skip_if_not_installed("ape")
  skip_if_not_installed("ggtree")
  motifs <- lapply(c("CACGTG", "TGACGT", "AAATTT", "GCGCGC"),
                   function(s) create_motif(s, name = s))
  score.mat <- compare_motifs2(motifs, matrix.out = "score")
  d <- as.dist((1 - score.mat) / 2)
  tr <- suppressWarnings(motif_tree2(d, labels = "name"))
  expect_s3_class(tr, "ggtree")
})

test_that("motif_tree2 distance values land in [0, 1]", {
  motifs <- lapply(c("CACGTG", "TGACGT", "AAATTT", "GCGCGC"),
                   function(s) create_motif(s, name = s))
  score.mat <- compare_motifs2(motifs, matrix.out = "score")
  dist.mat <- (1 - score.mat) / 2
  expect_true(all(dist.mat >= 0))
  expect_true(all(dist.mat <= 1))
  ## Self-similarity is 1, so self-distance is 0.
  expect_true(all(diag(dist.mat) < 1e-8))
})

test_that("motif_tree2 with linecol = 'name' colours branches", {
  skip_if_not_installed("ape")
  skip_if_not_installed("ggtree")
  motifs <- lapply(c("CACGTG", "TGACGT", "AAATTT", "GCGCGC", "TTGACATA"),
                   function(s) create_motif(s, name = s))
  ## Stamp distinct family values so linecol has something to colour by.
  fams <- c("ebox", "ebox-like", "AT", "GC", "creb")
  for (i in seq_along(motifs)) motifs[[i]]["family"] <- fams[i]
  tr <- suppressWarnings(motif_tree2(motifs, linecol = "family",
                                     labels = "name",
                                     layout = "rectangular"))
  expect_s3_class(tr, "ggtree")
})

test_that("motif_tree2 rejects illegal layout", {
  motifs <- list(create_motif("ACGT"), create_motif("TGCA"))
  expect_error(motif_tree2(motifs, layout = "not-a-layout"), "layout")
})

test_that("motif_tree and motif_tree2 both run on the same DNA fixture (v1/v2 parity smoke)", {
  skip_if_not_installed("ape")
  skip_if_not_installed("ggtree")
  motifs <- list(create_motif("ACGTAC", name = "a"),
                 create_motif("ACGTAC", name = "b"),
                 create_motif("GGGCCC", name = "c"))
  ## Both should produce a non-null result; exact shape (ggplot vs. tree
  ## object) differs by version, so only existence is checked.
  t1 <- tryCatch(motif_tree(motifs), error = function(e) NULL)
  t2 <- tryCatch(motif_tree2(motifs), error = function(e) NULL)
  expect_false(is.null(t1))
  expect_false(is.null(t2))
})
