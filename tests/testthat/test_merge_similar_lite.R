context("merge_similar_lite()")

test_that("returns a list of universalmotif S4 objects", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("GGGCCCCC", name = "c")
  r <- merge_similar_lite(list(m1, m2, m3), qvalue = 0.05)
  expect_type(r, "list")
  expect_true(all(vapply(r, function(x) is(x, "universalmotif"),
                         logical(1))))
})

test_that("clusters related motifs together; keeps unrelated separate", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  m4 <- create_motif("GGGCCCCC", name = "unrelated")
  expect_warning(r <- merge_similar_lite(list(m1, m2, m3, m4), qvalue = 0.05),
                 "cluster 1 holds 3 of the 4 motifs")
  expect_equal(length(r), 2L)
})

test_that("no-similarity input returns same number of motifs", {
  m1 <- create_motif("AAAATTTT", name = "a")
  m2 <- create_motif("GGGGCCCC", name = "b")
  m3 <- create_motif("CGCGCGCG", name = "c")
  r <- merge_similar_lite(list(m1, m2, m3), qvalue = 1e-30)
  expect_equal(length(r), 3L)
})

test_that("very loose qvalue collapses all motifs into one cluster", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  expect_warning(r <- merge_similar_lite(list(m1, m2, m3), qvalue = 1),
                 "unlike motifs were joined")
  expect_equal(length(r), 1L)
})

test_that("single-motif input passes through unchanged", {
  m1 <- create_motif("TTGACATA", name = "solo")
  r <- merge_similar_lite(list(m1))
  expect_equal(length(r), 1L)
  expect_equal(r[[1]]@name, "solo")
})

test_that("empty input is handled cleanly", {
  expect_equal(length(merge_similar_lite(list())), 0L)
  df <- merge_similar_lite(list(), return.clusters = TRUE)
  expect_s3_class(df, "universalmotif_df")
  expect_equal(nrow(df), 0L)
  expect_true(all(c("motif", "name", "motif.i", "cluster") %in% names(df)))
})

test_that("return.clusters = TRUE yields a universalmotif_df", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  m4 <- create_motif("GGGCCCCC", name = "unrelated")
  df <- suppressWarnings(merge_similar_lite(list(m1, m2, m3, m4), qvalue = 0.05,
                                            return.clusters = TRUE))
  expect_s3_class(df, "universalmotif_df")
  expect_true(all(c("motif", "name", "motif.i", "cluster") %in% names(df)))
  expect_equal(nrow(df), 4L)
  expect_equal(df$motif.i, 1:4)
  expect_equal(df$name, c("a", "b", "c", "unrelated"))
  expect_true(is(df$motif[[1]], "universalmotif"))
  expect_setequal(df$cluster, c(1L, 2L))
})

test_that("determinism: same input twice -> identical result", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("GGGCCCCC", name = "c")
  r1 <- merge_similar_lite(list(m1, m2, m3), qvalue = 0.05)
  r2 <- merge_similar_lite(list(m1, m2, m3), qvalue = 0.05)
  expect_equal(length(r1), length(r2))
  expect_equal(vapply(r1, function(x) x@consensus, character(1)),
               vapply(r2, function(x) x@consensus, character(1)))
})

test_that("AA motifs are rejected", {
  m1 <- create_motif("AAYY", alphabet = "AA", name = "p1")
  m2 <- create_motif("AAYY", alphabet = "AA", name = "p2")
  expect_error(merge_similar_lite(list(m1, m2)), regexp = "DNA/RNA")
})

test_that("null = parametric also runs", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("GGGCCCCC", name = "c")
  r <- merge_similar_lite(list(m1, m2, m3), qvalue = 0.05, null = "parametric")
  expect_type(r, "list")
  expect_true(length(r) >= 1L && length(r) <= 3L)
})

# ---- Linkage ------------------------------------------------------------------
#
# Reference values for the fixture motifs of helper-lite-fixture.R were worked
# out independently. M1, M2 and M3 form a chain at q <= 0.01: M1 and M3 are
# both linked to M2 but not to each other.

cluster_warnings <- function(expr) {
  w <- character()
  val <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
  list(value = val, warnings = w)
}

test_that("complete linkage agrees with the reference values", {
  ms <- lite_fixture()
  ref <- list(
    "0.01" = list(single = c(1, 1, 1, 2, 3, 4, 5), complete = c(1, 1, 2, 3, 4, 5, 6)),
    "0.03" = list(single = c(1, 1, 1, 2, 3, 4, 4), complete = c(1, 1, 1, 2, 3, 4, 4)),
    "0.06" = list(single = c(1, 1, 1, 2, 3, 1, 1), complete = c(1, 1, 1, 2, 3, 4, 4)),
    "0.2"  = list(single = c(1, 1, 1, 2, 1, 1, 1), complete = c(1, 1, 1, 2, 3, 4, 4)))
  for (q in names(ref)) for (link in c("single", "complete")) {
    res <- suppressWarnings(merge_similar_lite(ms, qvalue = as.numeric(q),
      linkage = link, return.clusters = TRUE))
    expect_equal(res$cluster, as.integer(ref[[q]][[link]]), info = paste(q, link))
  }
  merged <- merge_similar_lite(ms, qvalue = 0.01, linkage = "complete")
  expect_equal(vapply(merged, function(m) m@name, character(1)),
               c("merged_M1+M2", "M3", "M4", "M5", "M6", "M7"))
  expect_equal(merged[[1]]@consensus, "CTTGACATAG")
  merged <- suppressWarnings(merge_similar_lite(ms, qvalue = 0.06))
  expect_equal(merged[[1]]@name, "merged_M1+M2+M3+M6+M7")
  expect_equal(merged[[1]]@consensus, "SMTRWGTCAWG")
})

test_that("single linkage is the default", {
  ms <- lite_fixture()
  a <- suppressWarnings(merge_similar_lite(ms, qvalue = 0.06, return.clusters = TRUE))
  b <- suppressWarnings(merge_similar_lite(ms, qvalue = 0.06, return.clusters = TRUE,
                                           linkage = "single"))
  expect_identical(a, b)
  expect_error(merge_similar_lite(ms, linkage = "average"), "should be one of")
})

test_that("a cluster holding over a quarter of the motifs brings a warning", {
  ms <- lite_fixture()
  w <- cluster_warnings(merge_similar_lite(ms, qvalue = 0.06))$warnings
  expect_length(w, 1L)
  expect_match(w, "cluster 1 holds 5 of the 7 motifs")
  expect_match(w, "complete linkage")
  w <- cluster_warnings(merge_similar_lite(ms, qvalue = 0.06, linkage = "complete",
                                           return.clusters = TRUE))$warnings
  expect_length(w, 1L)
  expect_match(w, "cluster 1 holds 3 of the 7 motifs")
  expect_no_match(w, "complete linkage")
  ## Clusters of two never bring it, however few the motifs.
  expect_silent(merge_similar_lite(ms, qvalue = 0.01, linkage = "complete"))
  expect_silent(merge_similar_lite(ms[c(1, 2, 4)], qvalue = 0.01))
})

test_that("gapped motifs are clustered and merged with their gaps", {
  ms <- lite_fixture()
  x <- c(ms[1:5],
         list(add_gap(ms[[6]], gaploc = 4, mingap = 0, maxgap = 3),
              add_gap(ms[[7]], gaploc = 3, mingap = 1, maxgap = 2)))
  res <- cluster_warnings(merge_similar_lite(x, qvalue = 0.3))
  expect_equal(vapply(res$value, function(m) m@name, character(1)),
               c("merged_M1+M2+M3+M5+M6", "M4", "M7"))
  expect_equal(res$value[[1]]@consensus, "NHTATGTCANGTG")
  expect_equal(res$value[[3]]@gapinfo@gaploc, 3)
  expect_length(res$warnings, 2L)
  expect_match(res$warnings[2],
               "^cluster 1: the gaps of M6 were filled with blank columns")
  cl <- suppressWarnings(merge_similar_lite(x, qvalue = 0.3, return.clusters = TRUE))
  expect_equal(cl$cluster, c(1L, 1L, 1L, 2L, 1L, 1L, 3L))
  cl <- suppressWarnings(merge_similar_lite(x, qvalue = 0.05, return.clusters = TRUE))
  expect_equal(cl$cluster, c(1L, 1L, 1L, 2L, 3L, 4L, 5L))
  ## Without their gaps, the motifs cluster as the ungapped fixture does.
  expect_identical(
    suppressWarnings(merge_similar_lite(x, qvalue = 0.06, use.gaps = FALSE,
                                        return.clusters = TRUE))$cluster,
    suppressWarnings(merge_similar_lite(ms, qvalue = 0.06,
                                        return.clusters = TRUE))$cluster)
})
