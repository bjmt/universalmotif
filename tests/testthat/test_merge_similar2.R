context("merge_similar2()")

test_that("returns a list of universalmotif S4 objects", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("GGGCCCCC", name = "c")
  r <- merge_similar2(list(m1, m2, m3), qvalue = 0.05)
  expect_type(r, "list")
  expect_true(all(vapply(r, function(x) is(x, "universalmotif"),
                         logical(1))))
})

test_that("clusters related motifs together; keeps unrelated separate", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  m4 <- create_motif("GGGCCCCC", name = "unrelated")
  r <- merge_similar2(list(m1, m2, m3, m4), qvalue = 0.05)
  expect_equal(length(r), 2L)
})

test_that("no-similarity input returns same number of motifs", {
  m1 <- create_motif("AAAATTTT", name = "a")
  m2 <- create_motif("GGGGCCCC", name = "b")
  m3 <- create_motif("CGCGCGCG", name = "c")
  r <- merge_similar2(list(m1, m2, m3), qvalue = 1e-30)
  expect_equal(length(r), 3L)
})

test_that("very loose qvalue collapses all motifs into one cluster", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  r <- merge_similar2(list(m1, m2, m3), qvalue = 1)
  expect_equal(length(r), 1L)
})

test_that("single-motif input passes through unchanged", {
  m1 <- create_motif("TTGACATA", name = "solo")
  r <- merge_similar2(list(m1))
  expect_equal(length(r), 1L)
  expect_equal(r[[1]]@name, "solo")
})

test_that("empty input is handled cleanly", {
  expect_equal(length(merge_similar2(list())), 0L)
  df <- merge_similar2(list(), return.clusters = TRUE)
  expect_s3_class(df, "universalmotif_df")
  expect_equal(nrow(df), 0L)
  expect_true(all(c("motif", "name", "motif.i", "cluster") %in% names(df)))
})

test_that("return.clusters = TRUE yields a universalmotif_df", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  m4 <- create_motif("GGGCCCCC", name = "unrelated")
  df <- merge_similar2(list(m1, m2, m3, m4), qvalue = 0.05,
                       return.clusters = TRUE)
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
  r1 <- merge_similar2(list(m1, m2, m3), qvalue = 0.05)
  r2 <- merge_similar2(list(m1, m2, m3), qvalue = 0.05)
  expect_equal(length(r1), length(r2))
  expect_equal(vapply(r1, function(x) x@consensus, character(1)),
               vapply(r2, function(x) x@consensus, character(1)))
})

test_that("AA motifs are rejected", {
  m1 <- create_motif("AAYY", alphabet = "AA", name = "p1")
  m2 <- create_motif("AAYY", alphabet = "AA", name = "p2")
  expect_error(merge_similar2(list(m1, m2)), regexp = "DNA/RNA")
})

test_that("null = parametric also runs", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("GGGCCCCC", name = "c")
  r <- merge_similar2(list(m1, m2, m3), qvalue = 0.05, null = "parametric")
  expect_type(r, "list")
  expect_true(length(r) >= 1L && length(r) <= 3L)
})
