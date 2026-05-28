context("merge_similar()")

test_that("returns a list of universalmotif objects", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("GGGCCCCC", name = "c")
  r <- merge_similar(list(m1, m2, m3), threshold = 0.7)
  expect_type(r, "list")
  expect_true(all(vapply(r, function(x) is(x, "universalmotif"),
                         logical(1))))
})

test_that("clusters related motifs together; keeps unrelated separate", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("TGACATAT", name = "c")
  m4 <- create_motif("GGGCCCCC", name = "unrelated")
  r <- merge_similar(list(m1, m2, m3, m4), threshold = 0.7)
  expect_true(length(r) < 4L)
  expect_true(length(r) >= 2L)
})

test_that("very strict threshold keeps all motifs separate", {
  m1 <- create_motif("AAAATTTT", name = "a")
  m2 <- create_motif("GGGGCCCC", name = "b")
  m3 <- create_motif("CGCGCGCG", name = "c")
  r <- merge_similar(list(m1, m2, m3), threshold = 0.999)
  expect_equal(length(r), 3L)
})

test_that("single-motif input passes through unchanged", {
  m1 <- create_motif("TTGACATA", name = "solo")
  r <- merge_similar(list(m1))
  if (is.list(r)) {
    expect_equal(length(r), 1L)
    expect_equal(r[[1]]@name, "solo")
  } else {
    expect_s4_class(r, "universalmotif")
    expect_equal(r@name, "solo")
  }
})

test_that("return.clusters = TRUE yields a list of cluster vectors", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  m3 <- create_motif("GGGCCCCC", name = "unrelated")
  r <- merge_similar(list(m1, m2, m3), threshold = 0.7, return.clusters = TRUE)
  expect_type(r, "list")
})
