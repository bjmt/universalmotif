context("view_logo()")

test_that("view_logo returns a ggplot with the expected data columns", {
  data(examplemotif)
  mat <- examplemotif["motif"]
  p <- view_logo(mat)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("x", "y", "letter.id", "group") %in% names(p$data)))
  ## One polygon path per letter glyph; data should be non-empty.
  expect_gt(nrow(p$data), 0L)
})

test_that("view_logo errors on a matrix without row names", {
  mat <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 4)
  expect_error(view_logo(mat), regexp = "row names")
})

test_that("view_motifs(return.raw = TRUE) returns a list of matrices", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  res <- view_motifs(list(m1, m2), return.raw = TRUE)
  expect_type(res, "list")
  expect_true(all(vapply(res, is.matrix, logical(1))) ||
              all(vapply(res, is.list, logical(1))))
})
