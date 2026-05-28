context("trim_motifs()")

test_that("motif trimming works", {

  m1 <- create_motif("NNCCCNN", nsites = 100)
  m2 <- trim_motifs(m1)

  expect_equal(ncol(m2), 3)

})

test_that("trim_motifs() preserves successfully-trimmed motifs when one is fully trimmed (regression: return(NULL) on partial trim)", {

  m1 <- create_motif("ACGTACGT", pseudocount = 1, nsites = 100)
  m2 <- create_motif(matrix(rep(0.25, 16), nrow = 4,
                             dimnames = list(c("A", "C", "G", "T"), NULL)),
                     alphabet = "DNA", type = "PPM")
  out <- suppressMessages(trim_motifs(c(m1, m2), min.ic = 0.5))
  out_list <- if (is.list(out)) out else list(out)
  expect_true(length(out_list) >= 1)
  expect_equal(out_list[[1]]@name, m1@name)

})

test_that("min.ic = 0 keeps every column", {
  m <- create_motif("NNCCCNN", nsites = 100)
  trimmed <- trim_motifs(m, min.ic = 0)
  expect_equal(ncol(trimmed), ncol(m))
})

test_that("min.ic above all column ICs errors with a clear message", {
  ## With min.ic well above the maximum possible IC (2 for DNA), every column
  ## should be trimmed away; the documented behaviour is to error rather than
  ## return a zero-column motif.
  m <- create_motif("NNCCCNN", nsites = 100)
  expect_error(suppressMessages(trim_motifs(m, min.ic = 5)),
                regexp = "completely trimmed")
})

test_that("trim_motifs() trims only from the left when trim.from = 'left'", {
  m <- create_motif("NNCCCNN", nsites = 100)
  trimmed <- trim_motifs(m, min.ic = 0.5, trim.from = "left")
  expect_lte(ncol(trimmed), ncol(m))
  ## Right-edge low-IC columns are kept when only trimming from the left.
  expect_gt(ncol(trimmed), 3)
})

test_that("trim_motifs() trims only from the right when trim.from = 'right'", {
  m <- create_motif("NNCCCNN", nsites = 100)
  trimmed <- trim_motifs(m, min.ic = 0.5, trim.from = "right")
  expect_lte(ncol(trimmed), ncol(m))
  expect_gt(ncol(trimmed), 3)
})
