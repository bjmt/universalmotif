context("geom_motif() / geom_logo()")

suppressWarnings(library(ggplot2))

dna_cols <- c(A = "#109648", C = "#255C99", G = "#F7B32B", T = "#D62839")

test_that("annotate_motif adds a layer of rescaled polygons inside the box", {
  data(examplemotif)
  p <- ggplot(data.frame(x = 1:10, y = 1:10), aes(x, y)) +
    geom_point() +
    annotate_motif(examplemotif, xmin = 2, xmax = 8, ymin = 3, ymax = 6)
  expect_s3_class(p, "ggplot")
  expect_length(p$layers, 2L)
  d <- layer_data(p, 2L)
  expect_true(all(c("x", "y", "group", "fill_colour") %in% names(d)))
  expect_gt(nrow(d), 0L)
  expect_true(all(d$x >= 2 - 1e-6 & d$x <= 8 + 1e-6))
  expect_true(all(d$y >= 3 - 1e-6 & d$y <= 6 + 1e-6))
})

test_that("geom_motif works across alphabets and use.types", {
  for (a in c("DNA", "RNA", "AA", "QWERTY")) {
    m <- create_motif(alphabet = a)
    p <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) +
      geom_motif(motif = m, width = 1, height = 2)
    expect_silent(d <- layer_data(p, 1L))
    expect_gt(nrow(d), 0L)
  }
  m <- create_motif("TATAAA")
  for (ty in c("PPM", "PCM", "PWM", "ICM")) {
    p <- ggplot() + annotate_motif(m, 1, 6, 0, 2, use.type = ty)
    expect_gt(nrow(layer_data(p, 1L)), 0L)
  }
})

test_that("anchored mode places the logo around x/y", {
  m <- create_motif("ACGT", name = "m")          # ncol == 4
  p <- ggplot(data.frame(x = 10, y = 5), aes(x, y)) +
    geom_motif(motif = m, width = 1, height = 2)
  d <- layer_data(p, 1L)
  expect_true(all(d$x >= 10 - 1e-6 & d$x <= 14 + 1e-6))
  expect_true(all(d$y >= 4 - 1e-6 & d$y <= 6 + 1e-6))
})

test_that("geom_motif selects a different motif per group from a named list", {
  m1 <- create_motif("TTTT", name = "a")
  m2 <- create_motif("GGGG", name = "b")
  df <- data.frame(x = 0, y = c(1, 2), label = c("a", "b"))
  p <- ggplot(df, aes(x = x, y = y, motif = label)) +
    geom_motif(motif = list(a = m1, b = m2), height = 0.8)
  d <- layer_data(p, 1L)
  expect_gt(nrow(d), 0L)
  expect_gte(length(unique(d$group)), 2L)
})

test_that("annotate_logo / geom_logo accept a bare row-named matrix", {
  mat <- matrix(c(0.7, 0.1, 0.1, 0.1, 0.1, 0.7, 0.1, 0.1), nrow = 4,
                dimnames = list(c("A", "C", "G", "T"), NULL))
  p <- ggplot() +
    annotate_logo(mat, xmin = 1, xmax = 3, ymin = 0, ymax = 1,
                  colour.scheme = dna_cols)
  expect_gt(nrow(layer_data(p, 1L)), 0L)

  ## Arbitrary / negative values are tolerated.
  mat2 <- mat
  mat2[1, 1] <- -0.5
  p2 <- ggplot() + annotate_logo(mat2, 1, 3, -1, 1)
  expect_gt(nrow(layer_data(p2, 1L)), 0L)
})

test_that("geom_logo errors on a matrix without row names", {
  mat <- matrix(c(1, 0, 0, 0, 0, 1, 0, 0), nrow = 4)
  p <- ggplot() + annotate_logo(mat, 1, 2, 0, 1)
  expect_error(ggplot_build(p), "row names")
})

test_that("colour.scheme must cover every letter", {
  data(examplemotif)
  p <- ggplot() +
    annotate_motif(examplemotif, 1, 5, 0, 2, colour.scheme = c(A = "red"))
  expect_error(ggplot_build(p), "colour.scheme")
})

test_that("an unknown motif key errors clearly", {
  m1 <- create_motif("TTTT", name = "a")
  df <- data.frame(x = 0, y = 1, label = "zzz")
  p <- ggplot(df, aes(x = x, y = y, motif = label)) +
    geom_motif(motif = list(a = m1))
  expect_error(ggplot_build(p), "motif")
})

test_that("supplying both a box and an anchor errors", {
  m <- create_motif("ACGT")
  df <- data.frame(x = 0, y = 0, xmin = 1, xmax = 2, ymin = 0, ymax = 1)
  p <- ggplot(df, aes(x = x, y = y, xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax)) +
    geom_motif(motif = m)
  expect_error(ggplot_build(p), "both")
})

test_that("a zero-information motif draws nothing without warning or error", {
  uniform <- create_motif(
    matrix(rep(0.25, 24), nrow = 4, dimnames = list(c("A", "C", "G", "T"), NULL)),
    alphabet = "DNA", type = "PPM")
  good <- create_motif("TATAAA", name = "good")
  ## Single uniform motif: empty, but builds cleanly.
  p <- ggplot(data.frame(x = 0, y = 0), aes(x, y)) + geom_motif(motif = uniform)
  expect_silent(d <- layer_data(p, 1L))
  expect_equal(nrow(d), 0L)
  ## Mixed list: the informative motif is still drawn for its row.
  df <- data.frame(x = 0, y = c(1, 2), label = c("u", "g"))
  p2 <- ggplot(df, aes(x, y, motif = label)) +
    geom_motif(motif = list(u = uniform, g = good))
  expect_silent(d2 <- layer_data(p2, 1L))
  expect_gt(nrow(d2), 0L)
})

test_that("align_motif_mats returns equal-width, named matrices", {
  m1 <- create_motif("TTGACATA", name = "a")
  m2 <- create_motif("CTTGACAT", name = "b")
  res <- universalmotif:::align_motif_mats(list(m1, m2))
  expect_type(res, "list")
  expect_length(res, 2L)
  w <- vapply(res, ncol, integer(1))
  expect_equal(length(unique(w)), 1L)
  expect_equal(names(res), c("a", "b"))
})
