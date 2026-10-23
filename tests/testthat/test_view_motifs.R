context("view_motifs()")

test_that("viewing single motifs works", {

  m.ppm <- create_motif(type = "PPM")
  m.pcm <- create_motif(type = "PCM")
  m.pwm <- create_motif(type = "PWM")
  m.icm <- create_motif(type = "ICM")
  m.rna <- create_motif(alphabet = "RNA")
  m.aa <- create_motif(alphabet = "AA")
  m.custom <- create_motif(alphabet = "QWERTY")

  r.ppm <- view_motifs(m.ppm, use.type = "PPM")
  r.pcm <- view_motifs(m.pcm, use.type = "PCM")
  r.pwm <- view_motifs(m.pwm, use.type = "PWM")
  r.icm <- view_motifs(m.icm, use.type = "ICM")
  r.rna <- view_motifs(m.rna)
  r.aa <- view_motifs(m.aa)
  r.custom <- view_motifs(m.custom)

  expect_is(r.ppm, "gg")
  expect_is(r.pcm, "gg")
  expect_is(r.pwm, "gg")
  expect_is(r.icm, "gg")
  expect_is(r.rna, "gg")
  expect_is(r.aa, "gg")
  expect_is(r.custom, "gg")

  # expect_equal(attributes(r.ppm$data)$class, "waiver")
  # expect_equal(length(r.ppm$data), 0)

})

test_that("viewing motif stacks works", {

  m1 <- create_motif(name = "mot1")
  m2 <- create_motif(name = "mot2")

  r <- view_motifs(c(m1, m2))

  expect_is(r, "gg")

})

test_that("all motifs must contain the requested multifreq matrix", {
  train <- Biostrings::DNAStringSet(rep("ACGTA", 20))
  m2 <- suppressMessages(create_motif(train, add.multifreq = 2,
                                      pseudocount = 1))
  m3 <- suppressMessages(create_motif(train, add.multifreq = 3,
                                      pseudocount = 1))
  expect_error(view_motifs(list(m2, m3), use.freq = 2),
               "not all motifs have corresponding multifreq matrix",
               fixed = TRUE)
})

test_that("logo axes keep their top break and suit CWMs", {
  yb <- function(p) ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$breaks
  expect_equal(yb(view_motifs(create_motif(6, alphabet = "QWERTY"))), c(0, 1, 2))
  expect_equal(signed_logo_axis(-3.2, 7.6, "PWM")$breaks, c(-3.3, -2, 0, 4, 7))
  cwm <- signed_logo_axis(-0.02, 0.004, "CWM")
  expect_equal(cwm$breaks, c(-0.02, -0.01, 0))
  expect_equal(cwm$limits, c(-0.02 - 0.05 * 0.024, 0.004))
})

test_that("small negative letters are drawn with the same gaps as positive ones", {
  font <- as.data.frame(fontDFroboto)
  d <- make_position_polygon_data(c(A = 0.5, C = 0.3, G = -0.4, T = -0.2),
    fontDF = font, min.height = 0.01, y.spacer = 0.05, position.id = 1)
  ext <- tapply(d$y, d$letter.id, range)
  expect_equal(ext[["1.A.1"]][1] - ext[["1.C.1"]][2], 0.025)
  expect_equal(ext[["1.G.1"]][1] - ext[["1.T.1"]][2], 0.025)
  expect_equal(ext[["1.T.1"]][1], -0.6)
  d <- make_position_polygon_data(c(A = 0.004, C = -0.02, G = 0, T = -0.012),
    fontDF = font, min.height = 0.001, position.id = 1)
  expect_true(all(c("1.C.1", "1.T.1") %in% d$letter.id))
})
