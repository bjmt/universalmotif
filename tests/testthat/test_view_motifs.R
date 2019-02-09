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

  expect_equal(attributes(r.ppm$data)$class, "waiver")
  expect_equal(length(r.ppm$data), 0)

})

test_that("viewing motif stacks works", {

  m1 <- create_motif(name = "mot1")
  m2 <- create_motif(name = "mot2")

  r <- view_motifs(c(m1, m2))

  expect_is(r, "gg")

})
