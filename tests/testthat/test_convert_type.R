context("convert_type()")

test_that("type conversion works", {

  m.ppm <- create_motif("TWWWWWWW", nsites = 10)
  m.pwm <- convert_type(m.ppm, "PWM")
  m.icm <- convert_type(m.ppm, "ICM")
  m.pcm <- convert_type(m.ppm, "PCM")

  expect_equal(m.ppm@motif[1, 2], 0.5)
  expect_equal(m.icm@motif[4, 1], 2)
  expect_equal(m.pwm@motif[2, 2], -Inf)
  expect_equal(m.pcm@motif[1, 2], 5)

  m.icm2 <- convert_type(m.pcm, "ICM")
  # m.icm3 <- convert_type(m.pcm, "ICM", nsize_correction = TRUE)
  m.icm4 <- convert_type(m.pwm, "ICM")

  expect_equal(m.icm2@motif[1, 2], 0.5)
  # expect_equal(round(m.icm3@motif[1, 2], digits = 2), 0.38)
  expect_equal(m.icm4@motif[1, 2], 0.5)

  m.ppm2 <- convert_type(m.icm, "PPM")
  m.pcm2 <- convert_type(m.icm, "PCM")
  m.pwm2 <- convert_type(m.icm, "PWM")

  expect_equal(m.ppm2@motif[1, 2], 0.5)
  expect_equal(m.pcm2@motif[1, 2], 5)
  expect_equal(m.pwm2@motif[2, 2], -Inf)

})
