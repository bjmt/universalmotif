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

})
