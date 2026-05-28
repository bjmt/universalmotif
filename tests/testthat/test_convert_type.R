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

test_that("PPM -> ICM at nsites = 0 falls back to default and stays finite", {
  m <- create_motif("ACGTACGT", nsites = 100)
  m@nsites <- numeric(0)
  icm <- convert_type(m, "ICM")
  expect_true(all(is.finite(icm@motif)))
  expect_equal(icm@type, "ICM")
})

test_that("PPM -> ICM with nsize_correction at nsites = 1 produces finite values", {
  m <- create_motif("ACGT", nsites = 1)
  icm <- convert_type(m, "ICM", nsize_correction = TRUE)
  expect_true(all(is.finite(icm@motif)))
})

test_that("high pseudocount smears PCM -> PPM toward uniform", {
  m <- create_motif("ACGTACGT", nsites = 10)
  m_pcm <- convert_type(m, "PCM")
  m_pcm@pseudocount <- 1000
  m_ppm <- convert_type(m_pcm, "PPM")
  ## With pseudocount >> nsites, columns should be near-uniform (0.25 each).
  expect_true(all(abs(m_ppm@motif - 0.25) < 0.05))
})

test_that("PPM -> PCM -> PPM round-trips at non-default nsites", {
  m <- create_motif("TWWWWWWW", nsites = 50)
  m_pcm <- convert_type(m, "PCM")
  m_back <- convert_type(m_pcm, "PPM")
  expect_equal(m@motif, m_back@motif, tolerance = 1e-3)
})

test_that("PWM -> PPM round-trips when bkg is uniform", {
  m <- create_motif("AAAATTTC", nsites = 100,
                    bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25))
  m_pwm <- convert_type(m, "PWM")
  m_back <- convert_type(m_pwm, "PPM")
  expect_equal(m@motif, m_back@motif, tolerance = 1e-3)
})
