library("universalmotif")
context("Test universalmotif-class related methods")

test_pcm <- matrix(c( 0,  5, 10, 15, 20,  4,  6, 11, 14,  9,
                      5, 10,  5,  2,  0,  6,  8,  9,  3,  1,
                      9,  4,  1,  2,  0, 10,  6,  0,  1,  0,
                      6,  1,  4,  1,  0,  0,  0,  0,  2, 10),
                   byrow = TRUE, nrow = 4)
test_ppm <- matrix(c(0.00, 0.25, 0.50, 0.75, 1, 0.2, 0.3, 0.55, 0.70, 0.45,
                     0.25, 0.50, 0.25, 0.10, 0, 0.3, 0.4, 0.45, 0.15, 0.05,
                     0.45, 0.20, 0.05, 0.10, 0, 0.5, 0.3, 0.00, 0.05, 0.00,
                     0.30, 0.05, 0.20, 0.05, 0, 0.0, 0.0, 0.00, 0.10, 0.50),
                   byrow = TRUE, nrow = 4)
umot_test_pcm <- universalmotif(name = "testmotif", motif = test_pcm)
umot_test_ppm <- universalmotif(name = "testmotif", motif = test_ppm,
                                nsites = 20) 
umot_test_ppm_from_pcm <- convert_type(umot_test_pcm, "PPM")
umot_test_pwm_from_pcm <- convert_type(umot_test_pcm, "PWM", pseudoweight = 0)
umot_test_pwm_from_ppm <- convert_type(umot_test_ppm, "PWM", pseudoweight = 0)
umot_test_ppm_from_pwm <- convert_type(umot_test_pwm_from_ppm, "PPM",
                                       pseudoweight = 0)
umot_test_icm_from_pcm <- convert_type(umot_test_pcm, "ICM")
umot_test_icm2_from_pcm <- convert_type(umot_test_pcm, "ICM",
                                       background = c(0.2, 0.3, 0.2, 0.3))

test_that("initalized slots are correct", {
  expect_true(all(motif_slots(umot_test_pcm, "consensus") == "NNNAASVMAW"))
  expect_equal(motif_slots(umot_test_pcm, "icscore"), c("icscore" = 6.473),
               tolerance = 0.00009)
  expect_identical(motif_slots(umot_test_pcm, "icscore"),
                   motif_slots(umot_test_ppm, "icscore"))
})

# the following tests are extremely important!
test_that("type manipulation functions are correct", {
  expect_identical(convert_type(umot_test_ppm_from_pcm, "PCM")@motif,
                   umot_test_pcm@motif)
  expect_identical(convert_type(umot_test_pwm_from_pcm, "PCM")@motif,
                   umot_test_pcm@motif)
  expect_identical(umot_test_pwm_from_pcm, umot_test_pwm_from_ppm)
  expect_identical(umot_test_ppm_from_pwm, umot_test_ppm)
  expect_true(all(umot_test_icm_from_pcm@motif <= 2))
  expect_true(any(umot_test_icm2_from_pcm@motif > 2))
})
