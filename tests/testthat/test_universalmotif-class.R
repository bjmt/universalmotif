library("universalmotif")
context("Test universalmotif-class related methods")

test_pcm <- matrix(c( 0,  5, 10, 15, 20,  4,  6, 11, 14,  9,
                      5, 10,  5,  2,  0,  6,  8,  9,  3,  1,
                      9,  4,  1,  2,  0, 10,  6,  0,  1,  0,
                      6,  1,  4,  1,  0,  0,  0,  0,  2, 10),
                   byrow = TRUE, nrow = 4)
umot_test_pcm <- new("universalmotif", name = "testmotif", motif = test_pcm,
                     type = "PCM")
umot_test_ppm <- convert_type(umot_test_pcm, "PPM")
umot_test_pwm <- convert_type(umot_test_pcm, "PWM")

test_that("initalized slots are correct", {
  expect_true(all(motif_slots(umot_test_pcm, "consensus") == 
                   c("N", "N", "N", "A", "A", "S", "V", "M", "A", "W")))
  expect_equal(motif_slots(umot_test_pcm, "icscore"), c("icscore" = 6.473),
               tolerance = 0.00009)
})

test_that("type manipulation functions are correct", {
  expect_identical(convert_type(umot_test_ppm, "PCM")@motif, umot_test_pcm@motif)
  expect_identical(convert_type(umot_test_pwm, "PCM")@motif, umot_test_pcm@motif)
})
