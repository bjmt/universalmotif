context("utils")

test_that("misc utilities work", {

  m <- matrix(c(0.99, 0.01, 0, 0, 0.25, 0.25, 0.25, 0.25), ncol = 2)
  m <- create_motif(m, alphabet = "DNA")

  m2 <- round_motif(m, pct.tolerance = 0.05)

  expect_true(m2@motif[1] == 1)
  expect_true(all(m2@motif[2:4] == 0))
  expect_true(all(m2@motif[, 2] == m@motif[, 2]))

})

test_that("summarise_motifs works", {

  m1 <- create_motif("GYNCHGCCKT", name = "motif-1", nsites = 212)
  m2 <- create_motif("GTYCTWYWAK", name = "motif-2", nsites = 147)
  m3 <- create_motif("TAYTAAMCAA", name = "motif-3", nsites = 117)

  res <- summarise_motifs(list(m1, m2, m3))

  expect_true(is.data.frame(res))
  expect_equal(res$name, c("motif-1", "motif-2", "motif-3"))
  expect_equal(res$consensus, c("GYNCHGCCKT", "GTYCTWYWAK", "TAYTAAMCAA"))
  expect_equal(res$alphabet, c("DNA", "DNA", "DNA"))
  expect_equal(res$strand, c("+-", "+-", "+-"))
  expect_equal(round(res$icscore, 2), c(14.42, 15.01, 18.00))
  expect_equal(res$nsites, c(212, 147, 117))

})

test_that("string utilities works", {

  expect_equal(collapse_cpp(c("A", "A", "A")), "AAA")

})

test_that("type utilities work", {

  expect_equal(universalmotif:::pcm_to_ppmC(c(1, 1)), c(0.5, 0.5))
  expect_equal(universalmotif:::pcm_to_ppmC(c(1, 0), 1), c(0.75, 0.25))
  expect_equal(universalmotif:::ppm_to_pcmC(c(0.25, 0.75), 10), c(3, 7))
  expect_equal(round(universalmotif:::ppm_to_pwmC(c(0.25, 0.75), numeric()), 3),
               c(-1.000, 0.585))
  expect_equal(round(universalmotif:::ppm_to_pwmC(c(0.25, 0.75), c(0.4, 0.6)), 3),
               c(-0.678, 0.322))
  expect_equal(round(universalmotif:::pwm_to_ppmC(c(-1, 0.5), numeric()), 4),
               c(0.2612, 0.7388))
  expect_equal(round(universalmotif:::ppm_to_icmC(c(0.25, 0.75), numeric()), 4),
               c(0.0472, 0.1415))
  expect_equal(round(universalmotif:::ppm_to_icmC(c(0.25, 0.75), numeric(), TRUE), 4),
               c(0.0000, 0.4387))
  expect_equal(round(universalmotif:::position_icscoreC(c(0.25, 0.75), numeric()), 3), 0.185)
  expect_equal(universalmotif:::icm_to_ppmC(c(0.5, 0.0, 0.3)),
               c(0.625, 0.000, 0.375))

})
