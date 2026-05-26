context("CWM motif type")

cwm_matrix <- function() {
  matrix(c( 0.10, -0.20,  0.80, -0.10,
           -0.30,  0.90, -0.20,  0.10,
            0.70, -0.10, -0.10, -0.40,
           -0.20, -0.10, -0.10,  0.90),
         nrow = 4, byrow = FALSE,
         dimnames = list(c("A","C","G","T"), NULL))
}

test_that("create_motif accepts type = 'CWM' with signed values", {
  m <- cwm_matrix()
  cwm <- create_motif(m, type = "CWM", name = "cwm_test")
  expect_equal(cwm@type, "CWM")
  expect_true(any(cwm@motif < 0))
  expect_equal(dim(cwm@motif), c(4L, 4L))
})

test_that("convert_type CWM -> PPM produces a valid PPM", {
  cwm <- create_motif(cwm_matrix(), type = "CWM", name = "x")
  ppm <- convert_type(cwm, "PPM")
  expect_equal(ppm@type, "PPM")
  expect_true(all(abs(colSums(ppm@motif) - 1) < 1e-8))
  expect_true(all(ppm@motif >= 0))
})

test_that("convert_type CWM -> PWM routes via PPM and yields signed log-odds", {
  cwm <- create_motif(cwm_matrix(), type = "CWM", name = "x")
  pwm <- suppressMessages(convert_type(cwm, "PWM"))
  expect_equal(pwm@type, "PWM")
  expect_true(any(pwm@motif < 0))
})

test_that("convert_type PPM -> CWM errors with a clear message", {
  cwm <- create_motif(cwm_matrix(), type = "CWM", name = "x")
  ppm <- convert_type(cwm, "PPM")
  expect_error(convert_type(ppm, "CWM"),
               "CWM cannot be derived from probabilistic motif types")
})

test_that("write_meme / read_meme round-trip with CWM = TRUE preserves values", {
  cwm <- create_motif(cwm_matrix(), type = "CWM", name = "rt_cwm")
  f <- tempfile(fileext = ".meme")
  on.exit(unlink(f), add = TRUE)
  write_meme(cwm, f, CWM = TRUE, overwrite = TRUE)
  re <- read_meme(f, CWM = TRUE)
  re_one <- if (is.list(re)) re[[1]] else re
  expect_equal(re_one@type, "CWM")
  expect_equal(unname(re_one@motif), unname(cwm@motif), tolerance = 1e-6)
})

test_that("write_meme with CWM = FALSE on a CWM converts to PPM first", {
  cwm <- create_motif(cwm_matrix(), type = "CWM", name = "norm_cwm")
  f <- tempfile(fileext = ".meme")
  on.exit(unlink(f), add = TRUE)
  write_meme(cwm, f, CWM = FALSE, overwrite = TRUE)
  re <- read_meme(f)
  re_one <- if (is.list(re)) re[[1]] else re
  expect_equal(re_one@type, "PPM")
  expect_true(all(abs(colSums(re_one@motif) - 1) < 1e-6))
})

test_that("view_motifs and view_motifs2 accept use.type = 'CWM'", {
  cwm <- create_motif(cwm_matrix(), type = "CWM", name = "x")
  expect_s3_class(view_motifs(cwm, use.type = "CWM"), "ggplot")
  expect_s3_class(view_motifs2(cwm, use.type = "CWM"), "ggplot")
})

test_that("use.type = 'CWM' rejects non-CWM motifs", {
  cwm <- create_motif(cwm_matrix(), type = "CWM", name = "x")
  ppm <- convert_type(cwm, "PPM")
  expect_error(view_motifs2(list(cwm, ppm), use.type = "CWM"),
               "type = \"CWM\"")
})

test_that("motif_rc preserves CWM type", {
  cwm <- create_motif(cwm_matrix(), type = "CWM", name = "x")
  rc <- motif_rc(cwm)
  expect_equal(rc@type, "CWM")
  expect_equal(ncol(rc@motif), ncol(cwm@motif))
})

test_that("trim_motifs works on CWM input", {
  m <- cbind(matrix(0.01, nrow = 4, ncol = 2,
                    dimnames = list(c("A","C","G","T"), NULL)),
             cwm_matrix(),
             matrix(0.01, nrow = 4, ncol = 2,
                    dimnames = list(c("A","C","G","T"), NULL)))
  cwm <- create_motif(m, type = "CWM", name = "flanked")
  trimmed <- trim_motifs(cwm)
  expect_equal(trimmed@type, "CWM")
  expect_lt(ncol(trimmed@motif), ncol(cwm@motif))
})
