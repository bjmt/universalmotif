context("Test motif peak finding")

test_that("density and Pvalue calculations work", {

  hits <- c(rep(10, 200), rep(15, 200), rep(20, 200))
  res <- motif_peaks(hits, 1000, 50)

  expect_true(all(res$Peaks$Peak == c(10, 15, 20)))

  res <- motif_peaks(hits, 1000, 50, 100)
  expect_true(all(res$Peaks$Peak == 15))

  res <- motif_peaks(hits, 1000, 50, peak.width = 50)
  expect_true(all(res$Peaks$Peak == 15))

})
