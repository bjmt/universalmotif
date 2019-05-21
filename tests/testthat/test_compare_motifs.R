context("compare_motifs()")

test_that("basic comparison works", {

  motif1 <- create_motif("TTTWWW")
  motif2 <- create_motif("GGGWWW")

  res.mpcc <- compare_motifs(list(motif1, motif2), method = "MPCC")
  res.pcc <- compare_motifs(list(motif1, motif2), method = "PCC")
  res.eucl <- compare_motifs(list(motif1, motif2), method = "EUCL")
  res.meucl <- compare_motifs(list(motif1, motif2), method = "MEUCL")
  res.sw <- compare_motifs(list(motif1, motif2), method = "SW")
  res.msw <- compare_motifs(list(motif1, motif2), method = "MSW")
  res.kl <- compare_motifs(list(motif1, motif2), method = "KL")
  res.mkl <- compare_motifs(list(motif1, motif2), method = "MKL")
  res.allr <- compare_motifs(list(motif1, motif2), method = "ALLR")
  res.mallr <- compare_motifs(list(motif1, motif2), method = "MALLR")

  expect_equal(round(as.vector(res.mpcc), digits = 4),
               c(1, 0.3333, 0.3333, 1))
  expect_equal(as.vector(res.pcc),
               c(6, 2, 2, 6))
  expect_equal(round(as.vector(res.eucl), 0),
               c(0, 4, 4, 0))
  expect_equal(round(as.vector(res.meucl), 1),
               c(0, 0.7, 0.7, 0))
  expect_equal(round(as.vector(res.sw), 0),
               c(12, 6, 6, 12))
  expect_equal(round(as.vector(res.msw), 0),
               c(2, 1, 1, 2))
  expect_equal(round(as.vector(res.kl), digits = 1),
               c(0, 13.8, 13.8, 0))
  expect_equal(round(as.vector(res.mkl), digits = 2),
               c(0, 2.31, 2.31, 0))
  expect_equal(round(as.vector(res.allr), digits = 3),
               c(5.929, -7.916, -7.916, 5.929))
  expect_equal(round(as.vector(res.mallr), digits = 3),
               c(0.988, -1.319, -1.319, 0.988))

  res.mpcc2 <- compare_motifs(list(motif1, motif2), method = "MPCC",
                              use.type = "ICM")
  res.pcc2 <- compare_motifs(list(motif1, motif2), method = "PCC",
                             use.type = "ICM")
  res.eucl2 <- compare_motifs(list(motif1, motif2), method = "EUCL",
                              use.type = "ICM")
  res.meucl2 <- compare_motifs(list(motif1, motif2), method = "MEUCL",
                               use.type = "ICM")
  res.sw2 <- compare_motifs(list(motif1, motif2), method = "SW",
                            use.type = "ICM")
  res.msw2 <- compare_motifs(list(motif1, motif2), method = "MSW",
                             use.type = "ICM")
  res.kl2 <- compare_motifs(list(motif1, motif2), method = "KL",
                            use.type = "ICM")
  res.mkl2 <- compare_motifs(list(motif1, motif2), method = "MKL",
                             use.type = "ICM")

  expect_equal(round(as.vector(res.mpcc2), digits = 3),
               c(1, 0.333, 0.333, 1))
  expect_equal(round(as.vector(res.pcc2), digits = 1),
               c(6, 2, 2, 6))
  expect_equal(round(as.vector(res.eucl2), 0),
               c(0, 8, 8, 0))
  expect_equal(round(as.vector(res.meucl2), 0),
               c(0, 1, 1, 0))
  expect_equal(round(as.vector(res.sw2), 0),
               c(12, -9, -9, 12))
  expect_equal(round(as.vector(res.msw2), 1),
               c(2, -1.5, -1.5, 2))
  expect_equal(round(as.vector(res.kl2), digits = 1),
               c(0, 27.8, 27.8, 0))
  expect_equal(round(as.vector(res.mkl2), digits = 2),
               c(0, 4.64, 4.64, 0))

  res.mpcc <- compare_motifs(list(motif1, motif2), method = "MPCC",
                             normalise.scores = TRUE)

  expect_equal(round(as.vector(res.mpcc), digits = 3),
               c(1, 0.333, 0.333, 1))

})

test_that("comparisons with p-values works", {

  motif1 <- create_motif("TTTWWW", name = "mot1")
  motif2 <- create_motif("GGGWWW", name = "mot2")

  res <- compare_motifs(list(motif1, motif2), 1:2, max.p = 1,
                        method = "MPCC")

  expect_equal(round(res$Pval, digits = 2), 0.19)

})

test_that("custom db scores are handled correctly", {

  motif1 <- create_motif("TTTWWW", name = "mot1")
  motif2 <- create_motif("GGGWWW", name = "mot2")

  db.scores <- suppressWarnings(make_DBscores(c(motif1, motif2), method = "PCC",
                             rand.tries = 50, widths = 5:7, progress = FALSE))

  res <- compare_motifs(c(motif1, motif2), 1, db.scores, method = "PCC",
                        max.p = 1)

  expect_true(res$Pval < 1)

})
