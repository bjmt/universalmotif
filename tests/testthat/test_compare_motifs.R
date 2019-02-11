context("compare_motifs()")

test_that("basic comparison works", {

  motif1 <- create_motif("TTTWWW")
  motif2 <- create_motif("GGGWWW")

  res.mpcc <- compare_motifs(list(motif1, motif2), method = "MPCC",
                             progress = FALSE)
  res.pcc <- compare_motifs(list(motif1, motif2), method = "PCC",
                            progress = FALSE)
  res.eucl <- compare_motifs(list(motif1, motif2), method = "EUCL",
                             progress = FALSE)
  res.meucl <- compare_motifs(list(motif1, motif2), method = "MEUCL",
                              progress = FALSE)
  res.sw <- compare_motifs(list(motif1, motif2), method = "SW",
                           progress = FALSE)
  res.msw <- compare_motifs(list(motif1, motif2), method = "MSW",
                            progress = FALSE)
  res.kl <- compare_motifs(list(motif1, motif2), method = "KL",
                           progress = FALSE)
  res.mkl <- compare_motifs(list(motif1, motif2), method = "MKL",
                            progress = FALSE)

  expect_equal(round(as.vector(res.mpcc), digits = 4),
               c(1, 0.3333, 0.3333, 1))
  expect_equal(as.vector(res.pcc),
               c(36, 12, 12, 36))
  expect_equal(as.vector(res.eucl),
               c(0, 3, 3, 0))
  expect_equal(as.vector(res.meucl),
               c(0, 0.5, 0.5, 0))
  expect_equal(as.vector(res.sw),
               c(0, 6, 6, 0))
  expect_equal(as.vector(res.msw),
               c(2, 1, 1, 2))
  expect_equal(round(as.vector(res.kl), digits = 2),
               c(0, 13.73, 13.73, 0))
  expect_equal(round(as.vector(res.mkl), digits = 3),
               c(0, 2.289, 2.289, 0))

  res.mpcc2 <- compare_motifs(list(motif1, motif2), method = "MPCC",
                              progress = FALSE, use.type = "ICM")
  res.pcc2 <- compare_motifs(list(motif1, motif2), method = "PCC",
                             progress = FALSE, use.type = "ICM")
  res.eucl2 <- compare_motifs(list(motif1, motif2), method = "EUCL",
                              progress = FALSE, use.type = "ICM")
  res.meucl2 <- compare_motifs(list(motif1, motif2), method = "MEUCL",
                               progress = FALSE, use.type = "ICM")
  res.sw2 <- compare_motifs(list(motif1, motif2), method = "SW",
                            progress = FALSE, use.type = "ICM")
  res.msw2 <- compare_motifs(list(motif1, motif2), method = "MSW",
                             progress = FALSE, use.type = "ICM")
  res.kl2 <- compare_motifs(list(motif1, motif2), method = "KL",
                            progress = FALSE, use.type = "ICM")
  res.mkl2 <- compare_motifs(list(motif1, motif2), method = "MKL",
                             progress = FALSE, use.type = "ICM")

  expect_equal(round(as.vector(res.mpcc2), digits = 4),
               c(1, 0.3846, 0.3846, 1))
  expect_equal(round(as.vector(res.pcc2), digits = 2),
               c(36, 13.85, 13.85, 36))
  expect_equal(as.vector(res.eucl2),
               c(0, 6, 6, 0))
  expect_equal(as.vector(res.meucl2),
               c(0, 1, 1, 0))
  expect_equal(as.vector(res.sw2),
               c(0, -9, -9, 0))
  expect_equal(as.vector(res.msw2),
               c(2, -1.5, -1.5, 2))
  expect_equal(round(as.vector(res.kl2), digits = 2),
               c(0, 24.23, 24.23, 0))
  expect_equal(round(as.vector(res.mkl2), digits = 3),
               c(0, 4.038, 4.038, 0))

  res.mpcc <- compare_motifs(list(motif1, motif2), method = "MPCC",
                             progress = FALSE, normalise.scores = TRUE)

  expect_equal(round(as.vector(res.mpcc), digits = 4),
               c(1, 0.3333, 0.3333, 1))

})

test_that("comparisons with p-values works", {

  motif1 <- create_motif("TTTWWW", name = "mot1")
  motif2 <- create_motif("GGGWWW", name = "mot2")

  res <- compare_motifs(list(motif1, motif2), 1:2, max.p = 1,
                        method = "MPCC", progress = FALSE)

  expect_equal(round(res$Pval, digits = 4), 0.1601)

})

test_that("custom db scores are handled correctly", {

  motif1 <- create_motif("TTTWWW", name = "mot1")
  motif2 <- create_motif("GGGWWW", name = "mot2")

  db.scores <- make_DBscores(c(motif1, motif2), method = "PCC",
                             rand.tries = 10, progress = FALSE)

  res <- compare_motifs(c(motif1, motif2), 1:2, db.scores, method = "PCC",
                        normalise.scores = TRUE, max.p = 1, progress = FALSE)

  expect_true(res$Pval < 1)

})
