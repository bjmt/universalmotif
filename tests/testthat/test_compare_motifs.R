context("compare_motifs()")

test_that("basic comparison works", {

  motif1 <- create_motif("TTTWWW")
  motif2 <- create_motif("GGGWWW")

  res.pcc <- compare_motifs(list(motif1, motif2), method = "PCC", score.strat = "sum")
  res.eucl <- compare_motifs(list(motif1, motif2), method = "EUCL", score.strat = "sum")
  res.sw <- compare_motifs(list(motif1, motif2), method = "SW", score.strat = "sum")
  res.kl <- compare_motifs(list(motif1, motif2), method = "KL", score.strat = "sum")
  res.allr <- compare_motifs(list(motif1, motif2), method = "ALLR", score.strat = "sum")

  expect_equal(as.vector(res.pcc),
               c(6, 2, 2, 6))
  expect_equal(round(as.vector(res.eucl), 0),
               c(0, 4, 4, 0))
  expect_equal(round(as.vector(res.sw), 0),
               c(12, 6, 6, 12))
  expect_equal(round(as.vector(res.kl), digits = 1),
               c(0, 13.8, 13.8, 0))
  expect_equal(round(as.vector(res.allr), digits = 3),
               c(5.929, -7.916, -7.916, 5.929))

  res.pcc2 <- compare_motifs(list(motif1, motif2), method = "PCC",
                             use.type = "ICM", score.strat = "sum")
  res.eucl2 <- compare_motifs(list(motif1, motif2), method = "EUCL",
                              use.type = "ICM", score.strat = "sum")
  res.sw2 <- compare_motifs(list(motif1, motif2), method = "SW",
                            use.type = "ICM", score.strat = "sum")
  res.kl2 <- compare_motifs(list(motif1, motif2), method = "KL",
                            use.type = "ICM", score.strat = "sum")

  expect_equal(round(as.vector(res.pcc2), digits = 1),
               c(6, 2, 2, 6))
  expect_equal(round(as.vector(res.eucl2), 0),
               c(0, 8, 8, 0))
  expect_equal(round(as.vector(res.sw2), 0),
               c(12, -9, -9, 12))
  expect_equal(round(as.vector(res.kl2), digits = 1),
               c(0, 27.8, 27.8, 0))

})

test_that("comparisons with p-values works", {

  motif1 <- create_motif("TTTWWW", name = "mot1")
  motif2 <- create_motif("GGGWWW", name = "mot2")

  res <- compare_motifs(list(motif1, motif2), 1:2, max.p = 1,
                        method = "PCC", score.strat = "a.mean")

  expect_equal(round(as.vector(res$Pval)[1], digits = 2), 0.20)

})

test_that("custom db scores are handled correctly", {

  motif1 <- create_motif("TTTWWW", name = "mot1")
  motif2 <- create_motif("GGGWWW", name = "mot2")

  db.scores <- suppressWarnings(make_DBscores(c(motif1, motif2), method = "PCC",
                             rand.tries = 100, widths = 5:7, progress = FALSE,
                             normalise.scores = FALSE, score.strat = "a.mean"))

  res <- compare_motifs(c(motif1, motif2), 1, db.scores, method = "PCC",
                        max.p = 1)

  expect_true(res$Pval[1] < 1)

})
