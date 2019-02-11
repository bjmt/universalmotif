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
  expect_equal(as.vector(res.pcc), c(36, 12, 12, 36))
  expect_equal(as.vector(res.eucl), c(0, 3, 3, 0))
  expect_equal(as.vector(res.meucl), c(0, 0.5, 0.5, 0))
  expect_equal(as.vector(res.sw), c(0, 6, 6, 0))
  expect_equal(as.vector(res.msw), c(2, 1, 1, 2))
  expect_equal(round(as.vector(res.kl), digits = 2),
               c(0, 13.73, 13.73, 0))
  expect_equal(round(as.vector(res.mkl), digits = 3),
               c(0, 2.289, 2.289, 0))

})

test_that("comparisons with p-values works", {

  motif1 <- create_motif("TTTWWW", name = "mot1")
  motif2 <- create_motif("GGGWWW", name = "mot2")

  res <- compare_motifs(list(motif1, motif2), 1:2, max.p = 1,
                        method = "MPCC", progress = FALSE)

  expect_equal(round(res$Pval, digits = 4), 0.1601)

})
