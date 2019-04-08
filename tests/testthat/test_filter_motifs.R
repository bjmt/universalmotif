context("filter_motifs()")

test_that("motif filtering works", {

  m1 <- create_motif("AAAAAAA", name = "motif1", altname = "asdf",
                     family = "zxcv", organism = "rtyu", strand = "-",
                     pval = 0.123, qval = 0.3456, eval = 0.6789)
  m2 <- create_motif("NNNNN", name = "motif2", altname = "qwer",
                     family = "sdfg", organism = "dfhg", strand = "+",
                     pval = 0.3456, qval = 0.123, eval = 0.345)
  m <- list(m1, m2)

  f <- filter_motifs(m, altname = "asdf", family = "zxcv", organism = "rtyu",
                     strand = "-", pval = 0.2)
  f2 <- filter_motifs(m, qval = 0.2, eval = 0.4)

  expect_equal(length(filter_motifs(m, name = "motif1")), 1)
  expect_equal(length(filter_motifs(m, width = 6)), 1)
  expect_equal(length(filter_motifs(m, icscore = 12)), 1)

  expect_equal(f[[1]], m1)
  expect_equal(f2[[1]], m2)

})
