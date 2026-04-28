context("shuffle_motifs()")

test_that("motif shuffling works", {

  m1 <- create_motif("ATATCGGATACGGCTGTC")
  m2 <- shuffle_motifs(m1)

  expect_s4_class(m2, "universalmotif")
  expect_true(m1@consensus != m2@consensus)

})

test_that("shuffle_motifs() preserves pseudocount and type (regression: both were reset to defaults)", {

  m <- create_motif("ATATCGGATACGGCTGTC", type = "PCM", nsites = 20L, pseudocount = 2)
  s <- shuffle_motifs(m)

  expect_equal(s@pseudocount, 2)
  expect_equal(s@type, "PCM")

})

test_that("shuffle_motifs() handles a list containing a length-1 motif (regression: missing drop=FALSE)", {

  m1 <- create_motif("A", nsites = 10, pseudocount = 1)
  m2 <- create_motif("ACGT", nsites = 10, pseudocount = 1)
  expect_error(shuffle_motifs(c(m1, m2)), regexp = NA)

})
