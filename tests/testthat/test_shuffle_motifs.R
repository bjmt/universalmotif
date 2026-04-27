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
