context("shuffle_motifs()")

test_that("motif shuffling works", {

  m1 <- create_motif("ATATCGGATACGGCTGTC")
  m2 <- shuffle_motifs(m1)

  expect_s4_class(m2, "universalmotif")
  expect_true(m1@consensus != m2@consensus)

})
