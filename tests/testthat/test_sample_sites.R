context("sample_sites()")

test_that("site sampling works", {

  m <- create_motif("AAAAA", nsites = 100)
  s <- sample_sites(m)

  m2 <- create_motif(create_sequences("RNA", seqlen = 10), add.multifreq = 2)
  s2 <- sample_sites(m2, use.freq = 2)

  expect_true(all(as.character(s) == "AAAAA"))

  expect_s4_class(s2, "RNAStringSet")

})
