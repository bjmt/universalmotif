context("create_sequences()")

test_that("sequence creation works", {

  expect_s4_class(create_sequences(), "DNAStringSet")
  expect_s4_class(create_sequences(alphabet = "RNA"), "RNAStringSet")
  expect_s4_class(create_sequences(alphabet = "AA"), "AAStringSet")
  expect_s4_class(create_sequences(alphabet = "QWER"), "BStringSet")

  d <- rep(1 / 16, 16)
  names(d) <- universalmotif:::DNA_DI

  r <- d
  names(r) <- gsub("T", "U", names(r))

  s1 <- create_sequences(freqs = d)
  s2 <- create_sequences(alphabet = "RNA", freqs = r)

  expect_s4_class(s1, "DNAStringSet")
  expect_s4_class(s2, "RNAStringSet")

  tri <- rep(1 / 64, 64)
  names(tri) <- get_klets(DNA_BASES, 3)

  tri.r <- tri
  names(tri.r) <- gsub("T", "U", names(tri.r))

  s3 <- create_sequences(freqs = tri)
  s4 <- create_sequences("RNA", freqs = tri.r)

  expect_s4_class(s3, "DNAStringSet")
  expect_s4_class(s4, "RNAStringSet")

})
