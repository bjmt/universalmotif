context("create_sequences()")

test_that("sequence creation works", {

  expect_s4_class(create_sequences(), "DNAStringSet")
  expect_s4_class(create_sequences(alphabet = "RNA"), "RNAStringSet")
  expect_s4_class(create_sequences(alphabet = "AA"), "AAStringSet")
  expect_s4_class(create_sequences(alphabet = "QWER"), "BStringSet")

  d <- c(AA = 0.0625, AC = 0.0625, AG = 0.0625, AT = 0.0625,
         CA = 0.0625, CC = 0.0625, CG = 0.0625, CT = 0.0625,
         GA = 0.0625, GC = 0.0625, GG = 0.0625, GT = 0.0625,
         TA = 0.0625, TC = 0.0625, TG = 0.0625, TT = 0.0625)

  r <- d
  names(r) <- gsub("T", "U", names(r))

  s1 <- create_sequences(difreqs = d)
  s2 <- create_sequences(alphabet = "RNA", difreqs = r)

  expect_s4_class(s1, "DNAStringSet")
  expect_s4_class(s2, "RNAStringSet")

  tri <- rep(1 / 64, 64)
  names(tri) <- universalmotif:::DNA_TRI

  tri.r <- tri
  names(tri.r) <- gsub("T", "U", names(tri.r))

  s3 <- create_sequences(trifreqs = tri)
  s4 <- create_sequences("RNA", trifreqs = tri.r)

  expect_s4_class(s3, "DNAStringSet")
  expect_s4_class(s4, "RNAStringSet")

})
