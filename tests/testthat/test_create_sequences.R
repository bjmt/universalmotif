context("Test sequence creation")
library(universalmotif)

test_that("sequence creation works", {

  expect_s4_class(create_sequences(), "DNAStringSet")
  expect_s4_class(create_sequences(alphabet = "RNA"), "RNAStringSet")
  expect_s4_class(create_sequences(alphabet = "AA"), "AAStringSet")
  expect_s4_class(create_sequences(alphabet = "QWER"), "BStringSet")

})
