context("Test motif reverse complement")

test_that("motif reverse complement works", {

  m1 <- create_motif("AAAAA")
  m2 <- motif_rc(m1)

  expect_equal(m2@consensus, "TTTTT")

})
