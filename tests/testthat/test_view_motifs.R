context("Test motif vis")

test_that("motif vis works", {

  m <- create_motif()
  r <- view_motifs(m)

  expect_is(r, "gg")
  expect_equal(attributes(r$data)$class, "waiver")
  expect_equal(length(r$data), 0)

})
