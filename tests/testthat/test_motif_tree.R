context("motif_tree()")

test_that("motif trees work", {

  m1 <- create_motif("GCGCGCGCWT", nsites = 100, name = "motif-1")
  m2 <- create_motif("AWAWAWAWAW", nsites = 100, name = "motif-2")
  m3 <- create_motif("NTGNTGTNGN", nsites = 100, name = "motif-3")

  m <- list(m1, m2, m3)
  r <- motif_tree(m, progress = FALSE)

  expect_is(r, "ggtree")
  expect_true(any(r$data$branch.length > 0))
  expect_true(length(r$data$branch.length) > 0)

})
