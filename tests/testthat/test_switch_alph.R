context("switch_alph()")

test_that("motif alphabet switching works", {

  m1 <- create_motif()
  m2 <- switch_alph(m1)

  expect_s4_class(m2, "universalmotif")
  expect_equal(m2@alphabet, "RNA")

})

test_that("DNA -> RNA -> DNA round-trips", {
  m <- create_motif("ACGT", alphabet = "DNA")
  back <- switch_alph(switch_alph(m))
  expect_equal(back@alphabet, "DNA")
  expect_equal(back@consensus, m@consensus)
})

test_that("AA motif cannot be switched", {
  m <- create_motif("AAYY", alphabet = "AA")
  expect_error(switch_alph(m), regexp = "will not convert|AA")
})

test_that("strand and name slots survive the switch", {
  m <- create_motif("ACGT", name = "named", alphabet = "DNA")
  m@strand <- "+"
  m2 <- switch_alph(m)
  expect_equal(m2@name, "named")
  expect_equal(m2@strand, "+")
})

test_that("row names follow the alphabet flip", {
  m <- create_motif("ACGT", alphabet = "DNA")
  m2 <- switch_alph(m)
  expect_equal(rownames(m2@motif), c("A", "C", "G", "U"))
  m3 <- switch_alph(m2)
  expect_equal(rownames(m3@motif), c("A", "C", "G", "T"))
})
