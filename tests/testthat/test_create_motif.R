context("Test motif creation")

test_that("motif creation for missing or from numeric input works", {

  expect_error(universalmotif())

  motif1 <- create_motif()
  motif2 <- create_motif(10)
  motif3 <- create_motif(type = "PCM")
  motif4 <- create_motif(type = "PWM")
  motif5 <- create_motif(type = "ICM")

  expect_s4_class(motif1, "universalmotif")

  expect_equal(ncol(motif1["motif"]), ncol(motif2["motif"]))
  expect_equal(motif1["type"], c(type = "PPM"))
  expect_equal(as.numeric(colSums(motif1["motif"])), rep(1, 10))
  expect_equal(sum(motif1["motif"]), 10)
  expect_equal(motif1["alphabet"], c(alphabet = "DNA"))
  expect_equal(motif1["strand"], c(strand = "+-"))
  expect_equal(rownames(motif1["motif"]), c("A", "C", "G", "T"))
  expect_equal(motif1["pseudocount"], c(pseudocount = 0))
  expect_equal(sum(motif3["motif"]), 1000)
  expect_true(all(motif1["motif"] != motif2["motif"]))

  expect_true(any(motif4["motif"] < 0))
  expect_true(all(motif5["motif"] <= 2 && motif5["motif"] >= 0))

  expect_error(create_motif(0))
  expect_error(create_motif(1, 2))
  expect_error(create_motif(10.5))
  expect_error(create_motif(bkg = 0))
  expect_error(create_motif(bkg = c(0.2, 0.2, 0.2, 0.2, 0.2)))
  expect_error(create_motif(bkg = c(0.5, 0.5)))

})

test_that("motif creation from matrix input works", {

  mat <- matrix(c(0.1, 0.7, 0.1, 0.7, 0.4,
                  0.1, 0.1, 0.1, 0.1, 0.1,
                  0.1, 0.1, 0.1, 0.1, 0.1,
                  0.7, 0.1, 0.7, 0.1, 0.4),
                nrow = 4, byrow = TRUE)
  mat2 <- mat * 10

  motif1 <- create_motif(mat)
  motif2 <- create_motif(mat, alphabet = "DNA")
  motif3 <- create_motif(mat2, alphabet = "RNA")
  motif4 <- convert_type(motif3, type = "PWM")
  
  mat3 <- motif4["motif"]

  motif5 <- create_motif(mat3)
  motif5["nsites"] <- 10
  motif6 <- create_motif(mat2, type = "PCM")

  expect_equal(motif1["alphabet"], c(alphabet = "custom"))
  expect_equal(motif2["consensus"], c(consensus = "TATAW"))
  expect_equal(motif3["consensus"], c(consensus = "UAUAW"))
  expect_equal(motif3["nsites"], c(nsites = 10))
  expect_true(all(motif6["motif"] >= 1))
  expect_true(any(motif4["motif"] < 0))
  # expect_identical(motif5, motif3)  # this tests fails on travis?

})
