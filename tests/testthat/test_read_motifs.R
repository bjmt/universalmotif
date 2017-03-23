library(universalmotif)
context("Test for correct read_motifs output")

test_meme <- system.file("extdata", "minimal.meme",
                         package = "universalmotif")
test_longmeme <- system.file("extdata", "full.meme",
                             package = "universalmotif")
test_homer <- system.file("extdata", "example.homer",
                          package = "universalmotif")
test_jaspar <- system.file("extdata", "example.jaspar",
                           package = "universalmotif")
test_transfac <- system.file("extdata", "transfac.txt",
                             package = "universalmotif")
test_uniprobe <- system.file("extdata", "uniprobe.txt",
                             package = "universalmotif")

test_that("read_motifs correctly autodetects motif format", {

  expect_equal(length(read_motifs(test_meme)), 2)
  expect_identical(class(read_motifs(test_meme)[[1]])[1], "universalmotif")
  expect_equal(length(read_motifs(test_longmeme)), 3)
  expect_identical(class(read_motifs(test_longmeme)[[1]])[1], "universalmotif")
  expect_equal(length(read_motifs(test_homer)), 3)
  expect_identical(class(read_motifs(test_homer)[[1]])[1], "universalmotif")
  expect_equal(length(read_motifs(test_jaspar)), 2)
  expect_identical(class(read_motifs(test_jaspar)[[1]])[1], "universalmotif")
  expect_equal(length(read_motifs(test_transfac)), 5)
  expect_identical(class(read_motifs(test_transfac)[[1]])[1], "universalmotif")
  expect_equal(length(read_motifs(test_uniprobe)), 5)
  expect_identical(class(read_motifs(test_uniprobe)[[1]])[1], "universalmotif")

})

test_that("read_meme correctly rejects wrong formats", {
  
  expect_error(read_meme(test_homer))
  expect_error(read_meme(test_jaspar))
  expect_error(read_meme(test_transfac))
  expect_error(read_meme(test_uniprobe))

})

test_that("read_homer correctly rejects wrong formats", {

  expect_error(read_homer(test_meme))
  expect_error(read_homer(test_longmeme))
  expect_error(read_homer(test_jaspar))
  expect_error(read_homer(test_transfac))
  expect_error(read_homer(test_uniprobe))

})

test_that("read_jaspar correctly rejects wrong formats", {

  expect_error(read_jaspar(test_meme))
  expect_error(read_jaspar(test_longmeme))
  expect_error(read_jaspar(test_homer))
  expect_error(read_jaspar(test_transfac))
  expect_error(read_jaspar(test_uniprobe))

})

test_that("read_transfac correctly rejects wrong formats", {

  expect_error(read_transfac(test_meme))
  expect_error(read_transfac(test_longmeme))
  expect_error(read_transfac(test_homer))
  expect_error(read_transfac(test_jaspar))
  expect_error(read_transfac(test_uniprobe))

})

test_that("read_uniprobe correctly rejects wrong formats", {

  expect_error(read_uniprobe(test_meme))
  expect_error(read_uniprobe(test_longmeme))
  expect_error(read_uniprobe(test_homer))
  expect_error(read_uniprobe(test_jaspar))
  expect_error(read_uniprobe(test_transfac))

})
