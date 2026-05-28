context("write_*()")

test_that("write functions are ok", {

  m <- create_motif("GTCAWTGAGA", name="test-motif", nsites=100)

  f.homer <- tempfile()
  f.jaspar <- tempfile()
  f.meme <- tempfile()
  f.transfac <- tempfile()
  f.universalmotif <- tempfile()
  f.matrix <- tempfile()

  write_homer(m, f.homer)
  write_jaspar(m, f.jaspar)
  write_meme(m, f.meme)
  write_transfac(m, f.transfac)
  write_motifs(m, f.universalmotif)
  write_matrix(m, f.matrix)

  homer <- read_homer(f.homer)
  jaspar <- read_jaspar(f.jaspar)
  meme <- read_meme(f.meme)
  transfac <- read_transfac(f.transfac)
  universalmotif <- read_motifs(f.universalmotif)
  matrix <- read_matrix(f.matrix)

  expect_true(homer@icscore >= 18 && homer@icscore <= 19)
  expect_true(jaspar@icscore >= 18 && jaspar@icscore <= 19)
  expect_true(meme@icscore >= 18 && meme@icscore <= 19)
  expect_true(transfac@icscore >= 18 && transfac@icscore <= 19)
  expect_true(universalmotif@icscore >= 18 && universalmotif@icscore <= 19)
  expect_true(matrix@icscore >= 18 && matrix@icscore <= 19)

  file.remove(f.homer)
  file.remove(f.jaspar)
  file.remove(f.meme)
  file.remove(f.transfac)
  file.remove(f.universalmotif)
  file.remove(f.matrix)

})

## Round-trip identity tests: read a canonical fixture, write to tempfile, read
## back, and confirm the round-tripped subset matches. Lossy formats only
## preserve part of the motif state; the comparisons below pick the slots each
## format is expected to round-trip cleanly.

extract_matrix <- function(m) {
  ## Comparable matrix representation (PPM, padded to width).
  m2 <- convert_type(m, "PPM")
  unname(m2@motif)
}

roundtrip_matrix <- function(m, writer, reader) {
  f <- tempfile()
  on.exit(file.remove(f), add = TRUE)
  writer(m, f)
  m2 <- reader(f)
  if (is.list(m2)) m2 <- m2[[1]]
  m2
}

test_that("MEME read/write round-trip preserves matrix and name", {
  m <- read_meme(system.file("extdata", "meme_minimal.txt",
                             package = "universalmotif"))[[1]]
  m2 <- roundtrip_matrix(m, write_meme, read_meme)
  expect_equal(extract_matrix(m), extract_matrix(m2), tolerance = 1e-4)
  expect_equal(m@name, m2@name)
  expect_equal(m@alphabet, m2@alphabet)
})

test_that("JASPAR read/write round-trip preserves matrix and name", {
  m <- read_jaspar(system.file("extdata", "jaspar.txt",
                                package = "universalmotif"))[[1]]
  m2 <- roundtrip_matrix(m, write_jaspar, read_jaspar)
  expect_equal(extract_matrix(m), extract_matrix(m2), tolerance = 1e-4)
  expect_equal(m@name, m2@name)
})

test_that("HOMER read/write round-trip preserves matrix and name", {
  m <- read_homer(system.file("extdata", "homer.txt",
                               package = "universalmotif"))[[1]]
  m2 <- roundtrip_matrix(m, write_homer, read_homer)
  expect_equal(extract_matrix(m), extract_matrix(m2), tolerance = 1e-3)
  expect_equal(m@name, m2@name)
})

test_that("TRANSFAC read/write round-trip preserves matrix and name", {
  m <- read_transfac(system.file("extdata", "transfac.txt",
                                  package = "universalmotif"))[[1]]
  m2 <- roundtrip_matrix(m, write_transfac, read_transfac)
  expect_equal(extract_matrix(m), extract_matrix(m2), tolerance = 1e-4)
  expect_equal(m@name, m2@name)
})

test_that("universalmotif read/write round-trip preserves the full object", {
  m <- create_motif("GTCAWTGAGA", name = "rt-test", nsites = 50)
  f <- tempfile(); on.exit(file.remove(f), add = TRUE)
  write_motifs(m, f)
  m2 <- read_motifs(f)
  expect_equal(extract_matrix(m), extract_matrix(m2), tolerance = 1e-8)
  expect_equal(m@name, m2@name)
  expect_equal(m@nsites, m2@nsites)
  expect_equal(m@alphabet, m2@alphabet)
})

test_that("matrix read/write round-trip preserves the matrix", {
  m <- create_motif("GTCAWTGAGA", name = "mat-test", nsites = 50)
  f <- tempfile(); on.exit(file.remove(f), add = TRUE)
  write_matrix(m, f)
  m2 <- read_matrix(f)
  expect_equal(extract_matrix(m), extract_matrix(m2), tolerance = 1e-4)
})
