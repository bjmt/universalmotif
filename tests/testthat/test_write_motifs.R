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
