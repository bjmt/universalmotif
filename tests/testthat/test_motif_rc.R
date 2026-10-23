context("motif_rc()")

test_that("motif reverse complement works", {

  m1 <- create_motif("AAAAA")
  m2 <- motif_rc(m1)

  expect_equal(m2@consensus, "TTTTT")

})

test_that("non-palindromic motif has a different RC", {
  m <- create_motif("AAAATTTC")
  rc <- motif_rc(m)
  expect_false(identical(m@consensus, rc@consensus))
  expect_equal(rc@consensus, "GAAATTTT")
})

test_that("motif_rc is self-inverse (round-trip identity on matrix)", {
  m <- create_motif("AAAATTTC")
  expect_equal(m@motif, motif_rc(motif_rc(m))@motif)
})

test_that("RNA motif_rc swaps A<->U and C<->G", {
  m <- create_motif("ACGU", alphabet = "RNA")
  rc <- motif_rc(m)
  ## RC of ACGU is ACGU (palindrome on the RNA alphabet by AU/CG pairing).
  expect_equal(rc@consensus, "ACGU")
})

test_that("motif_rc errors or no-ops on AA motifs", {
  m <- create_motif("AAYY", alphabet = "AA")
  expect_error(motif_rc(m))
})

test_that("motif_rc keeps gaps and higher-order matrices", {
  g <- add_gap(create_motif("AAAACCCCGG"), gaploc = 3, mingap = 1, maxgap = 3)
  rc <- motif_rc(g)
  expect_equal(rc@gapinfo@gaploc, 7)
  expect_equal(c(rc@gapinfo@mingap, rc@gapinfo@maxgap), c(1, 3))
  expect_true(rc@gapinfo@isgapped)
  seqs <- Biostrings::DNAStringSet(c("AACGT", "AACGT", "ATCGA"))
  m <- add_multifreq(create_motif(seqs), seqs, add.k = 2)
  rcs <- Biostrings::reverseComplement(seqs)
  m.rc <- add_multifreq(create_motif(rcs), rcs, add.k = 2)
  expect_equal(unname(motif_rc(m)@multifreq[["2"]]), unname(m.rc@multifreq[["2"]]))
  expect_equal(motif_rc(motif_rc(m))@multifreq, m@multifreq)
})
