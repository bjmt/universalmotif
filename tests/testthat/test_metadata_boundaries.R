context("motif metadata and format boundaries")

metadata_fixture <- function() {
  train <- Biostrings::DNAStringSet(c("ACGTACGT", "ACGTTCGT", "ACGTACGA"))
  m <- suppressMessages(create_motif(train, name = "metadata", add.multifreq = 2))
  m <- add_gap(m, gaploc = c(2, 5), mingap = c(1, 2), maxgap = c(3, 4))
  m@gapinfo@extrapvals <- c(.125, .25)
  m@extrainfo <- c(source = "fixture")
  m
}

test_that("native serialization preserves multiple gaps and higher-order matrices", {
  m <- metadata_fixture()
  for (minimal in c(FALSE, TRUE)) {
    f <- tempfile()
    on.exit(unlink(f), add = TRUE)
    expect_silent(write_motifs(m, f, minimal = minimal))
    restored <- read_motifs(f)
    expect_equal(restored@gapinfo, m@gapinfo)
    expect_identical(names(restored@multifreq), names(m@multifreq))
    expect_equal(unname(restored@multifreq[["2"]]), unname(m@multifreq[["2"]]),
                 tolerance = 1e-4)
    expect_identical(rownames(restored@multifreq[["2"]]), rownames(m@multifreq[["2"]]))
    expect_equal(restored@motif, m@motif, tolerance = 1e-4)
    if (!minimal) {
      expect_equal(restored@nsites, m@nsites)
      expect_identical(restored@extrainfo, m@extrainfo)
    }
  }
  # Verify that the location is not confused with the minimum gap width.
  constructed <- universalmotif_cpp(m@motif, isgapped = TRUE,
    gaploc = c(2, 5), mingap = c(1, 2), maxgap = c(3, 4))
  expect_equal(constructed@gapinfo@gaploc, c(2, 5))
})

test_that("native gap round trips preserve actual scanning results", {
  m <- add_gap(create_motif("ACGT", pseudocount = 1, nsites = 100),
               gaploc = 2, mingap = 1, maxgap = 2)
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  write_motifs(m, f)
  seqs <- Biostrings::DNAStringSet(c("ACAGTA", "ACAAGT"))
  scan <- function(m) suppressMessages(scan_sequences(m, seqs,
    threshold = .8, threshold.type = "logodds", RC = FALSE, verbose = 0))
  before <- scan(m)
  expect_gt(nrow(before), 0)
  expect_equal(before$match, c("AC.GT", "AC..GT"))
  pwm <- convert_type(m, "PWM")@motif
  expected <- sum(trunc(pwm[cbind(seq_len(4), seq_len(4))] * 1000)) / 1000
  expect_equal(before$score, rep(expected, 2))
  expect_equal(scan(read_motifs(f)), before)
})

test_that("gapped PWM expansion preserves scores with multiple and zero-width gaps", {
  pwm <- matrix(c(1, 0, -1, -1, 0, 1, -1, -1,
                  -1, -1, 1, 0, -1, -1, 0, 1), 4,
                dimnames = list(c("A", "C", "G", "T"), NULL))
  m <- add_gap(create_motif(pwm, type = "PWM"), gaploc = c(1, 3),
               mingap = c(0, 1), maxgap = c(1, 1))
  expanded <- ungap_single(m)
  expect_length(expanded, 2)
  for (x in expanded) {
    expect_identical(x@type, "PWM")
    gaps <- get_gaplocs(x)
    expect_true(all(x@motif[, gaps, drop = FALSE] == 0))
    expect_equal(unname(x@motif[, -gaps, drop = FALSE]), unname(pwm))
  }
})

test_that("old native files without gap fields still load", {
  m <- read_motifs(system.file("extdata", "universalmotif.txt",
                              package = "universalmotif"))
  expect_false(m@gapinfo@isgapped)
  expect_length(m@gapinfo@gaploc, 0L)
})

test_that("malformed native gap metadata is rejected", {
  m <- metadata_fixture()
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  write_motifs(m, f)
  lines <- readLines(f)
  # Use the single reader so an invalid field's diagnostic is not masked by
  # the multi-motif reader's intentionally tolerant partial-file handling.
  start <- which(lines == "---") + 1L
  parsed <- yaml::yaml.load(lines[start:length(lines)])
  parsed$extrapvals <- c(2, 3)
  expect_error(read_motifs_single(yaml::as.yaml(parsed)), "extrapvals")
})

test_that("exports warn when gap or higher-order metadata cannot be represented", {
  m <- metadata_fixture()
  writers <- list(write_meme, write_homer, write_jaspar, write_transfac, write_matrix)
  for (writer in writers) {
    f <- tempfile()
    on.exit(unlink(f), add = TRUE)
    expect_warning(writer(m, f), "cannot preserve gap definitions and higher-order matrices")
    expect_true(file.exists(f))
  }
  expect_identical(convert_motifs(m), m)
  expect_warning(convert_motifs(m, "Biostrings-PWM"), "cannot preserve")
})

test_that("TFBSTools conversion retains supported fields and warns for unsupported ones", {
  skip_if_not_installed("TFBSTools")
  m <- create_motif("ACGTACGT", name = "named", nsites = 40, pseudocount = 0,
                    extrainfo = c(source = "fixture"))
  converted <- convert_motifs(m, "TFBSTools-PFMatrix")
  restored <- convert_motifs(converted)
  expect_identical(restored@name, m@name)
  expect_equal(restored@nsites, m@nsites)
  # PFMatrix does not retain universalmotif's pseudocount setting.
  expect_equal(convert_type(restored, "PPM", pseudocount = 0)@motif,
               convert_type(m, "PPM", pseudocount = 0)@motif)
  expect_warning(convert_motifs(metadata_fixture(), "TFBSTools-PFMatrix"),
                 "cannot preserve")
})

test_that("CWM native round trips keep signs and do not invent site counts", {
  # First-column sum is zero, positive integral, or negative integral.
  for (offset in c(0, 1, -1)) {
    mat <- matrix(c(1 + offset, -1, 0, 0, .5, -.25, .125, -.375), 4,
                  dimnames = list(c("A", "C", "G", "T"), NULL))
    m <- create_motif(mat, type = "CWM")
    expect_length(m@nsites, 0L)
    f <- tempfile()
    on.exit(unlink(f), add = TRUE)
    write_motifs(m, f)
    restored <- read_motifs(f)
    expect_identical(restored@type, "CWM")
    expect_equal(restored@motif, m@motif)
    expect_length(restored@nsites, 0L)
    expect_length(create_motif(mat, type = "CWM", nsites = NA_real_)@nsites, 0L)
    expect_equal(create_motif(mat, type = "CWM", nsites = 12)@nsites, 12)
  }
})

test_that("site counts have an explicit missing and finite-integer contract", {
  mat <- matrix(rep(c(.5, .25, .125, .125), 2), 4,
                dimnames = list(c("A", "C", "G", "T"), NULL))
  for (n in list(numeric(), NA_real_))
    expect_length(create_motif(mat, type = "PPM", nsites = n)@nsites, 0L)
  for (n in c(0, -1, .5, 1.5, Inf, -Inf))
    expect_error(create_motif(mat, nsites = n), "nsites")
  expect_equal(create_motif(mat * 8, type = "PCM")@nsites, 8)
})
