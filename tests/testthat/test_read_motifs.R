context("read_*()")

test_that("read functions work ok", {

  homer <- read_homer(system.file("extdata", "homer.txt",
                                  package="universalmotif"))
  cisbp <- read_cisbp(system.file("extdata", "cisbp.txt",
                                  package="universalmotif"))
  jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
                                    package="universalmotif"))
  dreme <- read_meme(system.file("extdata", "dreme.txt",
                                    package="universalmotif"))
  meme_custom_alph <- read_meme(system.file("extdata", "meme_custom_alph.txt",
                                    package="universalmotif"))
  m <- system.file("extdata", "meme_full.txt", package = "universalmotif")
  meme <- read_meme(file = m)
  transfac <- read_transfac(system.file("extdata", "transfac.txt",
                                        package="universalmotif"))
  uniprobe <- read_uniprobe(system.file("extdata", "uniprobe_full.txt",
                                        package="universalmotif"))
  universalmotif <- read_motifs(system.file("extdata",
                                            "universalmotif.txt",
                                            package="universalmotif"))
  hocomoco <- read_matrix(system.file("extdata", "hocomoco.txt",
                                      package="universalmotif"),
                          headers = ">", alphabet = "DNA", positions = "rows")

  expect_equal(length(homer), 5)
  expect_equal(length(cisbp), 2)
  expect_equal(length(jaspar), 5)
  expect_equal(length(dreme), 3)
  expect_equal(length(meme), 3)
  expect_equal(length(transfac), 5)
  expect_equal(length(uniprobe), 3)
  expect_s4_class(universalmotif, "universalmotif")
  expect_s4_class(hocomoco, "universalmotif")

  meme2 <- read_meme(m, readsites = TRUE)

  expect_true(is.list(meme2))
  expect_equal(names(meme2), c("motifs", "sites"))
  expect_equal(length(meme2$sites), 3)
  expect_s4_class(meme2$sites[[1]], "DNAStringSet")

})

test_that("read_jaspar() derives nsites from PCM column sums (regression: was always numeric(0))", {

  jaspar <- read_jaspar(system.file("extdata", "jaspar.txt", package = "universalmotif"))
  if (!is.list(jaspar)) jaspar <- list(jaspar)
  expect_true(all(vapply(jaspar, function(m) length(m@nsites) == 1L, logical(1))))
  expect_equal(jaspar[[1]]@nsites, 100)

})

test_that("read_motifs() version dispatch uses numeric_version not lexicographic comparison (regression)", {

  # "1.1.100" is > "1.1.67" numerically but < "1.1.67" lexicographically
  # (because "1" < "6" at the third segment's first char).
  # Without numeric_version(), the pre-1.2.x reader path would be wrongly
  # selected for any version whose third segment exceeds one digit.
  expect_true(numeric_version("1.1.100") > numeric_version("1.1.67"))
  expect_true("1.1.100" < "1.1.67")

})

## Parser error-path tests: malformed files should fail with a clear message
## rather than silently parsing to an empty list or propagating -Inf / NA.

write_tmp <- function(text) {
  f <- tempfile()
  writeLines(text, f)
  f
}

test_that("read_jaspar errors on a file with no '>' headers", {
  f <- write_tmp(c("A 1 2 3", "C 1 2 3"))
  on.exit(file.remove(f), add = TRUE)
  expect_error(read_jaspar(f), regexp = "no motifs found")
})

test_that("read_jaspar errors on a header followed by an empty matrix", {
  f <- write_tmp(c(">empty",
                   "A [  ]",
                   "C [  ]",
                   "G [  ]",
                   "T [  ]"))
  on.exit(file.remove(f), add = TRUE)
  expect_error(read_jaspar(f), regexp = "empty motif")
})

test_that("read_meme errors on a file with no MOTIF blocks", {
  f <- write_tmp(c("MEME version 4", "ALPHABET= ACGT",
                   "strands: + -",
                   "Background letter frequencies",
                   "A 0.25 C 0.25 G 0.25 T 0.25"))
  on.exit(file.remove(f), add = TRUE)
  expect_error(read_meme(f), regexp = "no motifs found")
})

test_that("read_transfac errors on a file with no '//' separators", {
  f <- write_tmp(c("ID test", "P0 A C G T", "01 1 2 3 4"))
  on.exit(file.remove(f), add = TRUE)
  expect_error(read_transfac(f), regexp = "no motifs found")
})

test_that("read_cisbp errors on a file with no 'Pos' header", {
  f <- write_tmp(c("not a motif", "no header here"))
  on.exit(file.remove(f), add = TRUE)
  expect_error(read_cisbp(f), regexp = "no motifs found")
})

test_that("read_homer leaves nsites empty when the count fields are missing", {
  ## Regression coverage for R1: ifelse(is.na, numeric(0), x) used to return
  ## NA_real_ rather than an empty slot. With the fix, a malformed/missing
  ## count field stays as numeric(0) and the validity check passes.
  homer_text <- c(
    ">ACGTAC\tnamed-motif\t8.0\tNA\tNA,NA,NA",
    "0.97\t0.01\t0.01\t0.01",
    "0.01\t0.97\t0.01\t0.01",
    "0.01\t0.01\t0.97\t0.01",
    "0.01\t0.01\t0.01\t0.97",
    "0.97\t0.01\t0.01\t0.01",
    "0.01\t0.97\t0.01\t0.01"
  )
  f <- write_tmp(homer_text)
  on.exit(file.remove(f), add = TRUE)
  m <- read_homer(f)
  expect_s4_class(m, "universalmotif")
  expect_equal(length(m@nsites), 0L)
  expect_equal(length(m@bkgsites), 0L)
})

test_that("read_transfac skips header blocks, keeps whole fields and follows P0", {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(c(
    "VV  TRANSFAC MATRIX TABLE, Release 2024", "XX", "//",
    "ID  my_id", "XX", "AC  M00001", "XX", "NA  Some factor",
    "HC  Basic leucine zipper", "OS  Homo sapiens",
    "P0      T      G      C      A",
    "01      1      2      3      4",
    "02      0      0      0     10", "XX", "//"
  ), f)
  m <- read_transfac(f)
  expect_equal(m@name, "M00001")
  expect_equal(m@altname, "my_id")
  expect_equal(m@family, "Basic leucine zipper")
  expect_equal(m@organism, "Homo sapiens")
  expect_equal(unname(m@motif[, 1]), c(4, 3, 2, 1))
  expect_equal(unname(m@motif[, 2]), c(10, 0, 0, 0))
})

test_that("read_cisbp keeps whole field values", {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(c("TF\tT1", "TF Name\tSome name", "Family\tC2H2 ZF",
               "Species\tHomo sapiens", "Pos\tA\tC\tG\tT", "1\t0.1\t0.2\t0.3\t0.4"), f)
  m <- read_cisbp(f)
  expect_equal(c(m@name, m@altname, m@family, m@organism),
               c("T1", "Some name", "C2H2 ZF", "Homo sapiens"))
})

test_that("read_jaspar follows row labels and reads other alphabets", {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(c(">m1 one", "T [ 1 0 ]", "G [ 2 0 ]", "C [ 3 0 ]", "A [ 4 10 ]",
               ">m2 two", "X [ 1 2 ]", "Y [ 3 2 ]"), f)
  m <- read_jaspar(f)
  expect_equal(unname(m[[1]]@motif[, 1]), c(4, 3, 2, 1))
  expect_equal(m[[2]]@alphabet, "XY")
  writeLines(c(">m3", paste(rep(5, 20), collapse = " ")), f)
  expect_error(read_jaspar(f))
  writeLines(c(">m3", rep("[ 5 ]", 20)), f)
  expect_equal(read_jaspar(f)@alphabet, "AA")
})

test_that("read_homer keeps spaces in tab-separated names", {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(c(">ACGT\tmy motif name\t5.0\t-10\t0\tT:10(1%),B:5(0.5%),P:1e-5",
               "0.7\t0.1\t0.1\t0.1", "0.1\t0.7\t0.1\t0.1"), f)
  expect_equal(read_homer(f)@name, "my motif name")
})

test_that("read_meme matches custom backgrounds by letter", {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(c("MEME version 5", "", "ALPHABET \"Test\"", "Q", "R", "1", "2",
               "END ALPHABET", "", "strands: +", "",
               "Background letter frequencies", "Q 0.1 R 0.2 1 0.3 2 0.4", "",
               "MOTIF m1", "letter-probability matrix: alength= 4 w= 1 nsites= 10",
               " 0.7  0.1  0.1  0.1"), f)
  m <- read_meme(f)
  expect_equal(m@bkg[c("1", "2", "Q", "R")], c(`1` = .3, `2` = .4, Q = .1, R = .2))
  expect_equal(unname(m@motif["Q", 1]), .7)
})

test_that("read_meme reads sites at large positions and long block diagrams", {
  lines <- readLines(system.file("extdata", "meme_full.txt", package = "universalmotif"))
  i <- grep("^AT1G21910                \\(  744\\)", lines)[1]
  lines[i] <- sub("(  744)", "(12744)", lines[i], fixed = TRUE)
  j <- grep("^    \\[-3\\(2.60e-05\\)\\]_308_", lines)
  lines[j] <- "    [-3(2.60e-05)]_308_[+3(7.17e-05)]_129_\\"
  lines <- append(lines, "    [-2(7.96e-07)]_74_[+2(1.21e-05)]_56", after = j)
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(lines, f)
  m <- read_meme(f, readsites = TRUE, readsites.meta = TRUE)
  expect_equal(unname(as.character(m$sites[[1]][1])), "GAAAAAAAAAGAAAC")
  d <- m$sites.meta.combined
  expect_equal(d$Diagram[d$Sequence == "AT1G80840"],
               "186_[-3(2.60e-05)]_156_[-3(2.60e-05)]_308_[+3(7.17e-05)]_129_[-2(7.96e-07)]_74_[+2(1.21e-05)]_56")
})

test_that("rounded probabilities are read as a PPM", {
  f <- tempfile()
  on.exit(unlink(f))
  writeLines(c("0.25 0.25 0.25 0.23", "0.5 0.5 0 0.02"), f)
  m <- read_matrix(f, positions = "rows", headers = FALSE)
  expect_equal(m@type, "PPM")
  expect_equal(unname(colSums(m@motif)), c(1, 1))
})
