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
