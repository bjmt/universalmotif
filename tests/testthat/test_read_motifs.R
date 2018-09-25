context("Test read functions")
library(universalmotif)

test_that("read functions work ok", {

  homer <- read_homer(system.file("extdata", "homer.txt",
                                  package="universalmotif"))
  cisbp <- read_cisbp(system.file("extdata", "cisbp.txt",
                                  package="universalmotif"))
  jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
                                    package="universalmotif"))
  meme <- read_meme(system.file("extdata", "meme_full.txt",
                                package="universalmotif"))
  transfac <- read_transfac(system.file("extdata", "transfac.txt",
                                        package="universalmotif"))
  uniprobe <- read_uniprobe(system.file("extdata", "uniprobe_full.txt", 
                                        package="universalmotif"))
  universalmotif <- read_motifs(system.file("extdata", "universalmotif.txt",
                                            package="universalmotif"))

  expect_equal(length(homer), 5)
  expect_equal(length(cisbp), 2)
  expect_equal(length(jaspar), 5)
  expect_equal(length(meme), 3)
  expect_equal(length(transfac), 5)
  expect_equal(length(uniprobe), 3)
  expect_s4_class(universalmotif, "universalmotif")

})
