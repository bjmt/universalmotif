context("universalmotif-methods")

test_that("accessor functions work", {

  m <- create_motif("SGDGNTGGAY", pseudocount = 1, nsites = 88,
                    family = "asdf", organism = "qwer")
  
  expect_equal(unname(m["name"]), m@name)
  expect_equal(m["altname"], m@altname)
  expect_equal(unname(m["family"]), m@family)
  expect_equal(unname(m["organism"]), m@organism)
  expect_equal(m["motif"], m@motif)
  expect_equal(unname(m["alphabet"]), m@alphabet)
  expect_equal(unname(m["type"]), m@type)
  expect_equal(unname(m["icscore"]), m@icscore)
  expect_equal(unname(m["nsites"]), m@nsites)
  expect_equal(unname(m["pseudocount"]), m@pseudocount)
  expect_equal(m["bkg"], m@bkg)
  expect_equal(m["bkgsites"], m@bkgsites)
  expect_equal(unname(m["consensus"]), m@consensus)
  expect_equal(unname(m["strand"]), m@strand)
  expect_equal(m["pval"], m@pval)
  expect_equal(m["qval"], m@qval)
  expect_equal(m["eval"], m@eval)
  expect_equal(m["multifreq"], m@multifreq)
  expect_equal(unname(m["extrainfo"]), m@extrainfo)

  m["altname"] <- "zxcv"
  expect_equal(m@altname, "zxcv")
  expect_error(m["motif"] <- matrix(1:8, nrow = 2))

})

test_that("misc methods work", {

  m <- create_motif("SGDGNTGGAY", nsites=12)

  expect_equal(as.data.frame(m),
               data.frame(name="motif", altname=as.character(NA),
                          family=as.character(NA), organism=as.character(NA),
                          alphabet="DNA", icscore=m@icscore, nsites=12,
                          bkgsites=as.numeric(NA), consensus="SGDGNTGGAY",
                          strand="+-", pval=as.numeric(NA),
                          qval=as.numeric(NA), eval=as.numeric(NA),
                          stringsAsFactors=FALSE))
  expect_equal(subset(m, 1:4), create_motif("SGDG", nsites=12))
  expect_equal(round(rowMeans(m), 4), c(A=0.1583, C=0.1250, G=0.5083, T=0.2083))
  expect_equal(colMeans(m), c(S=0.25, G=0.25, D=0.25, G=0.25, N=0.25, T=0.25,
                              G=0.25, G=0.25, A=0.25, Y=0.25))
  expect_equal(colSums(m), c(S=1, G=1, D=1, G=1, N=1, T=1, G=1, G=1, A=1, Y=1))
  # expect_equal(round(rowSums(m), 4), c(A=1.583, C=1.250, G=5.083, T=2.083))
  expect_equal(nrow(m), 4)
  expect_equal(ncol(m), 10)
  expect_equal(rownames(m), c("A", "C", "G", "T"))
  expect_equal(colnames(m), c("S", "G", "D", "G", "N", "T", "G", "G", "A", "Y"))
  expect_equal(ncol(cbind(m, m)), 20)

})
