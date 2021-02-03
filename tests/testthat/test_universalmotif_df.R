context("universalmotif_df")

test_that("to_df(), update_motifs(), to_list() work", {
  m1 <- create_motif()
  m2 <- create_motif()
  m3 <- create_motif()
  m4 <- create_motif()
  mydf <- to_df(c(m1, m2, m3, m4))
  expect_equal(mydf$name[1], "motif")
  mydf$name <- LETTERS[1:4]
  expect_equal(mydf$name[1], "A")
  mydf <- update_motifs(mydf)
  expect_equal(mydf$name[1], "A")
  m <- to_list(mydf)
  expect_equal(m[[1]]["name"], "A")
})

test_that("extrainfo gets moved around correctly", {
  m1 <- create_motif()
  m2 <- create_motif()
  m1["extrainfo"] <- c(A = "B", C = "D", E = "F")
  m2["extrainfo"] <- c(C = "G", E = "H", I = "J")
  mydf1 <- to_df(c(m1, m2))
  mydf2 <- to_df(c(m1, m2), extrainfo = TRUE)
})
