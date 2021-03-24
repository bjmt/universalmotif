context("universalmotif_df")

test_that("to_df(), update_motifs(), to_list(), requires_update() work", {
  m1 <- create_motif()
  m2 <- create_motif()
  m3 <- create_motif()
  m4 <- create_motif()
  mydf <- to_df(c(m1, m2, m3, m4))
  expect_false(requires_update(mydf, extrainfo = FALSE))
  expect_equal(mydf$name[1], "motif")
  mydf$name <- LETTERS[1:4]
  expect_true(requires_update(mydf, extrainfo = FALSE))
  expect_equal(mydf$name[1], "A")
  mydf <- update_motifs(mydf, extrainfo = FALSE)
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
  mydf2$A[1] <- "K"
  expect_message(requires_update(mydf2, extrainfo = FALSE))
  expect_false(suppressMessages(requires_update(mydf2, extrainfo = FALSE)))
  expect_true(requires_update(mydf2, TRUE))

  # Check that holdout works
  # character isn't heldout, but factor & list should be
  mydf3 <- update_motifs(mydf2, extrainfo = FALSE)
  
  mydf3$list_column <- list("hello" = "darkness", "my old" = 0x667269656e64)
  mydf3$factor_column <- as.factor(letters[1:nrow(mydf3)])
  mydf3$char_column <- letters[1:nrow(mydf3)]
  mydf3$logical_column <- rep(TRUE, nrow(mydf3))

  mydf3_update <- update_motifs(mydf3, extrainfo = TRUE)
  expect_equal(mydf3$list_column, mydf3_update$list_column)
  expect_equal(mydf3$factor_column, mydf3_update$factor_column)
  expect_equal(mydf3$char_column, mydf3_update$char_column)
  expect_equal(mydf3$logical_column, mydf3_update$logical_column)
  # Ensure class is preserved
  expect_s3_class(mydf3, "universalmotif_df")
  expect_s3_class(mydf3_update, "universalmotif_df")
  expect_s3_class(update_motifs(mydf3, extrainfo = FALSE), "universalmotif_df")
  
  # test that force works
  forced <- update_motifs(mydf3, TRUE, TRUE)
  expect_type(forced$list_column, "character")
  expect_type(forced$factor_column, "character")
  expect_type(forced$char_column, "character")
  expect_type(forced$logical_column, "character")
  
  # to_list() and update_motifs() send message when necessary
  expect_message(update_motifs(mydf3, extrainfo = FALSE), "Discarding")
  expect_message(to_list(mydf3, extrainfo = FALSE), "Discarding")
  expect_message(update_motifs(mydf3, extrainfo = TRUE), NA)
  expect_message(to_list(mydf3, extrainfo = TRUE), NA)
})

test_that("update works", {
  m <- create_motif()
  m <- c(m, m)
  m[[1]]["extrainfo"] <- c(A = "2")
  m[[2]]["extrainfo"] <- c(B = "4")
  # Ensures rows aren't duplicated
  # Without extrainfo
  expect_equal(2, nrow(update_motifs(to_df(m, FALSE), FALSE)))
  expect_equal(2, nrow(update_motifs(to_df(m, TRUE), FALSE)))
  expect_equal(2, nrow(update_motifs(to_df(m, FALSE), TRUE)))
  # with extrainfo
  expect_equal(2, nrow(update_motifs(to_df(m, TRUE), TRUE)))
  
})
