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
  expect_message(requires_update(mydf2, extrainfo = FALSE), "Discarding")
  expect_true(suppressMessages(requires_update(mydf2, extrainfo = FALSE)))
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
  
  # Check that objects are correctly coerced to their character labels
  # In particular for factors to ensure they don't become characters of their numeric values.
  expect_equal(as.character(mydf3$list_column), forced$list_column)
  expect_equal(as.character(mydf3$factor_column), forced$factor_column)
  expect_equal(as.character(mydf3$char_column), forced$char_column)
  expect_equal(as.character(mydf3$logical_column), forced$logical_column)
  
  # to_list() and update_motifs() send message when necessary
  expect_message(update_motifs(mydf3, extrainfo = FALSE), "Discarding")
  expect_message(to_list(mydf3, extrainfo = FALSE), "Discarding")
  expect_message(update_motifs(mydf3, extrainfo = TRUE), NA)
  expect_message(to_list(mydf3, extrainfo = TRUE), NA)
  # Message when adding missing slots
  expect_message(update_motifs(mydf3[, -which(names(mydf3) == "altname")]), "Restoring")
})

test_that("to_list() rename with blank altname errors cleanly, not silently writing '0' (regression)", {

  m <- create_motif("ATCG", name = "mymotif")
  df <- to_df(m)
  expect_true(is.na(df$altname))
  df2 <- df
  df2$altname <- df$name   # altname <- "mymotif"
  df2$name <- NA_character_ # name <- NA (invalid: required slot)
  expect_error(suppressMessages(to_list(df2)), "name must be length 1")

})

test_that("to_list() correctly updates altname when both name and altname are present", {

  m <- create_motif("ATCG", name = "mymotif", altname = "myalt")
  df <- to_df(m)
  df$altname <- "newalt"
  m2 <- suppressMessages(to_list(df))
  if (is.list(m2)) m2 <- m2[[1]]
  expect_equal(m2@altname, "newalt")
  expect_equal(m2@name, "mymotif")

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

test_that("edits to the motif column are kept by update_motifs()", {
  m <- c(create_motif("AACCGGTT", name = "m1"),
         create_motif("AAACCCGT", name = "m2"))
  df <- to_df(m)
  ## A function that returns a new motif: its columns are refreshed quietly.
  df2 <- df
  df2$motif <- I(lapply(seq_along(m), function(i)
    if (i == 2) motif_rc(df$motif[[i]]) else df$motif[[i]]))
  expect_silent(u <- update_motifs(df2))
  expect_equal(u$consensus, c("AACCGGTT", "ACGGGTTT"))
  expect_equal(to_list(df2)[[2]]@consensus, "ACGGGTTT")
  ## Functions that modify the motif in place are no longer undone by the
  ## old values of its columns.
  df2$motif <- I(lapply(df$motif, convert_type, type = "PWM", nsize_correction = FALSE))
  expect_equal(vapply(suppressMessages(to_list(df2)), function(x) x@type,
                      character(1)), c("PWM", "PWM"))
  df2$motif <- I(lapply(df$motif, function(x) { x["name"] <- "new"; x }))
  expect_equal(update_motifs(df2)$name, c("new", "new"))
  ## Changes to the columns still reach the motifs, also in rows whose motif
  ## changed too.
  df2 <- df
  df2$name <- c("a", "b")
  df2$motif <- I(lapply(df$motif, trim_motifs, min.ic = 0.5))
  expect_equal(vapply(to_list(df2), function(x) x@name, character(1)), c("a", "b"))
  ## A column of a replaced motif that disagrees with it comes from the motif.
  df2 <- df
  df2$motif <- I(lapply(df$motif, function(x) { x <- motif_rc(x); x["name"] <- "z"; x }))
  expect_message(u <- update_motifs(df2), "'name'")
  expect_equal(u$name, c("z", "z"))
  ## Changing an unmodifiable column still warns.
  df2 <- df
  df2$consensus <- "AAAA"
  expect_warning(update_motifs(df2), "unmodifiable")
  ## Motifs leave the data.frame as they entered it.
  expect_identical(to_list(df), m)
  expect_null(attributes(to_list(df)[[1]])$universalmotif_df.row)
})
