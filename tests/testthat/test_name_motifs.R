context("name_motifs()")

## Reference values for the fixture motifs of helper-lite-fixture.R were
## worked out independently.

test_that("name_motifs() agrees with the reference values", {
  ms <- lite_fixture()
  db <- ms[c(2, 3, 4, 5, 7)]
  db[[1]]@altname <- "FAC2"
  res <- name_motifs(ms[c(1, 6)], db, top = 3, qvalue = 1)
  expect_identical(colnames(res),
    c("motif", "motif.i", "rank", "match", "altname", "score", "pvalue",
      "qvalue", "strand", "offset", "consensus", "match.consensus"))
  expect_equal(res$motif, rep(c("M1", "M6"), each = 3))
  expect_equal(res$motif.i, rep(1:2, each = 3))
  expect_equal(res$rank, rep(1:3, 2))
  expect_equal(res$match, c("M2", "M3", "M7", "M7", "M3", "M2"))
  expect_equal(res$altname, c("FAC2", NA, NA, NA, NA, "FAC2"))
  expect_equal(res$score, c(0.936792, 0.78106, 0.66503, 0.785782, 0.610347,
                            0.471998), tolerance = 1e-5)
  expect_equal(res$pvalue, c(0.000155432, 0.00176302, 0.00599699, 0.001479,
                             0.0385364, 0.174444), tolerance = 1e-4)
  expect_equal(res$qvalue, c(0.000777, 0.004408, 0.009995, 0.007395, 0.096341,
                             0.29074), tolerance = 1e-3)
  expect_equal(res$strand, c("+", "+", "+", "-", "-", "-"))
  expect_equal(res$offset, c(-1L, 1L, 1L, 1L, 1L, -1L))
  expect_equal(res$consensus[1], "WKGACATA")
  expect_equal(res$match.consensus[1], "TYGACATA")
})

test_that("name_motifs() corrects q-values for the size of the database only", {
  ms <- lite_fixture()
  db <- ms[c(2, 3, 4, 5, 7)]
  one <- name_motifs(ms[1], db, qvalue = 1, top = 5)
  two <- name_motifs(ms[c(1, 6)], db, qvalue = 1, top = 5)
  expect_equal(nrow(one), 5L)
  p <- one$pvalue
  expect_equal(one$qvalue, pmin(rev(cummin(rev(p * 5 / seq_along(p)))), 1))
  expect_equal(sum(two$motif.i == 1L), 5L)
})

test_that("name_motifs() keeps at most `top` matches passing `qvalue`", {
  ms <- lite_fixture()
  db <- ms[c(2, 3, 4, 5, 7)]
  res <- name_motifs(ms[c(1, 6)], db)
  expect_equal(res$match, c("M2", "M3", "M7", "M7"))
  expect_true(all(res$qvalue <= 0.05))
  res <- name_motifs(ms[c(1, 6)], db, top = 1, qvalue = 1)
  expect_equal(res$match, c("M2", "M7"))
  none <- name_motifs(ms[4], db[-3], qvalue = 0.01)
  expect_equal(nrow(none), 0L)
  expect_identical(colnames(none), colnames(res))
})

test_that("name_motifs() checks its arguments", {
  ms <- lite_fixture()
  expect_error(name_motifs(ms[1], list()), "at least one")
  expect_error(name_motifs(list(), ms), "at least one")
  expect_error(name_motifs(ms[1], ms, top = 0), "`top`")
  expect_error(name_motifs(ms[1], ms, qvalue = 2), "`qvalue`")
})

test_that("name_motifs() names motifs sharing names with the database", {
  db <- list(create_motif("TGACTCA", name = "MA0099", altname = "FOS::JUN"),
             create_motif("AGATAA", name = "MA0035", altname = "GATA1"),
             create_motif("CACGTG", name = "MA0004", altname = "Arnt"))
  found <- list(create_motif("ATGAGTCAT", name = "MA0035"),
                create_motif("CTTATCT", name = "pos_2"))
  res <- name_motifs(found, db, top = 1, qvalue = 1)
  expect_equal(res$motif, c("MA0035", "pos_2"))
  expect_equal(res$match, c("MA0099", "MA0035"))
  expect_equal(res$altname, c("FOS::JUN", "GATA1"))
})
