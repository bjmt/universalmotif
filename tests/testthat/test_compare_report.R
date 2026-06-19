context("compare_motifs()/compare_motifs_lite() HTML report")

## ---- pure helpers (no rmarkdown needed) -----------------------------------

test_that(".cr_escape escapes markdown metacharacters", {
  expect_equal(universalmotif:::.cr_escape("a_b*c[d]#e"),
               "a\\_b\\*c\\[d\\]\\#e")
  expect_equal(universalmotif:::.cr_escape("plain"), "plain")
  expect_equal(universalmotif:::.cr_escape(NA), "")
})

test_that(".cr_check_output validates the path", {
  good <- file.path(tempdir(), "r.html")
  expect_equal(universalmotif:::.cr_check_output(good), good)
  expect_error(universalmotif:::.cr_check_output(file.path(tempdir(), "r.txt")),
               regexp = "html")
  expect_error(universalmotif:::.cr_check_output("/no/such/dir/r.html"),
               regexp = "does not exist")
})

test_that(".cr_build_table drops index columns and rounds stats", {
  m <- data.frame(query = "q", query.i = 3L, target = "t", target.i = 7L,
                  offset = 2L, strand = "+", score = 0.123456,
                  pvalue = 1.23e-5, qvalue = 2.34e-4,
                  query.i.local = 1L, target.i.local = 2L,
                  stringsAsFactors = FALSE)
  tab <- universalmotif:::.cr_build_table(m)
  expect_false(any(c("query.i", "target.i", "query.i.local",
                     "target.i.local") %in% names(tab)))
  expect_equal(tab$score, signif(0.123456, 3L))
})

test_that(".cr_build_rmd produces the expected structure", {
  m <- data.frame(query = "q", query.i = 3L, target = "t_1", target.i = 7L,
                  offset = 2L, strand = "+", score = 0.5,
                  query.i.local = 1L, target.i.local = 2L,
                  stringsAsFactors = FALSE)
  rmd <- universalmotif:::.cr_build_rmd(
    quote(compare_motifs_lite(x, compare.to = 1)),
    list(Method = "PCC", "Total matches" = 1), m, "view_motifs_lite",
    list(use.type = "ICM", tryRC = TRUE, min.overlap = 6L), "/tmp/x.rds")
  txt <- paste(rmd, collapse = "\n")
  expect_true(grepl("df_print: paged", txt, fixed = TRUE))
  expect_true(grepl("view_motifs_lite(.report$motifs[c(1, 2)]", txt, fixed = TRUE))
  expect_true(grepl(paste0("t", "\\", "_1"), txt, fixed = TRUE))   # escaped name
  expect_equal(sum(startsWith(rmd, "```{r logo-")), 1L)
})

test_that(".cr_build_rmd collapses a long wrapped call without gaps", {
  long_call <- quote(compare_motifs_lite(motifs = m, compare.to = 1, qvalue = 1,
                                     output.report = "r.html",
                                     output.report.max.print = 7))
  m <- data.frame(query = "q", query.i = 1L, target = "t", target.i = 2L,
                  offset = 0L, strand = "+", score = 1,
                  query.i.local = 1L, target.i.local = 2L,
                  stringsAsFactors = FALSE)
  rmd <- universalmotif:::.cr_build_rmd(long_call, list(), m, "view_motifs_lite",
                                        list(), "/tmp/x.rds")
  call_line <- grep("Call:", rmd, value = TRUE)
  expect_length(call_line, 1L)
  expect_false(grepl("  ", call_line, fixed = TRUE))   # no indentation gap
})

test_that(".cr_offsets_via_prep returns query-relative offset and strand", {
  skip_if_not_installed("GenomicRanges")  # only to gate optional infra; harmless
  q  <- create_motif("AAAACCCC", name = "q")
  ti <- create_motif("AAAACCCC", name = "same")
  rc <- create_motif("GGGGTTTT", name = "rc")     # = revcomp(AAAACCCC)
  mlist <- list(q, ti, rc)
  mats <- lapply(mlist, function(x) convert_type(x, "PPM")@motif)
  pr <- list(method = "PCC", tryRC = TRUE, min.overlap = 6, min.mean.ic = 0,
             min.position.ic = 0, relative_entropy = FALSE,
             normalise.scores = FALSE, score.strat = "a.mean")
  out <- universalmotif:::.cr_offsets_via_prep(
    mats, get_bkgs(mlist, 1), get_nsites(mlist), rownames(mats[[1]]),
    data.frame(qi = c(1L, 1L), ti = c(2L, 3L)), pr)
  expect_equal(out$offset[1], 0L)        # identical: aligned at 0
  expect_equal(out$strand[1], "+")
  expect_equal(out$strand[2], "-")       # reverse-complement match
})

## ---- end-to-end render (guarded) ------------------------------------------

report_motifs <- function() {
  list(create_motif("TTGACATA",   name = "query"),
       create_motif("TTGACATT",   name = "near1"),
       create_motif("ATTGACATAG", name = "near2"),
       create_motif("GGGGCCCC",   name = "far1"))
}

test_that("compare_motifs() writes a self-contained HTML report", {
  skip_on_cran()
  skip_if_not_installed("knitr")
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  out <- tempfile(fileext = ".html")
  on.exit(unlink(out))
  suppressWarnings(suppressMessages(
    compare_motifs(report_motifs(), compare.to = 1, max.p = 1, max.e = Inf,
                   output.report = out, output.report.max.print = 3)))
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_true(grepl("query", html))
  expect_true(grepl("offset", html))
  ## figures embedded as base64 data-URIs (truly self-contained), not
  ## broken external *_files references
  expect_match(html, 'src="data:image', fixed = TRUE)
  expect_false(grepl('_files/figure', html, fixed = TRUE))
  ## the call appears as clean inline code, not markdown-escaped
  expect_match(html, "compare_motifs(", fixed = TRUE)
})

test_that("compare_motifs_lite() writes a report with native offset/strand", {
  skip_on_cran()
  skip_if_not_installed("knitr")
  skip_if_not_installed("rmarkdown")
  skip_if_not(rmarkdown::pandoc_available())
  out <- tempfile(fileext = ".html")
  on.exit(unlink(out))
  suppressWarnings(suppressMessages(
    compare_motifs_lite(report_motifs(), compare.to = 1, qvalue = 1,
                    output.report = out, output.report.max.print = 3)))
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_true(grepl("query", html))
  expect_true(grepl("strand", html))
  expect_match(html, 'src="data:image', fixed = TRUE)
  expect_false(grepl('_files/figure', html, fixed = TRUE))
  expect_match(html, "compare_motifs_lite(", fixed = TRUE)
})

## ---- guards ----------------------------------------------------------------

test_that("output.report requires compare.to (both functions)", {
  m <- report_motifs()
  expect_error(
    suppressWarnings(compare_motifs(m, output.report = tempfile(fileext = ".html"))),
    regexp = "compare\\.to")
  expect_error(
    compare_motifs_lite(m, output.report = tempfile(fileext = ".html")),
    regexp = "compare\\.to")
})

test_that("output.report rejects a non-html path", {
  m <- report_motifs()
  expect_error(
    compare_motifs(m, compare.to = 1, max.p = 1, output.report = "x.txt"),
    regexp = "html")
})
