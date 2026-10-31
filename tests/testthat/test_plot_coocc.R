context("plot_coocc()")

## A table shaped like motif_coocc() output: A and B co-occur strongly, A and
## D never miss each other (an infinite odds ratio), C and D never co-occur,
## and B and D were not tested.
fake_coocc <- function() {
  data.frame(
    motif_a = c("A", "A", "A", "B", "B", "C"),
    motif_b = c("B", "C", "D", "C", "D", "D"),
    odds_ratio = c(8, 0.5, Inf, 1, NA, 0),
    pvalue = c(1e-6, 0.9, 1e-3, 0.5, NA, 1),
    qvalue = c(1e-5, 0.9, 0.05, 0.7, NA, 1),
    stringsAsFactors = FALSE)
}

layer_of <- function(p, geom) {
  i <- which(vapply(p$layers, function(l) class(l$geom)[1], "") == geom)
  ggplot2::ggplot_build(p)$data[[i[1]]]
}

cell <- function(p, a, b) {
  d <- p$data
  d$value[as.character(d$a) == a & as.character(d$b) == b]
}

test_that("cells hold log2 odds ratios, symmetric and grey where untested", {
  res <- fake_coocc()
  p <- plot_coocc(res, order = "input")
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$a), c("A", "B", "C", "D"))
  ## The first motif is at the top.
  expect_equal(levels(p$data$b), c("D", "C", "B", "A"))
  expect_equal(nrow(p$data), 16L)
  expect_equal(cell(p, "A", "B"), 3)
  expect_equal(cell(p, "B", "A"), 3)
  expect_equal(cell(p, "A", "C"), -1)
  expect_true(is.na(cell(p, "B", "D")))
  expect_true(is.na(cell(p, "A", "A")))
  ## Infinite odds ratios sit at the ends of the scale.
  expect_equal(cell(p, "A", "D"), 3)
  expect_equal(cell(p, "C", "D"), -3)
  expect_equal(p$scales$get_scales("fill")$limits, c(-3, 3))
  expect_equal(p$scales$get_scales("fill")$name, "log2(odds ratio)")
  tiles <- layer_of(p, "GeomTile")
  expect_equal(nrow(tiles), 16L)
  expect_equal(sum(tiles$fill == "#D0D0D0"), 6L)
})

test_that("pairs at or below the cutoff are marked", {
  res <- fake_coocc()
  p <- plot_coocc(res, order = "input")
  pts <- layer_of(p, "GeomPoint")
  ## A-B and A-D, on both sides of the diagonal.
  expect_equal(nrow(pts), 4L)
  marked <- p$layers[[2]]$data
  expect_setequal(paste(marked$a, marked$b),
    c("A B", "B A", "A D", "D A"))
  ## Dots are white on the dark cells at the end of the scale.
  expect_true(all(marked$colour == "#FFFFFF"))
  p <- plot_coocc(res, cutoff = 1e-4, order = "input")
  marked <- p$layers[[2]]$data
  expect_setequal(paste(marked$a, marked$b), c("A B", "B A"))
  p <- plot_coocc(res, cutoff = 1e-9)
  expect_equal(nrow(p$layers[[2]]$data), 0L)
  expect_no_error(ggplot2::ggplot_build(p))
  ## Dots on light cells are dark.
  expect_equal(coocc_dot_colours(c(0, 3, NA), c(-3, 3),
    c("#053061", "#F7F7F7", "#67001F")), c("#1A1A1A", "#FFFFFF", "#1A1A1A"))
})

test_that("cells can show q-values or P-values", {
  res <- fake_coocc()
  p <- plot_coocc(res, value = "qvalue", order = "input")
  expect_equal(cell(p, "A", "B"), 5)
  expect_equal(p$scales$get_scales("fill")$limits, c(0, 5))
  expect_equal(p$scales$get_scales("fill")$name, "-log10(q-value)")
  p <- plot_coocc(res, value = "pvalue", order = "input")
  expect_equal(cell(p, "A", "B"), 6)
  expect_equal(p$scales$get_scales("fill")$name, "-log10(P-value)")
  p <- plot_coocc(res, colours = c("white", "black"), value = "pvalue",
    order = "input")
  tiles <- layer_of(p, "GeomTile")
  ## A-B has the smallest P-value, and C-D a P-value of one.
  expect_equal(toupper(tiles$fill[tiles$x == 1 & tiles$y == 3]), "#000000")
  expect_equal(toupper(tiles$fill[tiles$x == 3 & tiles$y == 1]), "#FFFFFF")
})

test_that("motifs can be ordered by clustering, input or name", {
  ## Motifs given apart in the table, two pairs found together.
  res <- data.frame(
    motif_a = c("A", "A", "A", "B", "B", "C"),
    motif_b = c("B", "C", "D", "C", "D", "D"),
    odds_ratio = c(0.5, 16, 0.6, 0.4, 20, 0.7),
    pvalue = c(0.9, 1e-8, 0.8, 0.9, 1e-9, 0.7),
    qvalue = c(0.9, 1e-7, 0.9, 0.9, 1e-8, 0.9))
  p <- plot_coocc(res)
  ord <- levels(p$data$a)
  expect_equal(abs(diff(match(c("A", "C"), ord))), 1L)
  expect_equal(abs(diff(match(c("B", "D"), ord))), 1L)
  expect_equal(levels(plot_coocc(res, order = "input")$data$a),
    c("A", "B", "C", "D"))
  p <- plot_coocc(res, order = c("D", "A", "C"))
  expect_equal(levels(p$data$a), c("D", "A", "C"))
  expect_equal(nrow(p$data), 9L)
  expect_error(plot_coocc(res, order = c("A", "nope")), "no motifs named nope")
})

test_that("bad input is refused", {
  res <- fake_coocc()
  expect_error(plot_coocc(as.list(res)), "data.frame")
  expect_error(plot_coocc(res[, -5]), "no `qvalue` column")
  expect_error(plot_coocc(res, value = "both"))
  expect_error(plot_coocc(res, cutoff = 0), "cutoff")
  expect_error(plot_coocc(res, order = 1:3), "order")
  expect_error(plot_coocc(res, colours = "red"), "colours")
  expect_error(plot_coocc(res[0, ]), "no rows")
})

test_that("motif_coocc() results can be plotted", {
  data(ArabidopsisPromoters)
  data(ArabidopsisMotif)
  motifs <- list(ArabidopsisMotif,
    create_motif("TTGACY", name = "W-box", nsites = 20),
    create_motif("CACGTG", name = "G-box", nsites = 20))
  res <- motif_coocc(motifs, ArabidopsisPromoters, pvalue = 1e-3,
    nthreads = 1)
  p <- plot_coocc(res)
  expect_s3_class(p, "ggplot")
  expect_setequal(levels(p$data$a), c(ArabidopsisMotif@name, "W-box",
    "G-box"))
  expect_no_error(ggplot2::ggplot_build(p))
})
