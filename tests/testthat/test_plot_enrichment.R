context("plot_enrichment()")

## A table shaped like enrich_motifs_lite() output, with a group of strongly
## enriched motifs and a cloud of unenriched ones.
fake_enrich <- function(n = 40, seed = 1) {
  set.seed(seed)
  k <- n %/% 4
  x <- c(rnorm(k, 2, 0.3), rnorm(n - k, 0, 0.4))
  lq <- -c(runif(k, 2, 6), runif(n - k, 0, 0.9))
  data.frame(motif = paste0("m", seq_len(n) - 1L), motif.i = seq_len(n),
    target.seq.hits = sample.int(79, n, TRUE),
    target.site.hits = sample.int(119, n, TRUE),
    log2.enrichment = x, pvalue = 10^(lq - 0.5), qvalue = 10^lq,
    stringsAsFactors = FALSE)
}

layer_of <- function(p, geom) {
  i <- which(vapply(p$layers, function(l) class(l$geom)[1], "") == geom)
  if (!length(i)) return(NULL)
  ggplot2::ggplot_build(p)$data[[i[1]]]
}

label_data <- function(p) {
  d <- layer_of(p, "GeomEnrichLabels")
  if (is.null(d)) return(character(0))
  d[!is.na(d$label), , drop = FALSE]
}

find_grobs <- function(g, cls) {
  if (inherits(g, cls)) return(list(g))
  kids <- if (inherits(g, "gtable")) g$grobs
    else if (inherits(g, "gTree")) g$children else list()
  unlist(lapply(kids, find_grobs, cls), recursive = FALSE)
}

## Place the labels of a plot in a panel of the given size (in points) and
## return the grobs drawn.
drawn_labels <- function(p, width = 400, height = 300) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  lab <- find_grobs(ggplot2::ggplotGrob(p), "um_enrich_labels")[[1]]
  grid::pushViewport(grid::viewport(width = grid::unit(width, "pt"),
    height = grid::unit(height, "pt")))
  kids <- makeContent.um_enrich_labels(lab)$children
  boxes <- lapply(kids, function(k) {
    if (inherits(k, "text")) {
      w <- grid::convertWidth(grid::grobWidth(k), "pt", valueOnly = TRUE)
      x <- grid::convertX(k$x, "pt", valueOnly = TRUE)
      y <- grid::convertY(k$y, "pt", valueOnly = TRUE)
      h <- 1.2 * k$gp$fontsize
      c(x, y - h / 2, x + w, y + h / 2)
    } else if (inherits(k, "polygon")) {
      x <- grid::convertX(k$x, "pt", valueOnly = TRUE)
      y <- grid::convertY(k$y, "pt", valueOnly = TRUE)
      c(min(x), min(y), max(x), max(y))
    } else NULL
  })
  type <- vapply(kids, function(k) class(k)[1], "")
  list(kids = kids, type = type, boxes = boxes)
}

overlaps <- function(a, b) a[1] < b[3] && b[1] < a[3] && a[2] < b[4] &&
  b[2] < a[4]

test_that("points are the table values", {
  res <- fake_enrich()
  p <- plot_enrichment(res)
  expect_s3_class(p, "ggplot")
  d <- layer_of(p, "GeomPoint")
  expect_equal(nrow(d), nrow(res))
  got <- d[order(d$x), c("x", "y")]
  want <- data.frame(x = res$log2.enrichment, y = -log10(res$qvalue))
  want <- want[order(want$x), ]
  expect_equal(unname(as.matrix(got)), unname(as.matrix(want)))
  expect_equal(p$labels$x, "log2(enrichment)")
  expect_equal(p$labels$y, "-log10(q-value)")
})

test_that("P-values are drawn on request and underflow stays finite", {
  res <- fake_enrich()
  res$qvalue[1] <- 0
  p <- plot_enrichment(res)
  d <- layer_of(p, "GeomPoint")
  expect_equal(max(d$y), -log10(5e-324))
  expect_gt(p$coordinates$limits$y[2], max(d$y))
  p <- plot_enrichment(res, y = "pvalue")
  d <- layer_of(p, "GeomPoint")
  expect_equal(sort(d$y), sort(-log10(res$pvalue)))
  expect_equal(p$labels$y, "-log10(P-value)")
})

test_that("enrich_motifs() tables are accepted", {
  res <- fake_enrich()
  old <- S4Vectors::DataFrame(motif = res$motif,
    target.hits = res$target.site.hits,
    target.seq.hits = res$target.seq.hits,
    Pval = res$pvalue, Qval = res$qvalue,
    target.enrichment = 2^res$log2.enrichment)
  p <- plot_enrichment(old, size.by = "target.site.hits")
  d <- layer_of(p, "GeomPoint")
  expect_equal(sort(d$x), sort(res$log2.enrichment))
  expect_equal(sort(d$y), sort(-log10(res$qvalue)))
  expect_equal(p$scales$get_scales("size")$name, "Target hits")
})

test_that("significant points are coloured and the cutoff is marked", {
  res <- fake_enrich()
  p <- plot_enrichment(res, cutoff = 0.01, colour = "#FF0000")
  d <- layer_of(p, "GeomPoint")
  expect_true(all(d$fill[d$y >= 2] == "#FF0000"))
  expect_true(all(d$fill[d$y < 2] == "#B4B4B4"))
  h <- layer_of(p, "GeomHline")
  expect_equal(h$yintercept, 2)
  expect_equal(h$linetype, "dashed")
  txt <- layer_of(p, "GeomText")
  expect_equal(txt$label, "q = 0.01")
  p <- plot_enrichment(res, y = "pvalue", cutoff = 0.001)
  expect_equal(layer_of(p, "GeomText")$label, "p = 0.001")
})

test_that("point sizes follow the size column", {
  res <- fake_enrich()
  p <- plot_enrichment(res, size.by = "target.site.hits")
  d <- layer_of(p, "GeomPoint")
  hits <- res$target.site.hits[match(round(d$x, 10),
    round(res$log2.enrichment, 10))]
  o <- order(hits)
  expect_true(all(diff(d$size[o]) >= -1e-12))
  expect_lt(d$size[which.min(hits)], d$size[which.max(hits)])
  expect_equal(p$scales$get_scales("size")$name, "Target hits")
  p <- plot_enrichment(res, size.by = NULL)
  d <- layer_of(p, "GeomPoint")
  expect_length(unique(d$size), 1L)
  expect_null(p$scales$get_scales("size"))
  p <- plot_enrichment(res)
  expect_equal(p$scales$get_scales("size")$name,
    "Target sequences\nwith hits")
})

test_that("the most significant motifs are labelled", {
  res <- fake_enrich()
  lab <- label_data(plot_enrichment(res, labels = 5))
  top <- res$motif[order(res$qvalue)][1:5]
  expect_setequal(lab$label, top)
  expect_equal(lab$label[order(lab$rank)], top)
  ## Each label points at its motif.
  i <- match(lab$label, res$motif)
  expect_equal(lab$x, res$log2.enrichment[i])
  expect_equal(lab$y, -log10(res$qvalue[i]))
  p <- plot_enrichment(res, labels = 0)
  expect_null(layer_of(p, "GeomEnrichLabels"))
  ## Only significant motifs are labelled by number.
  lab <- label_data(plot_enrichment(res, labels = 100, cutoff = 1e-3))
  expect_equal(nrow(lab), sum(res$qvalue <= 1e-3))
})

test_that("motifs can be labelled by name", {
  res <- fake_enrich()
  expect_setequal(label_data(plot_enrichment(res,
    labels = c("m30", "m31")))$label, c("m30", "m31"))
  expect_equal(label_data(plot_enrichment(res, labels = "m30"))$label, "m30")
  expect_error(plot_enrichment(res, labels = c("m1", "nope")),
    "no motifs named nope")
})

test_that("labels do not overlap each other or leave the panel", {
  res <- fake_enrich(n = 60, seed = 3)
  out <- drawn_labels(plot_enrichment(res, labels = 12))
  boxes <- out$boxes[out$type == "text"]
  expect_length(boxes, 12L)
  for (i in seq_along(boxes)) {
    for (j in seq_along(boxes)) {
      if (i < j) expect_false(overlaps(boxes[[i]], boxes[[j]]))
    }
    expect_true(boxes[[i]][1] >= 0 && boxes[[i]][3] <= 400 &&
      boxes[[i]][2] >= 0 && boxes[[i]][4] <= 300)
  }
})

test_that("label placement nudges labels apart", {
  ## Three points on top of each other with labels of the same size.
  px <- py <- rep(100, 3)
  spots <- enrich_place_labels(px, py, rep(3, 3), rep(30, 3), rep(10, 3),
    1:3, c(0, 0, 300, 300))
  boxes <- lapply(1:3, function(i) c(spots$dx[i], spots$dy[i],
    spots$dx[i] + 30, spots$dy[i] + 10))
  expect_false(overlaps(boxes[[1]], boxes[[2]]))
  expect_false(overlaps(boxes[[1]], boxes[[3]]))
  expect_false(overlaps(boxes[[2]], boxes[[3]]))
  ## A lone label sits next to its point, with no leader line.
  one <- enrich_place_labels(100, 100, 3, 30, 10, 1L, c(0, 0, 300, 300))
  expect_false(one$far)
  expect_lt(sqrt(one$dx^2 + (one$dy + 5)^2), 10)
  ## Every label not right beside its point gets a leader line, so that it
  ## is not read as belonging to a nearer point.
  spots <- enrich_place_labels(rep(100, 12), rep(100, 12), rep(3, 12),
    rep(30, 12), rep(10, 12), 1:12, c(0, 0, 300, 300))
  expect_true(any(spots$far) && !all(spots$far))
})

test_that("far labels get leader lines", {
  res <- fake_enrich(n = 60, seed = 3)
  out <- drawn_labels(plot_enrichment(res, labels = 12), 200, 150)
  expect_true(any(out$type == "segments"))
})

test_that("motifs with infinite enrichment sit at the edges", {
  res <- fake_enrich()
  res$log2.enrichment[36] <- -Inf
  res$log2.enrichment[2] <- Inf
  p <- plot_enrichment(res)
  pts <- layer_of(p, "GeomPoint")
  edges <- layer_of(p, "GeomEnrichEdges")
  expect_equal(nrow(edges), 2L)
  expect_equal(nrow(pts), nrow(res) - 2L)
  lim <- p$coordinates$limits$x
  lo <- edges[edges$direction == -1, ]
  hi <- edges[edges$direction == 1, ]
  expect_true(lim[1] < lo$x && lo$x < min(pts$x))
  expect_true(max(pts$x) < hi$x && hi$x < lim[2])
  expect_equal(lo$y, -log10(res$qvalue[36]))
  expect_equal(hi$y, -log10(res$qvalue[2]))
  expect_true(all(is.finite(c(edges$x, lim))))
})

test_that("logos are drawn beside their labels", {
  res <- fake_enrich()
  motifs <- lapply(0:9, function(i) create_motif("TTGACCA",
    name = paste0("m", i), nsites = 20))
  out <- drawn_labels(plot_enrichment(res, labels = 4, motifs = motifs))
  logos <- out$boxes[out$type == "polygon"]
  texts <- out$boxes[out$type == "text"]
  expect_length(logos, 4L)
  for (b in logos) {
    expect_gt(b[3] - b[1], 0)
    expect_false(any(vapply(texts, overlaps, logical(1), b)))
    ## Each logo sits level with one label, close beside it.
    level <- vapply(texts, function(t) {
      if (t[2] < b[4] && b[2] < t[4]) max(t[1] - b[3], b[1] - t[3]) else Inf
    }, numeric(1))
    expect_true(min(level) > 0 && min(level) < 25)
  }
  expect_error(plot_enrichment(res, labels = 4, motifs = motifs[1:2]),
    "has no motif named")
  ## A named list is matched by its names.
  top <- res$motif[which.min(res$qvalue)]
  named <- stats::setNames(list(create_motif("CACGTG")), top)
  out <- drawn_labels(plot_enrichment(res, labels = 1, motifs = named))
  expect_equal(sum(out$type == "polygon"), 1L)
})

test_that("bad input is refused", {
  res <- fake_enrich()
  expect_error(plot_enrichment(as.list(res)), "data.frame")
  expect_error(plot_enrichment(res[, c("motif", "qvalue")]),
    "log2.enrichment")
  expect_error(plot_enrichment(res[, setdiff(names(res), "qvalue")]),
    "no `qvalue` column")
  expect_error(plot_enrichment(res, y = "evalue"))
  expect_error(plot_enrichment(res, cutoff = 0), "cutoff")
  expect_error(plot_enrichment(res, cutoff = 2), "cutoff")
  expect_error(plot_enrichment(res, size.by = "nope"), "no `nope` column")
  expect_error(plot_enrichment(res, labels = -1), "labels")
  expect_error(plot_enrichment(res[0, ]), "no rows")
  res$log2.enrichment <- NA
  expect_error(plot_enrichment(res), "no enrichment values")
})

test_that("enrich_motifs_lite() results can be plotted", {
  data(ArabidopsisPromoters)
  data(ArabidopsisMotif)
  motifs <- list(ArabidopsisMotif,
    create_motif("TTGACY", name = "W-box", nsites = 20),
    create_motif("CACGTG", name = "G-box", nsites = 20))
  bkg <- shuffle_sequences(ArabidopsisPromoters, k = 2, rng.seed = 1)
  res <- enrich_motifs_lite(motifs, ArabidopsisPromoters, bkg,
    pvalue = 1e-3, qvalue = 1, nthreads = 1)
  p <- plot_enrichment(res, labels = 2, motifs = motifs)
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
  out <- drawn_labels(p)
  expect_equal(sum(out$type == "text"), nrow(label_data(p)))
})
