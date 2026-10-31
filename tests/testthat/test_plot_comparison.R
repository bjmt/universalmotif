context("plot_comparison()")

PC_MOTIFS <- read_jaspar(system.file("extdata", "jaspar.txt",
  package = "universalmotif"))
PC_NAMES <- vapply(PC_MOTIFS, function(m) m@name, character(1))

## The data of every layer drawn with the given geom.
pc_layers <- function(p, geom) {
  b <- ggplot2::ggplot_build(p)
  keep <- vapply(p$layers, function(l) class(l$geom)[1] == geom, logical(1))
  b$data[keep]
}

pc_tiles <- function(p) pc_layers(p, "GeomTile")[[1]]

## Row labels, top to bottom.
pc_rows <- function(p) {
  for (d in pc_layers(p, "GeomText")) {
    if (all(d$angle == 0) && all(d$hjust == 0) && !is.null(d$label) &&
        nrow(d) == length(PC_MOTIFS))
      return(as.character(d$label[order(-d$y)]))
  }
  NULL
}

pc_tree <- function(motifs = PC_MOTIFS, method = "complete") {
  s <- compare_motifs_lite(motifs, min.overlap = 6, matrix.out = "score")
  hclust(as.dist((1 - s) / 2), method = method)
}

pc_boxes <- function(p) {
  r <- pc_layers(p, "GeomRect")
  if (!length(r)) return(data.frame(start = integer(0), size = integer(0)))
  r <- r[[1]]
  out <- data.frame(start = round(r$xmin + 0.5),
    size = round(r$xmax - r$xmin))
  out[order(out$start), , drop = FALSE]
}

pc_runs <- function(tree, h) {
  g <- cutree(tree, h = h)[tree$order]
  runs <- rle(unname(g))
  ends <- cumsum(runs$lengths)
  keep <- runs$lengths > 1
  data.frame(start = (ends - runs$lengths + 1)[keep],
    size = runs$lengths[keep])
}

test_that("the order is that of the clustered scores", {
  for (method in c("complete", "average")) {
    tree <- pc_tree(method = method)
    p <- plot_comparison(PC_MOTIFS, method = method)
    expect_equal(pc_rows(p), PC_NAMES[tree$order])
  }
})

test_that("cells hold the comparison", {
  s <- compare_motifs_lite(PC_MOTIFS, min.overlap = 6, matrix.out = "score")
  ord <- pc_tree()$order
  p <- plot_comparison(PC_MOTIFS)
  n <- length(PC_MOTIFS)
  d <- p$layers[[1]]$data
  m <- matrix(NA_real_, n, n)
  m[cbind(n + 1 - d$y, d$x)] <- d$value
  expect_equal(m, unname(s[ord, ord]), tolerance = 1e-12)
  expect_equal(nrow(pc_tiles(p)), n * n)
  ## Every cell is written out, with two decimals.
  ann <- Filter(function(d) nrow(d) == n * n, pc_layers(p, "GeomText"))
  expect_length(ann, 1)
  expect_true(all(grepl("^-?[0-9]+\\.[0-9]{2}$", ann[[1]]$label)))
})

test_that("P-values are drawn as -log10 without the diagonal", {
  pv <- compare_motifs_lite(PC_MOTIFS, min.overlap = 6, matrix.out = "pvalue")
  ord <- pc_tree()$order
  p <- plot_comparison(PC_MOTIFS, value = "pvalue")
  n <- length(PC_MOTIFS)
  d <- p$layers[[1]]$data
  m <- matrix(NA_real_, n, n)
  m[cbind(n + 1 - d$y, d$x)] <- d$value
  want <- -log10(unname(pv[ord, ord]))
  diag(want) <- NA
  expect_equal(m, want, tolerance = 1e-12)
  ## The order still comes from the scores.
  expect_equal(pc_rows(p), PC_NAMES[ord])
  expect_equal(p$scales$get_scales("fill")$name, "-log10(P-value)")
  expect_equal(p$scales$get_scales("fill")$limits[1], 0)
})

test_that("a matrix can be given with or without the motifs", {
  s <- compare_motifs_lite(PC_MOTIFS, min.overlap = 6, matrix.out = "score")
  p <- plot_comparison(PC_MOTIFS, matrix = s)
  expect_equal(pc_rows(p), PC_NAMES[pc_tree()$order])
  p <- plot_comparison(matrix = s, logos = FALSE)
  expect_equal(pc_rows(p), colnames(s)[pc_tree()$order])
  expect_length(pc_layers(p, "GeomMotif"), 0)
  ## A P-value matrix alone is clustered on its symmetric -log10 values.
  pv <- compare_motifs_lite(PC_MOTIFS, min.overlap = 6, matrix.out = "pvalue")
  logs <- -log10(pv)
  sym <- (logs + t(logs)) / 2
  diag(sym) <- NA
  dist <- 1 - sym / max(sym, na.rm = TRUE)
  diag(dist) <- 0
  tree <- hclust(as.dist(dist))
  p <- plot_comparison(matrix = pv, value = "pvalue")
  expect_equal(pc_rows(p), colnames(pv)[tree$order])
})

test_that("cluster can be turned off or given as a tree", {
  p <- plot_comparison(PC_MOTIFS, cluster = FALSE)
  expect_equal(pc_rows(p), PC_NAMES)
  expect_length(pc_layers(p, "GeomSegment"), 0)
  tree <- pc_tree(method = "single")
  p <- plot_comparison(PC_MOTIFS, cluster = tree)
  expect_equal(pc_rows(p), PC_NAMES[tree$order])
})

test_that("dendrograms are brackets at the merge heights", {
  n <- length(PC_MOTIFS)
  tree <- pc_tree()
  p <- plot_comparison(PC_MOTIFS, dendrogram = "both")
  segs <- pc_layers(p, "GeomSegment")
  expect_length(segs, 2)
  expect_true(all(vapply(segs, nrow, 1L) == 3 * (n - 1)))
  ## The top tree's crossbars sit at heights in the order of the merges.
  top <- segs[[2]]
  across <- top[top$y == top$yend & top$x != top$xend, ]
  expect_equal(rank(across$y), rank(tree$height))
  expect_length(pc_layers(plot_comparison(PC_MOTIFS, dendrogram = "top"),
    "GeomSegment"), 1)
  expect_length(pc_layers(plot_comparison(PC_MOTIFS, dendrogram = "none"),
    "GeomSegment"), 0)
})

test_that("a cut boxes the groups and is drawn across the tree", {
  tree <- pc_tree()
  want <- pc_runs(tree, 0.3)
  expect_gt(nrow(want), 0)
  p <- plot_comparison(PC_MOTIFS, cut = 0.4, dendrogram = "both")
  expect_equal(pc_boxes(p), want, ignore_attr = TRUE)
  n.groups <- length(unique(cutree(tree, h = 0.3)))
  labs <- unlist(lapply(pc_layers(p, "GeomText"), function(d) d$label))
  expect_true(paste0("score 0.4, ", n.groups, " groups") %in% labs)
  ## The left line sits at the cut height and the top one across it.
  segs <- pc_layers(p, "GeomSegment")
  expect_length(segs, 4)
  ## A tree given as cluster is cut at a height, or at a score.
  expect_equal(pc_boxes(plot_comparison(PC_MOTIFS, cluster = tree,
    cut.height = 0.3)), want, ignore_attr = TRUE)
  expect_equal(pc_boxes(plot_comparison(PC_MOTIFS, cluster = tree,
    cut = 0.4)), want, ignore_attr = TRUE)
  ## Groups are boxed without a dendrogram too, and none at a score of 1.
  expect_equal(pc_boxes(plot_comparison(PC_MOTIFS, cut = 0.4,
    dendrogram = "none")), want, ignore_attr = TRUE)
  expect_equal(nrow(pc_boxes(plot_comparison(PC_MOTIFS, cut = 1))), 0)
})

test_that("a cut of a P-value matrix is a P-value", {
  pv <- compare_motifs_lite(PC_MOTIFS, min.overlap = 6, matrix.out = "pvalue")
  p <- plot_comparison(matrix = pv, value = "pvalue", cut = 0.05)
  labs <- unlist(lapply(pc_layers(p, "GeomText"), function(d) d$label))
  expect_true(any(grepl("^P-value 0.05, [0-9]+ groups?$", labs)))
  ## With motifs, the tree is one of scores and so is the cut.
  p <- plot_comparison(PC_MOTIFS, value = "pvalue", cut = 0.4)
  labs <- unlist(lapply(pc_layers(p, "GeomText"), function(d) d$label))
  expect_true(any(grepl("^score 0.4, ", labs)))
})

test_that("logos, annotation and the plot size", {
  p <- plot_comparison(PC_MOTIFS)
  expect_length(pc_layers(p, "GeomMotif"), 1)
  s <- attr(p, "size")
  expect_named(s, c("width", "height"))
  expect_true(all(s > 1 & s < 20))
  p <- plot_comparison(PC_MOTIFS, logos = FALSE, annotate = FALSE)
  expect_length(pc_layers(p, "GeomMotif"), 0)
  n <- length(PC_MOTIFS)
  expect_false(any(vapply(pc_layers(p, "GeomText"), nrow, 1L) == n * n))
  expect_lt(attr(p, "size")[["width"]], s[["width"]])
})

test_that("labels can be a slot or given", {
  p <- plot_comparison(PC_MOTIFS, labels = paste0("m", 1:5), cluster = FALSE)
  expect_equal(pc_rows(p), paste0("m", 1:5))
  dup <- lapply(PC_MOTIFS, function(m) { m@name <- "same"; m })
  p <- plot_comparison(dup, matrix = diag(5), cluster = FALSE)
  expect_equal(pc_rows(p), c("same", paste0("same.", 1:4)))
})

test_that("bad input is refused", {
  expect_error(plot_comparison(), "give `motifs`, a `matrix` or both")
  expect_error(plot_comparison(PC_MOTIFS, value = "evalue"))
  expect_error(plot_comparison(PC_MOTIFS, dendrogram = "right"))
  expect_error(plot_comparison(matrix = matrix(0, 2, 3)), "must be square")
  expect_error(plot_comparison(PC_MOTIFS[1:3], matrix = diag(4)),
    "one row and column per motif")
  expect_error(plot_comparison(PC_MOTIFS[1]), "at least two motifs")
  expect_error(plot_comparison(PC_MOTIFS[1:3], labels = c("a", "b")),
    "one label per motif")
  expect_error(plot_comparison(matrix = diag(3), logos = TRUE),
    "needs the motifs")
  expect_error(plot_comparison(PC_MOTIFS[1:4], cluster = pc_tree(PC_MOTIFS[1:3])),
    "one tip per motif")
  expect_error(plot_comparison(PC_MOTIFS, cluster = "yes"), "hclust object")
  expect_error(plot_comparison(PC_MOTIFS, cut = 0.5, cut.height = 0.2),
    "not both")
  expect_error(plot_comparison(PC_MOTIFS, cluster = FALSE, cut = 0.5),
    "cutting needs a tree")
  expect_error(plot_comparison(matrix = matrix(0.5, 3, 3), value = "pvalue",
    cut = 0), "above 0 and at most 1")
})
