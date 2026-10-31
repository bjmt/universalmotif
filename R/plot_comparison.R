#' Draw motif similarities as a clustered heatmap.
#'
#' Motifs are compared all against all with [compare_motifs_lite()], or the
#' comparison is given as `matrix`, a square matrix such as
#' `compare_motifs_lite()` returns without `compare.to`. Rows are queries and
#' columns targets, both in the order of a hierarchical clustering drawn as a
#' dendrogram beside the heatmap. Scores are drawn as they are, and P- and
#' q-values as their -log10, with the diagonal left blank since it compares
#' each motif with itself.
#'
#' @param motifs See [convert_motifs()] for accepted motif formats. The motifs
#'   to compare, or the motifs `matrix` was made from, which give the labels
#'   and logos.
#' @param matrix `matrix` A precomputed square matrix of the kind chosen by
#'   `value`. Its column names become the labels when no motifs are given.
#'   Scores are clustered as correlations, so a matrix of other similarity
#'   scores needs a tree of its own as `cluster`, or no clustering.
#' @param labels `character` A motif slot to label rows and columns with,
#'   such as `"name"` or `"altname"`, or one label per motif.
#' @param value `character(1)` What is drawn: `"score"`, the mean column
#'   correlation, or `"pvalue"` or `"qvalue"`.
#' @param min.overlap `integer(1)` Passed to [compare_motifs_lite()]. The
#'   default of 6 follows [motif_tree_lite()], so the scores can differ
#'   slightly from those of `compare_motifs_lite()` with its default of 5.
#' @param RC `logical(1)` Passed to [compare_motifs_lite()].
#' @param cluster `logical(1)` or `hclust` Order the motifs by clustering them
#'   with `method`, or keep them in input order with `FALSE`. A tree from
#'   [stats::hclust()] can be given instead. Motifs are clustered on their
#'   scores, with distances of `(1 - score) / 2`, so that the order is that of
#'   [motif_tree_lite()] with the same `min.overlap`. A matrix of P- or
#'   q-values given as `matrix` without motifs is clustered on the mean
#'   -log10 value of both directions, scaled to distances from 0 to 1.
#' @param method `character(1)` The linkage passed to [stats::hclust()].
#' @param dendrogram `character(1)` Where to draw the tree: `"left"`,
#'   `"top"`, `"both"` or `"none"`.
#' @param cut,cut.height `numeric(1)` Show where cutting the tree would split
#'   the motifs into groups, as a dashed line across the dendrogram and a box
#'   around each group of two or more motifs on the heatmap. `cut` is in the
#'   units clustered on, a score such as 0.8 (or a P- or q-value for such a
#'   matrix given without motifs), and with complete linkage every motif in a
#'   group is then at least that similar to every other. A tree given as
#'   `cluster` is taken to be one of scores; `cut.height` instead cuts at a
#'   height of the tree itself, as [stats::cutree()] does.
#' @param colours `character` The colours of the gradient, from low to high.
#'   The default is the ColorBrewer "Blues" palette.
#' @param limits `numeric(2)` The range of the colour scale, which by default
#'   spans the values drawn, from 0 for P- and q-values. Values outside it
#'   take the colour at the nearest end.
#' @param annotate `logical(1)` Write each value in its cell, by default for
#'   up to 20 motifs.
#' @param logos `logical(1)` Draw each motif's logo beside its row label, by
#'   default when motifs are given and there are up to 60 of them.
#' @param text.size `numeric(1)` The largest font size of labels, in points.
#'   Labels shrink with the cells of larger heatmaps and are left out once
#'   they would be unreadable.
#' @param nthreads `numeric(1)` Passed to [compare_motifs_lite()].
#'
#' @return A `ggplot` object, with the width and height in inches it is drawn
#'   for as its `"size"` attribute (see Details).
#'
#' @details
#' Everything is drawn in a single panel with a fixed aspect ratio, where
#' each cell is one unit across, so that the dendrograms, logos and labels
#' stay aligned with the rows and columns whatever the size of the device.
#' Text sizes are fixed in points, however, and are chosen for cells of about
#' 0.45 inches for up to 15 motifs, shrinking to about 0.15 inches for 80 or
#' more. The size of device this needs, labels and colour bar included, is
#' given by `attr(p, "size")`, which can be passed on to [ggplot2::ggsave()]
#' as in the examples.
#'
#' @examples
#' motifs <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' p <- plot_comparison(motifs)
#' p
#' \dontrun{
#' ggplot2::ggsave("comparison.pdf", p, width = attr(p, "size")[1],
#'   height = attr(p, "size")[2])
#' }
#'
#' ## P-values, with the tree cut at a mean column correlation of 0.4:
#' plot_comparison(motifs, value = "pvalue", dendrogram = "both", cut = 0.4)
#'
#' @seealso [compare_motifs_lite()], [motif_tree_lite()], [geom_motif()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_comparison <- function(motifs = NULL, matrix = NULL, labels = "name",
  value = c("score", "pvalue", "qvalue"), min.overlap = 6L, RC = TRUE,
  cluster = TRUE, method = "complete",
  dendrogram = c("left", "top", "both", "none"), cut = NULL,
  cut.height = NULL, colours = PLOT_COMPARISON_BLUES, limits = NULL,
  annotate = NULL, logos = NULL, text.size = 10, nthreads = 1) {

  value <- match.arg(value)
  dendrogram <- match.arg(dendrogram)
  input <- comparison_input(motifs, matrix, labels, value, min.overlap, RC,
    nthreads)
  motifs <- input$motifs
  vals <- input$vals
  scores <- input$scores
  labs <- input$labels
  n <- length(labs)
  if (isTRUE(logos) && is.null(motifs))
    stop("`logos` needs the motifs", call. = FALSE)

  tree <- NULL
  if (inherits(cluster, "hclust")) {
    if (length(cluster$order) != n)
      stop("the tree must have one tip per motif", call. = FALSE)
    tree <- cluster
  } else if (!is.logical(cluster) || length(cluster) != 1 || is.na(cluster)) {
    stop("`cluster` must be TRUE, FALSE or an hclust object", call. = FALSE)
  } else if (cluster) {
    tree <- stats::hclust(stats::as.dist(comparison_distances(vals, scores)),
      method = method)
  }
  ord <- if (is.null(tree)) seq_len(n) else tree$order

  cut.at <- NULL
  if (!is.null(cut) && !is.null(cut.height))
    stop("give `cut` or `cut.height`, not both", call. = FALSE)
  if (!is.null(cut) || !is.null(cut.height)) {
    if (is.null(tree))
      stop("cutting needs a tree, so `cluster` cannot be FALSE", call. = FALSE)
    ## The cutoff is in the units clustered on, which are scores unless only
    ## a matrix of P- or q-values was given.
    cut.unit <- if (!is.null(scores) || inherits(cluster, "hclust")) "score"
      else value
    cut.at <- if (!is.null(cut.height)) as.numeric(cut.height)
      else comparison_cut_height(as.numeric(cut), vals, cut.unit)
    groups <- stats::cutree(tree, h = cut.at)
  }

  shown <- if (value == "score") vals else neg_log10(vals)
  shown <- shown[ord, ord, drop = FALSE]
  if (value != "score") diag(shown) <- NA
  finite <- shown[is.finite(shown)]
  if (is.null(limits)) {
    lo <- if (value != "score") 0 else if (length(finite)) min(finite) else 0
    hi <- if (length(finite)) max(finite) else 1
    if (hi <= lo) hi <- lo + 1
    limits <- c(lo, hi)
  }

  ## Cells shrink as motifs are added, down to a floor, and the labels with
  ## them. Sizes below are in cells, with `cell` inches to a cell.
  cell <- max(min(0.45, 7 / n), min(0.15, 12 / n))
  label.size <- min(text.size, 0.75 * cell * 72)
  show.labels <- label.size >= 4
  if (is.null(annotate)) annotate <- n <= 20
  if (is.null(logos)) logos <- !is.null(motifs) && n <= 60
  left.tree <- !is.null(tree) && dendrogram %in% c("left", "both")
  top.tree <- !is.null(tree) && dendrogram %in% c("top", "both")
  tree.w <- min(max(0.15 * n * cell, 0.6), 1.5) / cell
  gap <- 0.04 / cell
  ordered <- labs[ord]

  ## Row k of the ordering is drawn at y = n - k + 1, so the first is on top.
  cells <- data.frame(x = rep(seq_len(n), each = n), y = rep(n:1, times = n),
    value = as.vector(shown))
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(ggplot2::aes(x = .data$x, y = .data$y,
        fill = .data$value), data = cells,
      colour = if (n <= 30) "white" else NA, linewidth = 0.3) +
    ggplot2::scale_fill_gradientn(colours = colours, limits = limits,
      oob = function(x, range) pmin(pmax(x, range[1]), range[2]),
      na.value = "#E6E6E6", name = PLOT_COMPARISON_LABELS[[value]],
      guide = ggplot2::guide_colourbar(barheight = grid::unit(
        min(n * cell, 2.5), "in"), barwidth = grid::unit(0.12, "in")))

  if (annotate) {
    ann <- cells[is.finite(cells$value), , drop = FALSE]
    scaled <- (pmin(pmax(ann$value, limits[1]), limits[2]) - limits[1]) /
      (limits[2] - limits[1])
    rgb <- grDevices::colorRamp(colours)(scaled) / 255
    dark <- 0.2126 * rgb[, 1] + 0.7152 * rgb[, 2] + 0.0722 * rgb[, 3] < 0.5
    ann$label <- formatC(ann$value, format = "f",
      digits = if (value == "score") 2 else 1)
    ann$label[ann$label == "-0.0"] <- "0.0"
    ann$colour <- ifelse(dark, "white", "black")
    p <- p + ggplot2::geom_text(ggplot2::aes(x = .data$x, y = .data$y,
        label = .data$label, colour = .data$colour), data = ann,
      size = min(0.8 * text.size, 0.3 * cell * 72) / ggplot2::.pt) +
      ggplot2::scale_colour_identity()
  }

  xmin <- 0.5
  xmax <- n + 0.5
  ymin <- 0.5
  ymax <- n + 0.5

  if (!is.null(cut.at)) {
    n.groups <- length(unique(groups))
    unit <- c(score = "score", pvalue = "P-value", qvalue = "q-value")[[cut.unit]]
    what <- if (!is.null(cut)) paste(unit, format(cut)) else
      paste("height", format(cut.at))
    cut.label <- paste0(what, ", ", n.groups, " group",
      if (n.groups != 1) "s" else "")
    ## Groups take consecutive rows, since the order is the tree's.
    runs <- rle(unname(groups[ord]))
    ends <- cumsum(runs$lengths)
    starts <- ends - runs$lengths + 1
    keep <- runs$lengths > 1
    if (any(keep)) {
      boxes <- data.frame(xmin = starts[keep] - 0.5, xmax = ends[keep] + 0.5,
        ymin = n - ends[keep] + 0.5, ymax = n - starts[keep] + 1.5)
      p <- p + ggplot2::geom_rect(ggplot2::aes(xmin = .data$xmin,
          xmax = .data$xmax, ymin = .data$ymin, ymax = .data$ymax),
        data = boxes, fill = NA, colour = PLOT_COMPARISON_CUT, linewidth = 0.6)
    }
  }

  if (left.tree || top.tree) {
    segs <- dendrogram_segments(tree, ord)
    top <- max(segs$top, if (is.null(cut.at)) 0 else cut.at) * 1.02
    if (left.tree) {
      ## Heights grow leftwards from the heatmap.
      hx <- function(h) 0.5 - gap - h / top * tree.w
      p <- p + ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$xend,
          y = .data$y, yend = .data$yend),
        data = data.frame(x = hx(segs$h), xend = hx(segs$hend),
          y = n + 1 - segs$pos, yend = n + 1 - segs$posend),
        linewidth = 0.3, lineend = "square")
      xmin <- 0.5 - gap - tree.w
      if (!is.null(cut.at)) {
        p <- p + ggplot2::annotate("segment", x = hx(cut.at), xend = hx(cut.at),
            y = 0.5, yend = n + 0.5, colour = PLOT_COMPARISON_CUT,
            linetype = "dashed", linewidth = 0.4) +
          ggplot2::annotate("text", x = hx(cut.at), y = n + 0.5 + 0.04 / cell,
            label = cut.label, vjust = 0, size = 0.75 * text.size / ggplot2::.pt,
            colour = PLOT_COMPARISON_CUT)
        ymax <- max(ymax, n + 0.5 + (0.08 + 0.75 * text.size / 72) / cell)
        ## Text is sized in points, so its room in cells is an estimate.
        xmin <- min(xmin, hx(cut.at) - nchar(cut.label) * 0.3 * 0.75 * text.size /
          72 / cell)
      }
    }
    if (top.tree) {
      hy <- function(h) n + 0.5 + gap + h / top * tree.w
      p <- p + ggplot2::geom_segment(ggplot2::aes(x = .data$x, xend = .data$xend,
          y = .data$y, yend = .data$yend),
        data = data.frame(x = segs$pos, xend = segs$posend, y = hy(segs$h),
          yend = hy(segs$hend)),
        linewidth = 0.3, lineend = "square")
      ymax <- max(ymax, n + 0.5 + gap + tree.w)
      if (!is.null(cut.at)) {
        p <- p + ggplot2::annotate("segment", x = 0.5, xend = n + 0.5,
          y = hy(cut.at), yend = hy(cut.at), colour = PLOT_COMPARISON_CUT,
          linetype = "dashed", linewidth = 0.4)
        if (!left.tree) {
          p <- p + ggplot2::annotate("text", x = n + 0.5 + 0.06 / cell,
            y = hy(cut.at), label = cut.label, hjust = 0,
            size = 0.75 * text.size / ggplot2::.pt, colour = PLOT_COMPARISON_CUT)
        }
      }
    }
  }

  ## The room taken right of the heatmap by a cut label over the top tree.
  right <- n + 0.5
  if (!is.null(cut.at) && top.tree && !left.tree)
    right <- n + 0.5 + (0.06 + nchar(cut.label) * 0.6 * 0.75 * text.size / 72) / cell
  if (logos) {
    widths <- vapply(motifs, function(m) ncol(m@motif), integer(1))
    use.type <- if (all(vapply(motifs, function(m) m@type, "") == "CWM")) "CWM"
      else "ICM"
    keys <- paste0("motif", seq_len(n))
    logo.list <- stats::setNames(motifs, keys)
    x0 <- n + 0.5 + 0.08 / cell
    p <- p + geom_motif(ggplot2::aes(x = .data$x, y = .data$y,
        motif = .data$motif),
      data = data.frame(x = x0, y = n:1, motif = keys[ord]),
      motif = logo.list, use.type = use.type, width = 0.3, height = 0.8)
    right <- max(right, x0 + 0.3 * max(widths))
  }
  if (show.labels) {
    x0 <- (if (logos) x0 + 0.3 * max(widths) + 0.04 / cell else n + 0.5) + 0.06 / cell
    p <- p + ggplot2::annotate("text", x = x0, y = n:1, label = ordered,
        hjust = 0, size = label.size / ggplot2::.pt) +
      ggplot2::annotate("text", x = seq_len(n), y = 0.5 - 0.03 / cell,
        label = ordered, angle = 90, hjust = 1,
        size = label.size / ggplot2::.pt)
    ## Text is sized in points, so its room in cells is an estimate.
    label.w <- max(nchar(ordered)) * 0.6 * label.size / 72 / cell
    right <- max(right, x0 + label.w)
    ymin <- 0.5 - 0.03 / cell - label.w
  }

  p <- p + ggplot2::coord_fixed(xlim = c(xmin, right), ylim = c(ymin, ymax),
      expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_size = 0.8 * text.size) +
    ggplot2::theme(plot.margin = grid::unit(rep(0.1, 4), "in"),
      legend.justification = "top")
  ## The size in inches the plot is drawn for, with the colour bar beside it.
  legend.w <- 0.25 + max(0.8, nchar(PLOT_COMPARISON_LABELS[[value]]) * 0.55 *
    0.8 * text.size / 72)
  attr(p, "size") <- c(width = (right - xmin) * cell + legend.w + 0.2,
    height = max((ymax - ymin) * cell, min(n * cell, 2.5) + 0.5) + 0.2)
  p

}

PLOT_COMPARISON_BLUES <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1",
  "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
PLOT_COMPARISON_CUT <- "#C0392B"
PLOT_COMPARISON_LABELS <- list(score = "Similarity score",
  pvalue = "-log10(P-value)", qvalue = "-log10(q-value)")

## The motifs, the matrix of values as given or computed, the score matrix
## to cluster on when there is one, and the labels.
comparison_input <- function(motifs, matrix, labels, value, min.overlap, RC,
  nthreads) {

  if (!is.null(motifs)) {
    motifs <- convert_motifs(motifs)
    if (!is.list(motifs)) motifs <- list(motifs)
  }
  if (!is.null(matrix)) {
    vals <- as.matrix(matrix)
    storage.mode(vals) <- "double"
    if (nrow(vals) != ncol(vals))
      stop("`matrix` must be square", call. = FALSE)
    if (!is.null(motifs) && length(motifs) != nrow(vals))
      stop("`matrix` must have one row and column per motif", call. = FALSE)
    scores <- if (value == "score") vals else NULL
  } else {
    if (is.null(motifs))
      stop("give `motifs`, a `matrix` or both", call. = FALSE)
    if (length(motifs) < 2)
      stop("at least two motifs are needed for a comparison", call. = FALSE)
    scores <- compare_motifs_lite(motifs, min.overlap = min.overlap, RC = RC,
      matrix.out = "score", nthreads = nthreads)
    vals <- if (value == "score") scores else
      compare_motifs_lite(motifs, min.overlap = min.overlap, RC = RC,
        matrix.out = value, nthreads = nthreads)
    scores <- unname(as.matrix(scores))
    vals <- unname(as.matrix(vals))
  }
  n <- nrow(vals)
  if (n < 2)
    stop("at least two motifs are needed for a comparison", call. = FALSE)

  if (is.character(labels) && length(labels) == 1 && n != 1 &&
      (is.null(motifs) || labels %in% methods::slotNames(motifs[[1]]))) {
    if (!is.null(motifs)) {
      labs <- vapply(motifs, function(m) {
        x <- methods::slot(m, labels)
        if (length(x) == 0 || is.na(x[1])) "" else as.character(x[1])
      }, character(1))
    } else {
      labs <- colnames(vals)
      if (is.null(labs)) labs <- as.character(seq_len(n))
    }
    if (labels == "name" || is.null(motifs)) labs <- make.unique(labs)
  } else {
    labs <- as.character(labels)
    if (length(labs) != n)
      stop("give one label per motif", call. = FALSE)
  }
  list(motifs = motifs, vals = unname(vals), scores = scores, labels = labs)

}

## Distances to cluster on: (1 - score) / 2 for scores, and otherwise one
## minus the mean -log10 value of both directions over its largest.
comparison_distances <- function(vals, scores) {

  if (!is.null(scores)) {
    dist <- (1 - (scores + t(scores)) / 2) / 2
  } else {
    logs <- neg_log10(vals)
    sym <- (logs + t(logs)) / 2
    diag(sym) <- NA
    top <- comparison_log_top(sym)
    dist <- if (top > 0) pmin(pmax(1 - sym / top, 0), 1) else sym * 0 + 1
  }
  dist[!is.finite(dist)] <- 1
  diag(dist) <- 0
  dist

}

comparison_log_top <- function(sym) {
  finite <- sym[is.finite(sym)]
  if (length(finite)) max(finite) else 0
}

## The tree height of a cutoff given in the units clustered on, the way
## comparison_distances() turns values into distances.
comparison_cut_height <- function(cut, vals, unit) {

  if (length(cut) != 1 || !is.finite(cut))
    stop("`cut` must be a single number", call. = FALSE)
  if (unit == "score") return((1 - cut) / 2)
  if (cut <= 0 || cut > 1)
    stop("a P- or q-value `cut` must be above 0 and at most 1", call. = FALSE)
  logs <- neg_log10(vals)
  sym <- (logs + t(logs)) / 2
  diag(sym) <- NA
  top <- comparison_log_top(sym)
  if (top > 0) min(max(1 - neg_log10(cut) / top, 0), 1) else 1

}

neg_log10 <- function(p) 0 - log10(pmax(p, .Machine$double.xmin))

## Dendrogram segments of an hclust tree as positions along the heatmap and
## heights, with the leaves at height 0, and the largest height.
dendrogram_segments <- function(tree, ord) {

  n <- length(ord)
  node.pos <- numeric(n)
  node.pos[ord] <- seq_len(n)
  merge.pos <- numeric(n - 1)
  merge.h <- numeric(n - 1)
  pos.of <- function(i) if (i < 0) node.pos[-i] else merge.pos[i]
  h.of <- function(i) if (i < 0) 0 else merge.h[i]
  segs <- vector("list", n - 1)
  for (k in seq_len(n - 1)) {
    a <- tree$merge[k, 1]
    b <- tree$merge[k, 2]
    h <- tree$height[k]
    pa <- pos.of(a)
    pb <- pos.of(b)
    ## A bracket: up from each child to the merge height, and across.
    segs[[k]] <- data.frame(pos = c(pa, pa, pb), posend = c(pa, pb, pb),
      h = c(h.of(a), h, h.of(b)), hend = c(h, h, h))
    merge.pos[k] <- (pa + pb) / 2
    merge.h[k] <- h
  }
  segs <- do.call(rbind, segs)
  top <- max(c(merge.h[is.finite(merge.h)], 0))
  segs$top <- if (top > 0) top else 1
  segs

}
