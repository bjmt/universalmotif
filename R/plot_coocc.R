#' Plot motif co-occurrence results as a heatmap.
#'
#' Draw the results of [motif_coocc()] as a heatmap with a cell for every
#' pair of motifs, filled the same on both sides of the diagonal with the
#' log2 odds ratio of the pair (on a scale centred on zero) or with its
#' `-log10` q-value or P-value. Pairs with q-values at or below `cutoff`
#' are marked with a dot.
#'
#' @param results `data.frame` returned by [motif_coocc()], with the
#'   `motif_a`, `motif_b`, `odds_ratio`, `pvalue` and `qvalue` columns.
#' @param value `character(1)` One of `c("odds_ratio", "qvalue", "pvalue")`,
#'   the value the cells are filled with.
#' @param cutoff `numeric(1)` The q-value at or below which pairs are marked,
#'   in `(0, 1]`.
#' @param order `character` `"cluster"` to order the motifs by
#'   average-linkage clustering of their rows of shown values, `"input"` to
#'   keep them in the order they first appear in `results`, or the names of
#'   the motifs in the order to draw them, which can leave some out.
#' @param colours `character` or `NULL` The colours of the fill gradient,
#'   from lowest to highest. By default a diverging blue to red scale for
#'   odds ratios and a sequential blue scale for q-values and P-values.
#' @param text.size `numeric(1)` The font size of the motif names and the
#'   legends, in points.
#'
#' @return A `ggplot` object.
#'
#' @details
#' Untested pairs (those with `NA` P-values, which [motif_coocc()] gives to
#' pairs that co-occur in fewer than `min.coocc` sequences and to self
#' pairs) and the diagonal are filled grey, as are tested pairs whose odds
#' ratio is undefined. Odds ratios of zero or infinity are drawn at the ends
#' of the colour scale, whose limits are set by the largest finite log2 odds
#' ratio.
#'
#' The dots are drawn in white on dark cells and in black on light ones.
#' When clustering, untested pairs count as having no association (a log2
#' odds ratio or `-log10` P-value of zero), infinite log2 odds ratios as the
#' largest finite one, and a missing diagonal as the strongest association,
#' so that motifs found together end up next to each other.
#'
#' @examples
#' motifs <- list(create_motif("TTGACCAA", name = "W-box", nsites = 20),
#'   create_motif("CACGTGGC", name = "G-box", nsites = 20),
#'   create_motif("TATAAATA", name = "TATA-box", nsites = 20),
#'   create_motif("CCGGAAGT", name = "ETS", nsites = 20))
#' seqs <- create_sequences(seqnum = 200, seqlen = 200, rng.seed = 1)
#' ## Plant W-box and G-box sites together, and TATA-box sites alone.
#' set.seed(1)
#' seqs[1:60] <- implant_motifs(motifs[[1]], seqs[1:60])
#' seqs[1:60] <- implant_motifs(motifs[[2]], seqs[1:60])
#' seqs[61:120] <- implant_motifs(motifs[[3]], seqs[61:120])
#' res <- motif_coocc(motifs, seqs, pvalue = 1e-4, nthreads = 1)
#' plot_coocc(res)
#'
#' ## Fill by q-value, with the motifs in a given order.
#' plot_coocc(res, value = "qvalue", order = c("W-box", "G-box", "TATA-box"))
#'
#' @seealso [motif_coocc()], [plot_enrichment()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_coocc <- function(results, value = c("odds_ratio", "qvalue", "pvalue"),
  cutoff = 0.1, order = "cluster", colours = NULL, text.size = 9) {

  value <- match.arg(value)
  if (methods::is(results, "DataFrame")) results <- as.data.frame(results)
  if (!is.data.frame(results))
    stop("`results` must be a data.frame from motif_coocc()", call. = FALSE)
  for (col in c("motif_a", "motif_b", "odds_ratio", "pvalue", "qvalue")) {
    if (!col %in% names(results))
      stop("`results` has no `", col, "` column; was it made by ",
        "motif_coocc()?", call. = FALSE)
  }
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) ||
      cutoff <= 0 || cutoff > 1)
    stop("`cutoff` must be a single number in (0, 1]", call. = FALSE)
  if (!is.character(order) || !length(order) || anyNA(order))
    stop("`order` must be \"cluster\", \"input\" or motif names",
      call. = FALSE)
  if (!is.null(colours) && (!is.character(colours) || length(colours) < 2L))
    stop("`colours` must be NULL or at least two colours", call. = FALSE)
  if (!is.numeric(text.size) || length(text.size) != 1L || is.na(text.size) ||
      text.size <= 0)
    stop("`text.size` must be a single positive number", call. = FALSE)
  if (!nrow(results))
    stop("`results` has no rows to plot", call. = FALSE)

  mats <- coocc_matrix(results, value)
  motifs <- mats$motifs
  vals <- mats$vals
  qs <- mats$qs

  if (length(order) == 1L && order %in% c("cluster", "input")) {
    idx <- seq_along(motifs)
    if (order == "cluster" && length(motifs) > 2L)
      idx <- coocc_cluster_order(vals)
  } else {
    wanted <- unique(order)
    missing <- setdiff(wanted, motifs)
    if (length(missing))
      stop("`results` has no motifs named ",
        paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
    idx <- match(wanted, motifs)
  }
  motifs <- motifs[idx]
  vals <- vals[idx, idx, drop = FALSE]
  qs <- qs[idx, idx, drop = FALSE]
  n <- length(motifs)

  finite <- is.finite(vals)
  if (value == "odds_ratio") {
    lim <- if (any(finite)) max(abs(vals[finite])) else 1
    if (lim == 0) lim <- 1
    limits <- c(-lim, lim)
    if (is.null(colours))
      colours <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
        "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
    title <- "log2(odds ratio)"
  } else {
    top <- if (any(finite)) max(vals[finite]) else 1
    if (top == 0) top <- 1
    limits <- c(0, top)
    if (is.null(colours))
      colours <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
        "#4292C6", "#2171B5", "#08519C", "#08306B")
    title <- if (value == "qvalue") "-log10(q-value)" else "-log10(P-value)"
  }
  shown <- pmin(pmax(vals, limits[1]), limits[2])

  cells <- data.frame(
    a = factor(rep(motifs, times = n), levels = motifs),
    b = factor(rep(motifs, each = n), levels = rev(motifs)),
    value = as.vector(shown), qvalue = as.vector(qs))
  marked <- cells[!is.na(cells$qvalue) & cells$qvalue <= cutoff, ,
    drop = FALSE]
  marked$colour <- coocc_dot_colours(marked$value, limits, colours)
  marked$mark <- rep("mark", nrow(marked))

  ## Dots shrink as the motifs grow in number, in square points.
  dot <- sqrt(min(16, 900 / n^1.2)) / ggplot2::.pt
  cut.label <- formatC(cutoff, format = "g")

  ggplot2::ggplot(cells, ggplot2::aes(x = .data$a, y = .data$b)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value), colour = "white",
      linewidth = 0.8) +
    ggplot2::geom_point(data = marked,
      ggplot2::aes(colour = .data$colour, shape = .data$mark),
      size = dot, stroke = 0) +
    ggplot2::scale_fill_gradientn(colours = colours, limits = limits,
      na.value = "#D0D0D0", name = title,
      guide = ggplot2::guide_colourbar(order = 2L,
        theme = ggplot2::theme(legend.key.height =
          grid::unit(min(5 + 0.6 * n, 16), "lines")))) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_shape_manual(values = c(mark = 16), limits = "mark",
      labels = as.expression(bquote(q <= .(cut.label))), name = NULL,
      guide = ggplot2::guide_legend(order = 1L,
        override.aes = list(colour = "#1A1A1A", size = 4 / ggplot2::.pt))) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = text.size, colour = "black"),
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
        vjust = 0.5),
      legend.title = ggplot2::element_text(size = text.size),
      legend.text = ggplot2::element_text(size = text.size)
    )

}

#-----------------------------------------------------------------------------
# Internal helpers

## The motifs, in order of first appearance, with the shown values and
## q-values of each pair in symmetric matrices, NA where untested.
coocc_matrix <- function(results, value) {
  a <- as.character(results$motif_a)
  b <- as.character(results$motif_b)
  a[is.na(a)] <- "NA"
  b[is.na(b)] <- "NA"
  motifs <- unique(c(a, b))
  n <- length(motifs)
  shown <- if (value == "odds_ratio") log2(as.numeric(results$odds_ratio))
    else -log10(pmax(as.numeric(results[[value]]), 5e-324))
  q <- as.numeric(results$qvalue)
  tested <- !is.na(results$pvalue)
  vals <- qs <- matrix(NA_real_, n, n)
  i <- match(a, motifs)[tested]
  j <- match(b, motifs)[tested]
  vals[cbind(i, j)] <- vals[cbind(j, i)] <- shown[tested]
  qs[cbind(i, j)] <- qs[cbind(j, i)] <- q[tested]
  ## NaN odds ratios count as missing.
  vals[is.nan(vals)] <- NA_real_
  list(motifs = motifs, vals = vals, qs = qs)
}

## Untested pairs count as zero, which is no association for log odds ratios
## and a P-value of 1 otherwise, infinite log odds as the largest finite one,
## and a missing diagonal as the strongest association, so that two motifs
## found together end up close.
coocc_cluster_order <- function(vals) {
  finite <- is.finite(vals)
  big <- if (any(finite)) max(abs(vals[finite])) else 1
  rows <- vals
  rows[is.na(rows)] <- 0
  rows <- pmin(pmax(rows, -big), big)
  diagonal <- is.na(diag(vals))
  diag(rows)[diagonal] <- big
  stats::hclust(stats::dist(rows), method = "average")$order
}

## Dark dots on light cells and white dots on dark ones, with the cell
## colours worked out from the same gradient as the fill scale.
coocc_dot_colours <- function(values, limits, colours) {
  if (!length(values)) return(character(0))
  ramp <- grDevices::colorRamp(colours, space = "Lab")
  t <- (values - limits[1]) / (limits[2] - limits[1])
  rgb <- matrix(as.vector(grDevices::col2rgb("#D0D0D0")), length(values), 3L,
    byrow = TRUE)
  ok <- !is.na(t)
  if (any(ok)) rgb[ok, ] <- ramp(t[ok])
  light <- drop(rgb %*% c(0.299, 0.587, 0.114)) / 255 > 0.55
  ifelse(light, "#1A1A1A", "#FFFFFF")
}
