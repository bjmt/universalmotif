#' Plot motif enrichment results as a volcano plot.
#'
#' Draw the results of [enrich_motifs_lite()] or [enrich_motifs()] with each
#' motif as a point at its log2 enrichment and `-log10` q-value (or P-value
#' with `y = "pvalue"`), sized by the number of target sequences with hits.
#' Motifs at or below `cutoff` are coloured, a dashed line marks the cutoff,
#' and the most significant motifs are labelled, optionally with their logos
#' beside the labels.
#'
#' @param results `data.frame` or `DataFrame` returned by
#'   [enrich_motifs_lite()] or [enrich_motifs()]. It needs a `motif` column,
#'   either a `log2.enrichment` or a `target.enrichment` column, and the
#'   P-value or q-value column chosen with `y` (`pvalue`/`qvalue` from
#'   [enrich_motifs_lite()], `Pval`/`Qval` from [enrich_motifs()]).
#' @param y `character(1)` One of `c("qvalue", "pvalue")`, the value drawn on
#'   the vertical axis and compared to `cutoff`.
#' @param cutoff `numeric(1)` The q-value or P-value at or below which motifs
#'   are significant, in `(0, 1]`.
#' @param labels `numeric(1)` or `character` The number of most significant
#'   motifs to label, or the names of the motifs to label whether they are
#'   significant or not. Labels are moved apart so that they do not cover
#'   each other or other points.
#' @param size.by `character(1)` or `NULL` A column of `results` to size the
#'   points by, or `NULL` for points of one size. The names of the hit count
#'   columns of [enrich_motifs_lite()] and [enrich_motifs()] can be used for
#'   either table (for example `"target.site.hits"` also finds the
#'   `target.hits` column of [enrich_motifs()]).
#' @param motifs See [convert_motifs()] for accepted motif formats. The
#'   motifs given to the enrichment function, to draw a logo beside each
#'   label. Motifs are matched to the rows of `results` by the names of the
#'   list when it has them, and otherwise by motif name (with duplicates
#'   renamed by [make.unique()], as [enrich_motifs_lite()] does).
#' @param logo.height `numeric(1)` or `NULL` The height of the logos in
#'   points, by default twice `text.size`.
#' @param colour `character(1)` The colour of significant points.
#' @param text.size `numeric(1)` The font size of the labels and legend, in
#'   points.
#'
#' @return A `ggplot` object.
#'
#' @details
#' Motifs without any target hits have an enrichment of zero, and so an
#' infinitely negative log2 enrichment; these are drawn as triangles at the
#' left edge of the plot instead. Likewise motifs with hits in the targets
#' but none in the background ([enrich_motifs()] can give these an infinite
#' enrichment) are drawn as triangles at the right edge. P-values and
#' q-values too small to be stored as doubles are drawn at the smallest
#' representable value.
#'
#' The labels are placed when the plot is drawn, so that they stay next to
#' their points and apart from each other whatever size the plot is drawn
#' at. Each label in turn, from the most significant, tries spots at growing
#' distances around its point and takes the one that covers the fewest of
#' the labels placed before it, the other points and the edge of the
#' panel. Labels that end up away from their points are joined to them by
#' a thin line.
#'
#' [enrich_motifs_lite()] only returns motifs with q-values of 0.1 or less
#' by default; use `qvalue = 1` to plot every motif, as the plot is
#' easier to read with the motifs that fall short of the cutoff in it too.
#'
#' @examples
#' data(ArabidopsisPromoters)
#' data(ArabidopsisMotif)
#' motifs <- list(ArabidopsisMotif,
#'   create_motif("TTGACY", name = "W-box", nsites = 20),
#'   create_motif("CACGTG", name = "G-box", nsites = 20),
#'   create_motif("TATAWAWA", name = "TATA-box", nsites = 20),
#'   create_motif("CCGGAAGT", name = "ETS", nsites = 20))
#' bkg <- shuffle_sequences(ArabidopsisPromoters, k = 2, rng.seed = 1)
#' res <- enrich_motifs_lite(motifs, ArabidopsisPromoters, bkg,
#'   pvalue = 1e-3, qvalue = 1, nthreads = 1)
#' plot_enrichment(res, labels = 3, motifs = motifs)
#'
#' ## Label motifs by name, and draw P-values with points of one size.
#' plot_enrichment(res, y = "pvalue", cutoff = 0.001,
#'   labels = c("W-box", "ETS"), size.by = NULL)
#'
#' @seealso [enrich_motifs_lite()], [enrich_motifs()], [geom_motif()],
#'   [plot_coocc()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_enrichment <- function(results, y = c("qvalue", "pvalue"), cutoff = 0.1,
  labels = 10, size.by = "target.seq.hits", motifs = NULL, logo.height = NULL,
  colour = "#255C99", text.size = 9) {

  y <- match.arg(y)
  if (methods::is(results, "DataFrame")) results <- as.data.frame(results)
  if (!is.data.frame(results))
    stop("`results` must be a data.frame from enrich_motifs_lite() or ",
      "enrich_motifs()", call. = FALSE)
  if (!is.numeric(cutoff) || length(cutoff) != 1L || is.na(cutoff) ||
      cutoff <= 0 || cutoff > 1)
    stop("`cutoff` must be a single number in (0, 1]", call. = FALSE)
  if (!is.numeric(text.size) || length(text.size) != 1L || is.na(text.size) ||
      text.size <= 1)
    stop("`text.size` must be a single number greater than 1", call. = FALSE)
  if (!is.null(logo.height) && (!is.numeric(logo.height) ||
      length(logo.height) != 1L || is.na(logo.height) || logo.height <= 0))
    stop("`logo.height` must be NULL or a single positive number",
      call. = FALSE)
  if (is.numeric(labels)) {
    if (length(labels) != 1L || is.na(labels) || labels < 0)
      stop("`labels` must be a non-negative number or motif names",
        call. = FALSE)
  } else if (!is.character(labels)) {
    stop("`labels` must be a non-negative number or motif names",
      call. = FALSE)
  }

  cols <- names(results)
  if (!"motif" %in% cols ||
      !any(c("log2.enrichment", "target.enrichment") %in% cols))
    stop("`results` needs `motif` and `log2.enrichment` (or ",
      "`target.enrichment`) columns; was it made by enrich_motifs_lite() ",
      "or enrich_motifs()?", call. = FALSE)
  ycol <- intersect(if (y == "qvalue") c("qvalue", "Qval") else
    c("pvalue", "Pval"), cols)
  if (!length(ycol))
    stop("`results` has no `", y, "` column", call. = FALSE)
  if (!is.null(size.by)) {
    if (!is.character(size.by) || length(size.by) != 1L)
      stop("`size.by` must be NULL or a single column name", call. = FALSE)
    size.by <- enrich_size_column(size.by, cols)
  }
  if (!nrow(results))
    stop("`results` has no rows to plot", call. = FALSE)

  nm <- as.character(results$motif)
  nm[is.na(nm)] <- "NA"
  x <- if ("log2.enrichment" %in% cols) as.numeric(results$log2.enrichment)
    else log2(as.numeric(results$target.enrichment))
  yv <- -log10(pmax(as.numeric(results[[ycol[1]]]), 5e-324))
  sv <- if (is.null(size.by)) rep(0, length(nm))
    else as.numeric(results[[size.by]])
  sv[is.na(sv)] <- 0
  keep <- !is.na(x) & !is.na(yv)
  if (!any(keep))
    stop("`results` has no enrichment values to plot", call. = FALSE)
  line <- -log10(cutoff)
  sig <- keep & yv >= line

  if (is.numeric(labels)) {
    chosen <- which(sig)
    chosen <- chosen[order(-yv[chosen], -x[chosen])]
    chosen <- utils::head(chosen, floor(labels))
  } else {
    missing <- setdiff(labels, nm)
    if (length(missing))
      stop("`results` has no motifs named ",
        paste(utils::head(missing, 5L), collapse = ", "), call. = FALSE)
    chosen <- unique(match(labels, nm))
    chosen <- chosen[keep[chosen]]
    chosen <- chosen[order(-yv[chosen], -x[chosen])]
  }

  logos <- list()
  if (!is.null(motifs) && length(chosen))
    logos <- enrich_label_logos(motifs, nm[chosen])

  ## Motifs with no target hits (or no background hits) have no finite
  ## enrichment and sit at the edges instead.
  edge.lo <- keep & x == -Inf
  edge.hi <- keep & x == Inf
  finite <- keep & is.finite(x)
  xlo <- if (any(finite)) min(x[finite], 0) else -1
  xhi <- if (any(finite)) max(x[finite], 0) else 1
  pad <- 0.08 * max(xhi - xlo, 1)
  left <- xlo - (if (any(edge.lo)) 2 else 1) * pad
  right <- xhi + (if (any(edge.hi)) 2 else 1) * pad
  x[edge.lo] <- left + 0.35 * pad
  x[edge.hi] <- right - 0.35 * pad
  top <- 1.12 * max(yv[keep], line)
  if (top == 0) top <- 1

  rank <- rep(NA_integer_, length(nm))
  rank[chosen] <- seq_along(chosen)
  d <- data.frame(motif = nm, x = x, y = yv, size = sv,
    fill = ifelse(sig, colour, "#B4B4B4"), significant = sig,
    edge = ifelse(edge.lo, -1, ifelse(edge.hi, 1, 0)),
    label = ifelse(is.na(rank), NA_character_, nm), rank = rank,
    stringsAsFactors = FALSE)[keep, , drop = FALSE]
  ## Significant points are drawn over the others.
  d <- d[order(d$significant), , drop = FALSE]

  lo <- min(d$size)
  hi <- max(d$size)
  stroke <- 0.4
  sized <- !is.null(size.by) && hi > lo
  fixed <- if (sized) list() else
    list(size = enrich_area_to_size(36, stroke))
  ## The point layers map size only when there is a size scale.
  with.size <- function(mapping) {
    if (sized) mapping$size <- ggplot2::aes(size = .data$size)$size
    mapping
  }

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_vline(xintercept = 0, colour = "#E2E2E2",
      linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = line, colour = "#7F7F7F",
      linewidth = 0.3, linetype = "dashed") +
    ggplot2::annotate("text", x = Inf, y = line, hjust = 1.02, vjust = -0.6,
      label = paste(if (y == "qvalue") "q =" else "p =",
        formatC(cutoff, format = "g")),
      size = (text.size - 1) / ggplot2::.pt, colour = "#5E5E5E") +
    do.call(ggplot2::geom_point, c(list(
      mapping = with.size(ggplot2::aes(fill = .data$fill)),
      data = d[d$edge == 0, , drop = FALSE], shape = 21,
      colour = "white", stroke = stroke), fixed))

  if (any(d$edge != 0))
    p <- p + do.call(layer_enrich_edges, c(list(
      mapping = with.size(ggplot2::aes(fill = .data$fill,
        direction = .data$edge)),
      data = d[d$edge != 0, , drop = FALSE], stroke = stroke), fixed))

  if (length(chosen))
    p <- p + do.call(layer_enrich_labels, c(list(
      mapping = with.size(ggplot2::aes(label = .data$label,
        rank = .data$rank)),
      data = d, logos = logos, text.size = text.size,
      logo.height = if (is.null(logo.height)) 2 * text.size else logo.height,
      stroke = stroke), fixed))

  p <- p +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_cartesian(xlim = c(left, right), ylim = c(0, top),
      expand = FALSE, clip = "off") +
    ggplot2::labs(x = "log2(enrichment)",
      y = if (y == "qvalue") "-log10(q-value)" else "-log10(P-value)") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(colour = "black"),
      axis.line.y.left = ggplot2::element_line(colour = "black"),
      legend.title = ggplot2::element_text(size = text.size),
      legend.text = ggplot2::element_text(size = text.size),
      legend.justification = c(0, 1),
      legend.key = ggplot2::element_blank()
    )

  if (sized) {
    breaks <- enrich_size_breaks(lo, hi)
    p <- p + ggplot2::continuous_scale("size",
      palette = function(t) enrich_area_to_size(12 + 108 * t, stroke),
      name = enrich_size_title(size.by), breaks = breaks,
      labels = format(breaks, trim = TRUE), limits = c(lo, hi),
      guide = ggplot2::guide_legend(override.aes = list(fill = "#B4B4B4")))
  }

  p

}

#-----------------------------------------------------------------------------
# Internal helpers

## Find the column to size points by, accepting the hit count column names of
## either enrichment function.
enrich_size_column <- function(size.by, cols) {
  if (size.by %in% cols) return(size.by)
  alias <- c(target.site.hits = "target.hits", target.hits = "target.site.hits",
    bkg.site.hits = "bkg.hits", bkg.hits = "bkg.site.hits")
  if (size.by %in% names(alias) && alias[[size.by]] %in% cols)
    return(alias[[size.by]])
  stop("`results` has no `", size.by, "` column to size points by",
    call. = FALSE)
}

enrich_size_title <- function(size.by) {
  switch(size.by,
    target.seq.hits = "Target sequences\nwith hits",
    target.site.hits = , target.hits = "Target hits",
    bkg.seq.hits = "Background sequences\nwith hits",
    bkg.site.hits = , bkg.hits = "Background hits",
    size.by)
}

## Point areas are given in square points, as for matplotlib, and turned
## into ggplot2 sizes (the diameter in mm, less the border).
enrich_area_to_size <- function(area, stroke) {
  pmax(sqrt(area) - stroke * ggplot2::.stroke / 2, 0.5) / ggplot2::.pt
}

## Two or three round numbers within the range of the sizes, whole numbers
## when the sizes are.
enrich_size_breaks <- function(lo, hi) {
  mag <- 10^(floor(log10(hi - lo)) + (-2):1)
  steps <- sort(as.vector(outer(c(1, 2, 2.5, 4, 5), mag)))
  if (lo == round(lo) && hi == round(hi)) steps <- steps[steps == round(steps)]
  for (s in steps) {
    b <- seq(ceiling(lo / s), floor(hi / s)) * s
    if (length(b) <= 3L) return(if (length(b) >= 2L) b else unique(c(lo, hi)))
  }
  unique(c(lo, hi))
}

## Logo polygons for the labelled motifs, rescaled into a unit box with the
## full height of the box standing for the most information a position can
## hold, so that weak positions stay small.
enrich_label_logos <- function(motifs, wanted) {
  given <- names(motifs)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  keys <- make.unique(vapply(motifs, function(m) m@name, character(1)))
  if (!is.null(given) && length(given) == length(motifs)) {
    named <- !is.na(given) & nzchar(given)
    idx <- match(wanted, given[named])
    idx <- which(named)[idx]
  } else {
    idx <- rep(NA_integer_, length(wanted))
  }
  idx[is.na(idx)] <- match(wanted[is.na(idx)], keys)
  if (anyNA(idx))
    stop("`motifs` has no motif named ",
      paste(utils::head(wanted[is.na(idx)], 5L), collapse = ", "),
      call. = FALSE)
  fontDF <- as.data.frame(fontDFroboto)
  out <- lapply(idx, function(i) {
    info <- motif_to_polygons(motifs[[i]], "ICM", fontDF, NULL, 0.01, 0.04,
      0.01, TRUE, TRUE, NULL, FALSE)
    poly <- info$poly
    top <- log2(nrow(motifs[[i]]@motif))
    list(ncol = info$ncol,
      x = (poly$x - info$xlo) / (info$xhi - info$xlo),
      y = pmin(poly$y / top, 1),
      id = match(poly$letter.id, unique(poly$letter.id)),
      fill = poly$fill_colour[!duplicated(poly$letter.id)])
  })
  names(out) <- wanted
  out
}

## Greedily place label boxes of the given sizes next to points. Everything
## is in points from the bottom left of the panel. Each label in `which`,
## taken in turn, tries spots on rings of growing distance around its point
## and takes the one that overlaps the labels placed before it least, then
## other points, then the edge of `frame`, preferring closer spots. Returns
## the bottom left corner of each box relative to its point, and whether it
## sits beyond the nearest ring, where it could be taken for a neighbour's
## label and so needs a leader line.
enrich_place_labels <- function(px, py, radius, w, h, which, frame) {
  placed <- matrix(numeric(0), ncol = 4L)
  dots <- cbind(px - radius, py - radius, px + radius, py + radius)
  gaps <- c(3, 9, 17, 27, 40, 56, 76, 100, 130)
  out <- data.frame(dx = numeric(length(which)), dy = numeric(length(which)),
    far = logical(length(which)))
  for (k in seq_along(which)) {
    i <- which[k]
    wk <- w[k]
    hk <- h[k]
    best <- NULL
    best.cost <- Inf
    for (ring in seq_along(gaps)) {
      r <- radius[i] + gaps[ring]
      dd <- r * 0.7
      dx <- c(r, -r - wk, dd, -dd - wk, dd, -dd - wk, -wk / 2, -wk / 2)
      dy <- c(-hk / 2, -hk / 2, dd, dd, -dd - hk, -dd - hk, r, -r - hk)
      for (s in seq_along(dx)) {
        box <- c(px[i] + dx[s], py[i] + dy[s], px[i] + dx[s] + wk,
          py[i] + dy[s] + hk)
        ow <- pmin(box[3], placed[, 3]) - pmax(box[1], placed[, 1])
        oh <- pmin(box[4], placed[, 4]) - pmax(box[2], placed[, 2])
        hit.labels <- sum(ifelse(ow > 0 & oh > 0, ow * oh, 0))
        on.dots <- pmin(box[3], dots[, 3]) > pmax(box[1], dots[, 1]) &
          pmin(box[4], dots[, 4]) > pmax(box[2], dots[, 2])
        hit.dots <- sum(on.dots, na.rm = TRUE) - isTRUE(on.dots[i])
        outside <- max(frame[1] - box[1], 0) + max(box[3] - frame[3], 0) +
          max(frame[2] - box[2], 0) + max(box[4] - frame[4], 0)
        cost <- 1e4 * hit.labels + 1e3 * outside + 5 * hit.dots + ring - 1
        if (cost < best.cost) {
          best <- c(dx[s], dy[s], ring > 1)
          best.cost <- cost
        }
      }
      if (best.cost < ring) break
    }
    placed <- rbind(placed, c(px[i] + best[1], py[i] + best[2],
      px[i] + best[1] + wk, py[i] + best[2] + hk))
    out$dx[k] <- best[1]
    out$dy[k] <- best[2]
    out$far[k] <- as.logical(best[3])
  }
  out
}

#-----------------------------------------------------------------------------
# Layers

## Triangles for motifs with infinite enrichment, pointing off the edge of
## the plot. They take the same sizes as the round points.
GeomEnrichEdges <- ggplot2::ggproto("GeomEnrichEdges", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(size = 1.5, fill = "#B4B4B4", direction = -1),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(data, panel_params, coord, stroke = 0.4) {
    coords <- coord$transform(data, panel_params)
    ## A little larger than the circles, so that they look as heavy.
    r <- 1.15 * (coords$size * ggplot2::.pt +
      stroke * ggplot2::.stroke / 2) / 2
    n <- nrow(coords)
    dx <- as.vector(rbind(-r, 0.6 * r, 0.6 * r) *
      rep(-coords$direction, each = 3L))
    dy <- as.vector(rbind(0, 0.9 * r, -0.9 * r))
    grid::polygonGrob(
      x = grid::unit(rep(coords$x, each = 3L), "npc") + grid::unit(dx, "pt"),
      y = grid::unit(rep(coords$y, each = 3L), "npc") + grid::unit(dy, "pt"),
      id = rep(seq_len(n), each = 3L),
      gp = grid::gpar(fill = coords$fill, col = "white",
        lwd = stroke * ggplot2::.stroke / 2)
    )
  }
)

layer_enrich_edges <- function(mapping, data, stroke, ...) {
  ggplot2::layer(geom = GeomEnrichEdges, stat = "identity", data = data,
    mapping = mapping, position = "identity", show.legend = FALSE,
    inherit.aes = TRUE, params = list(stroke = stroke, ...))
}

## Labels for the chosen points. The layer gets every point, so that labels
## can be kept off them, and the label placement is left to draw time when
## the size of the panel is known.
GeomEnrichLabels <- ggplot2::ggproto("GeomEnrichLabels", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(size = 1.5, label = NA, rank = NA),
  draw_key = ggplot2::draw_key_blank,
  draw_panel = function(data, panel_params, coord, logos = list(),
    text.size = 9, logo.height = 18, stroke = 0.4) {
    coords <- coord$transform(data, panel_params)
    grid::gTree(points = coords[, c("x", "y", "size", "label", "rank")],
      logos = logos, text.size = text.size, logo.height = logo.height,
      stroke = stroke, cl = "um_enrich_labels")
  }
)

layer_enrich_labels <- function(mapping, data, logos, text.size, logo.height,
  stroke, ...) {
  ggplot2::layer(geom = GeomEnrichLabels, stat = "identity", data = data,
    mapping = mapping, position = "identity", show.legend = FALSE,
    inherit.aes = TRUE, params = list(logos = logos, text.size = text.size,
      logo.height = logo.height, stroke = stroke, ...))
}

#' @exportS3Method grid::makeContent
#' @noRd
makeContent.um_enrich_labels <- function(x) {

  pts <- x$points
  W <- grid::convertWidth(grid::unit(1, "npc"), "pt", valueOnly = TRUE)
  H <- grid::convertHeight(grid::unit(1, "npc"), "pt", valueOnly = TRUE)
  px <- pts$x * W
  py <- pts$y * H
  radius <- (pts$size * ggplot2::.pt + x$stroke * ggplot2::.stroke / 2) / 2
  chosen <- which(!is.na(pts$rank))
  chosen <- chosen[order(pts$rank[chosen])]
  if (!length(chosen)) return(grid::setChildren(x, grid::gList()))

  fs <- x$text.size
  gp <- grid::gpar(fontsize = fs)
  labs <- as.character(pts$label[chosen])
  tw <- vapply(labs, function(l) grid::convertWidth(grid::grobWidth(
    grid::textGrob(l, gp = gp)), "pt", valueOnly = TRUE), numeric(1))
  th <- rep(1.2 * fs, length(labs))
  logo.h <- x$logo.height
  gap <- 0.3 * fs
  lw <- vapply(labs, function(l) {
    if (is.null(x$logos[[l]])) 0 else 0.4 * logo.h * x$logos[[l]]$ncol
  }, numeric(1))
  bw <- ifelse(lw > 0, tw + gap + lw, tw)
  bh <- ifelse(lw > 0, pmax(th, logo.h), th)

  spots <- enrich_place_labels(px, py, radius, bw, bh, chosen, c(0, 0, W, H))

  kids <- list()
  for (k in seq_along(chosen)) {
    i <- chosen[k]
    dx <- spots$dx[k]
    dy <- spots$dy[k]
    ## The text takes the side of the box nearer the point, and the logo
    ## the side further away.
    text.right <- lw[k] > 0 && dx + bw[k] / 2 < 0
    tx <- px[i] + dx + if (text.right) lw[k] + gap else 0
    ty <- py[i] + dy + bh[k] / 2
    if (spots$far[k]) {
      sx <- min(max(px[i], tx), tx + tw[k])
      sy <- min(max(py[i], ty - th[k] / 2), ty + th[k] / 2)
      len <- sqrt((px[i] - sx)^2 + (py[i] - sy)^2)
      shrink <- radius[i] + 1
      if (len > shrink + 1) {
        ex <- px[i] + (sx - px[i]) * shrink / len
        ey <- py[i] + (sy - py[i]) * shrink / len
        kids[[length(kids) + 1L]] <- grid::segmentsGrob(
          grid::unit(sx, "pt"), grid::unit(sy, "pt"),
          grid::unit(ex, "pt"), grid::unit(ey, "pt"),
          gp = grid::gpar(col = "#7F7F7F", lwd = 0.6))
      }
    }
    kids[[length(kids) + 1L]] <- grid::textGrob(labs[k],
      x = grid::unit(tx, "pt"), y = grid::unit(ty, "pt"), hjust = 0,
      vjust = 0.5, gp = gp)
    if (lw[k] > 0) {
      logo <- x$logos[[labs[k]]]
      lx <- px[i] + if (text.right) dx else dx + tw[k] + gap
      ly <- py[i] + dy + (bh[k] - logo.h) / 2
      kids[[length(kids) + 1L]] <- grid::polygonGrob(
        x = grid::unit(lx + logo$x * lw[k], "pt"),
        y = grid::unit(ly + logo$y * logo.h, "pt"),
        id = logo$id, gp = grid::gpar(fill = logo$fill, col = NA))
    }
  }

  grid::setChildren(x, do.call(grid::gList, kids))

}
