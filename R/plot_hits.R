#' Draw motif hits along sequences.
#'
#' Each sequence is drawn as a thin horizontal line (a track), and each hit
#' from [scan_sequences()] or [scan_sequences_lite()] is drawn on it from its
#' start to its end as an arrow pointing along its strand, or as a block when
#' it has none. Hits that overlap are stacked in lanes below the line so that
#' none hide another, and are coloured by motif.
#'
#' @param hits `DataFrame`, `data.frame` or `GRanges` Hits from
#'   [scan_sequences()] or [scan_sequences_lite()], with `motif`, `start` and
#'   `stop` (or `end`) columns and a `sequence` column (the sequence names of
#'   a `GRanges`). A `strand` column sets the direction of the arrows.
#'   Positions are 1-based and inclusive; hits on the reverse strand may have
#'   their start after their stop, as [scan_sequences()] reports them.
#' @param sequences `XStringSet` or named `numeric` The sequences scanned, or
#'   their lengths. Every sequence is then drawn at its full length in the
#'   order given, including those without hits. Otherwise only the sequences
#'   with hits are drawn, each up to its last hit. Unnamed sequences are named
#'   by their position in the set, as [scan_sequences()] names them. Not used
#'   together with `ranges`.
#' @param ranges `GRanges` or `data.frame` Genomic ranges to draw, for hits
#'   found on chromosomes (for example by scanning a `BSgenome`). A
#'   `data.frame` needs `seqnames` (or `chrom`), `start` and `end` columns,
#'   with positions 1-based and inclusive. Each range becomes a track showing
#'   the hits that overlap it, labelled by a `name` column (or the names of a
#'   `GRanges`) when there is one and as `chrom:start-end` otherwise. The
#'   x-axis shows genomic positions for a single range, and positions from
#'   the start of each range for several.
#' @param names `character` Draw only the tracks with these names, in this
#'   order.
#' @param max.tracks `numeric(1)` Draw at most this many tracks, or all of
#'   them with `NULL`. A note above the plot says when tracks were left out.
#' @param motifs `character`, `universalmotif` or `list` The motifs to draw,
#'   by name or as the motifs scanned (anything [convert_motifs()] accepts).
#'   This sets their order in the legend, and motifs given as objects have
#'   their logos added to it. Duplicated names are made unique as
#'   [scan_sequences_lite()] makes them.
#' @param colours `character` A vector of colours named by motif. By default
#'   motifs take the default \pkg{ggplot2} hues, spread over the motifs in
#'   the order they were scanned (or given in `motifs`).
#' @param shade.by,size.by `character(1)` A numeric column of `hits` to set
#'   the opacity or the height of the hits by, from the smallest value to the
#'   largest, such as `"score"` or `"pvalue"`. P-values, Q-values and E-values
#'   are taken as `-log10()`, so that stronger hits are always darker and
#'   taller.
#' @param labels `logical(1)` Name the motif on each hit long enough to hold
#'   it.
#' @param legend `logical(1)` Draw a legend of the motifs (and of `shade.by`
#'   and `size.by`) beside the plot.
#' @param text.size `numeric(1)` Size in points of the track names, axis text
#'   and legend. Labels on hits are a little smaller.
#'
#' @return A `ggplot` object.
#'
#' @details
#' Arrow heads, the narrowest hit drawn, the gap kept between hits in a lane
#' and which hits are wide enough to be labelled are worked out for a plotting
#' area about 6.5 inches wide, with each lane about 0.17 inches tall. The plot
#' reads best at roughly those proportions, so a height of about
#' `0.17 * (number of lanes) + 1` inches suits it when saving with
#' [ggplot2::ggsave()].
#'
#' With `ranges`, a hit that runs past the end of a range is cut at that end,
#' and a hit falling in several overlapping ranges is drawn in each of them.
#'
#' @examples
#' data(ArabidopsisPromoters)
#' data(ArabidopsisMotif)
#' motifs <- c(ArabidopsisMotif, create_motif("CACGTG", name = "G-box"))
#' hits <- scan_sequences(motifs, ArabidopsisPromoters[1:4], threshold = 1e-4,
#'   threshold.type = "pvalue", nthreads = 1)
#' plot_hits(hits, ArabidopsisPromoters[1:4], motifs = motifs)
#'
#' ## Lite hits, shaded by P-value and without labels:
#' hits2 <- scan_sequences_lite(motifs, ArabidopsisPromoters[1:4],
#'   pvalue = 1e-3, return.granges = FALSE)
#' plot_hits(hits2, ArabidopsisPromoters[1:4], shade.by = "pvalue",
#'   labels = FALSE)
#'
#' @seealso [scan_sequences()], [scan_sequences_lite()], [geom_motif()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_hits <- function(hits, sequences = NULL, ranges = NULL, names = NULL,
  max.tracks = 25, motifs = NULL, colours = NULL, shade.by = NULL,
  size.by = NULL, labels = TRUE, legend = TRUE, text.size = 9) {

  track.names <- names
  if (!is.null(max.tracks) && (!is.numeric(max.tracks) ||
      length(max.tracks) != 1L || is.na(max.tracks) || max.tracks < 1))
    stop("`max.tracks` must be at least 1, or NULL for every track",
      call. = FALSE)
  if (!is.null(track.names) && !is.character(track.names))
    stop("`names` must be a character vector", call. = FALSE)
  if (!isTRUEorFALSE(labels)) stop("`labels` must be TRUE or FALSE",
    call. = FALSE)
  if (!isTRUEorFALSE(legend)) stop("`legend` must be TRUE or FALSE",
    call. = FALSE)
  if (!is.numeric(text.size) || length(text.size) != 1L ||
      is.na(text.size) || text.size <= 0)
    stop("`text.size` must be a single positive number", call. = FALSE)

  hits <- plot_hits_table(hits)

  ## The motifs drawn, their order and their logos.
  logos <- list()
  if (!is.null(motifs)) {
    if (is.character(motifs)) {
      motif.names <- unique(motifs)
    } else {
      motifs <- convert_motifs(motifs)
      if (!is.list(motifs)) motifs <- list(motifs)
      given <- make.unique(vapply(motifs, function(x) x@name, character(1)))
      logos <- stats::setNames(motifs, given)
      motif.names <- given
    }
    palette <- plot_hits_colours(seq_along(motif.names))
    hits <- hits[hits$motif %in% motif.names, , drop = FALSE]
  } else {
    ## Colours follow the motifs' places in the scan, so that they stay the
    ## same whichever tracks are drawn.
    if ("motif.i" %in% colnames(hits) && !anyNA(hits$motif.i)) {
      first <- !duplicated(hits$motif)
      o <- order(hits$motif.i[first])
      motif.names <- hits$motif[first][o]
      palette <- plot_hits_colours(hits$motif.i[first][o],
        max(0L, hits$motif.i[first]))
    } else {
      motif.names <- unique(hits$motif)
      palette <- plot_hits_colours(seq_along(motif.names))
    }
  }
  base::names(palette) <- motif.names
  if (!is.null(colours)) {
    colours <- unlist(colours)
    missing.cols <- setdiff(motif.names, base::names(colours))
    if (length(missing.cols))
      stop("`colours` has no colours for: ",
        paste(missing.cols, collapse = ", "), call. = FALSE)
    palette <- stats::setNames(as.character(colours[motif.names]),
      motif.names)
  }

  tr <- plot_hits_tracks(hits, sequences, ranges)
  tracks <- tr$tracks
  if (!is.null(track.names)) {
    unknown <- setdiff(track.names, tracks)
    if (length(unknown))
      stop("no sequences or ranges named: ", paste(unknown, collapse = ", "),
        call. = FALSE)
    chosen <- match(track.names, tracks)
  } else {
    chosen <- seq_along(tracks)
  }
  total <- length(chosen)
  if (!is.null(max.tracks)) chosen <- utils::head(chosen, max.tracks)
  if (!length(chosen))
    stop("`hits` has no rows to plot; give `sequences` to draw tracks ",
      "without hits", call. = FALSE)
  keep <- tr$hit.track %in% chosen
  hits <- tr$hits[keep, , drop = FALSE]
  track <- match(tr$hit.track[keep], chosen)
  tracks <- tracks[chosen]
  offsets <- tr$offsets[chosen]
  lo <- tr$lo[chosen]
  hi <- tr$hi[chosen]
  ntracks <- length(tracks)
  if (is.null(motifs)) motif.names <- motif.names[motif.names %in% hits$motif]

  ## Hits running off the end of a range are cut at its ends.
  s <- pmax(hits$start - 0.5 - offsets[track], lo[track])
  e <- pmin(hits$end + 0.5 - offsets[track], hi[track])
  nhits <- nrow(hits)
  frac.shade <- frac.size <- rep(1, nhits)
  shade <- size <- NULL
  if (!is.null(shade.by)) {
    shade <- plot_hits_values(hits, shade.by, "shade.by")
    frac.shade <- shade$frac
  }
  if (!is.null(size.by)) {
    size <- plot_hits_values(hits, size.by, "size.by")
    frac.size <- size$frac
  }
  alpha <- 0.25 + 0.75 * frac.shade
  height <- HITS_HEIGHT * (0.35 + 0.65 * frac.size)

  xmin <- min(lo)
  xmax <- max(hi)
  if (xmax <= xmin) xmax <- xmin + 1
  xpad <- 0.01 * (xmax - xmin)
  per.pt <- (xmax - xmin + 2 * xpad) / HITS_PANEL_PTS

  ## Arrow heads, the narrowest hit drawn and the gap between hits in a lane
  ## are set in points, for a panel of the nominal width.
  min.width <- 3 * per.pt
  centre <- (s + e) / 2
  narrow <- e - s < min.width
  ds <- ifelse(narrow, centre - min.width / 2, s)
  de <- ifelse(narrow, centre + min.width / 2, e)
  lane <- integer(nhits)
  nlanes <- rep(1L, ntracks)
  for (k in unique(track)) {
    idx <- which(track == k)
    ln <- plot_hits_lanes(ds[idx], de[idx], 2 * per.pt)
    lane[idx] <- ln$lane
    nlanes[k] <- max(ln$n, 1L)
  }
  top <- c(0, cumsum(nlanes + HITS_TRACK_GAP))[seq_len(ntracks)]
  units <- top[ntracks] + nlanes[ntracks]

  ## Lanes run downwards, so y is drawn negated.
  y <- -(top[track] + lane)
  head <- pmin(0.5 * (de - ds), 0.45 * height * HITS_LANE_PTS * per.pt)
  fwd <- hits$strand == "+"
  bwd <- hits$strand == "-"
  y.lo <- y - height / 2
  y.hi <- y + height / 2
  body.s <- ifelse(bwd, ds + head, ds)
  body.e <- ifelse(fwd, de - head, de)
  tip.l <- ifelse(bwd, ds, body.s)
  tip.r <- ifelse(fwd, de, body.e)
  poly <- data.frame(
    id = rep(seq_len(nhits), each = 6L),
    x = as.vector(rbind(body.s, body.e, tip.r, body.e, body.s, tip.l)),
    y = as.vector(rbind(y.lo, y.lo, y, y.hi, y.hi, y)),
    motif = factor(rep(hits$motif, each = 6L), levels = motif.names)
  )
  if (!is.null(shade)) poly$shade <- rep(shade$value, each = 6L)

  lines <- data.frame(lo = lo, hi = hi, y = -top)

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = lines,
      ggplot2::aes(x = .data$lo, xend = .data$hi, y = .data$y,
        yend = .data$y),
      colour = "grey55", linewidth = 0.3)

  hit.mapping <- ggplot2::aes(x = .data$x, y = .data$y, group = .data$id,
    fill = .data$motif)
  if (!is.null(shade)) hit.mapping$alpha <- quote(.data$shade)
  key.motif <- plot_hits_key_swatch(text.size, logos)
  p <- p + ggplot2::geom_polygon(data = poly, mapping = hit.mapping,
    colour = NA, show.legend = c(fill = TRUE, alpha = FALSE),
    key_glyph = key.motif)

  if (labels && nhits) {
    fs <- min(0.85 * text.size, 0.7 * HITS_HEIGHT * HITS_LANE_PTS)
    room <- (body.e - body.s) / per.pt
    need <- 0.6 * fs * nchar(hits$motif) + 4
    tall <- height * HITS_LANE_PTS >= 1.05 * fs
    ## Text goes black or white, whichever reads better on the colour as it
    ## shows against the white page.
    rgb <- grDevices::col2rgb(palette[hits$motif]) / 255
    shown <- t(rgb) * alpha + (1 - alpha)
    light <- as.vector(shown %*% c(0.299, 0.587, 0.114)) > 0.55
    show <- which(room >= need & tall & fs >= 4)
    if (length(show)) {
      lab <- data.frame(x = (body.s[show] + body.e[show]) / 2, y = y[show],
        label = hits$motif[show])
      p <- p + ggplot2::geom_text(data = lab,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = fs / ggplot2::.pt, colour = ifelse(light[show], "black",
          "white"), inherit.aes = FALSE)
    }
  }

  p <- p + ggplot2::scale_fill_manual(values = palette, limits = motif.names,
    breaks = motif.names, name = NULL, drop = FALSE,
    guide = ggplot2::guide_legend(order = 1))

  ## The logos are drawn by the keys of a layer of invisible points whose
  ## colour scale maps each motif to its own name, which the key then looks
  ## the logo up by. The colour and fill legends share their title and labels,
  ## so they merge.
  if (length(logos)) {
    p <- p + ggplot2::geom_point(
        data = data.frame(x = xmin, y = 0,
          motif = factor(motif.names, levels = motif.names)),
        mapping = ggplot2::aes(x = .data$x, y = .data$y,
          colour = .data$motif),
        shape = NA, na.rm = TRUE, inherit.aes = FALSE,
        key_glyph = plot_hits_key_logo(text.size, logos)) +
      ggplot2::scale_colour_manual(
        values = stats::setNames(motif.names, motif.names),
        limits = motif.names, breaks = motif.names, name = NULL,
        drop = FALSE, guide = ggplot2::guide_legend(order = 1))
  }

  ## Value legends come from another layer of invisible points, so that their
  ## keys can be drawn as triangles rather than as hits.
  value.mapping <- ggplot2::aes(x = .data$x, y = .data$y)
  value.data <- data.frame(x = rep(xmin, nhits), y = rep(0, nhits))
  if (!is.null(shade) && shade$hi > shade$lo) {
    p <- p + ggplot2::scale_alpha_continuous(range = c(0.25, 1),
      limits = c(shade$lo, shade$hi), breaks = shade$breaks,
      labels = shade$labels, name = shade$title,
      guide = ggplot2::guide_legend(order = 2))
    value.mapping$alpha <- quote(.data$shade)
    value.data$shade <- shade$value
  } else if (!is.null(shade)) {
    p <- p + ggplot2::scale_alpha_continuous(range = c(1, 1), guide = "none")
  }
  if (!is.null(size) && size$hi > size$lo) {
    ## A legend for the same column as shade.by merges with it.
    same <- "alpha" %in% base::names(value.mapping) &&
      identical(shade.by, size.by)
    p <- p + ggplot2::scale_size_continuous(range = c(0.35, 1),
      limits = c(size$lo, size$hi), breaks = size$breaks,
      labels = size$labels, name = size$title,
      guide = ggplot2::guide_legend(order = if (same) 2 else 3))
    value.mapping$size <- quote(.data$size)
    value.data$size <- size$value
  }
  if (length(value.mapping) > 2L && nhits) {
    p <- p + ggplot2::geom_point(data = value.data, mapping = value.mapping,
      shape = NA, na.rm = TRUE, inherit.aes = FALSE,
      key_glyph = plot_hits_key_value(text.size,
        "size" %in% base::names(value.mapping)))
  }

  note <- NULL
  if (total > length(chosen))
    note <- paste0("Showing the first ", length(chosen), " of ", total, " ",
      if (is.null(ranges)) "sequences" else "ranges",
      "; set max.tracks to draw more")

  p +
    ggplot2::scale_x_continuous(labels = plot_hits_comma) +
    ggplot2::scale_y_continuous(breaks = -top, labels = tracks) +
    ggplot2::coord_cartesian(xlim = c(xmin - xpad, xmax + xpad),
      ylim = c(-(units - 0.4), 0.6), expand = FALSE) +
    ggplot2::labs(x = tr$xlab, y = NULL, subtitle = note) +
    ggplot2::theme_classic(base_size = text.size) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = text.size, colour = "black"),
      axis.title.x = ggplot2::element_text(size = text.size),
      plot.subtitle = ggplot2::element_text(size = 0.85 * text.size,
        colour = "grey40"),
      legend.position = if (legend) "right" else "none",
      legend.justification = "top",
      legend.box.just = "left",
      legend.text = ggplot2::element_text(size = text.size),
      legend.title = ggplot2::element_text(size = text.size),
      legend.key.size = grid::unit(0.7 * text.size, "pt"),
      legend.key.spacing.y = grid::unit(0.3 * text.size, "pt")
    )

}

#-----------------------------------------------------------------------------
# Internal helpers

## Each lane is one unit tall and hits fill most of it, tracks are kept apart
## by a gap, and sizes in points assume a panel of the nominal width with
## lanes of the nominal height.
HITS_HEIGHT <- 0.8
HITS_TRACK_GAP <- 1.2
HITS_LANE_PTS <- 0.17 * 72
HITS_PANEL_PTS <- 6.5 * 72

## Colours for motifs by their 1-based places among n, the default hues of
## ggplot2 (those of scales::hue_pal()).
plot_hits_colours <- function(places, n = length(places)) {
  if (n == 0) return(character(0))
  hues <- seq(15, 375, length.out = n + 1L)[seq_len(n)]
  grDevices::hcl(h = hues[as.integer(places)], c = 100, l = 65)
}

plot_hits_comma <- function(v) {
  out <- formatC(v, format = "f", digits = 0, big.mark = ",")
  out[is.na(v)] <- NA
  out
}

## Put hits from scan_sequences() or scan_sequences_lite() in a plain
## data.frame with 1-based start <= end and a strand of "+", "-" or "*".
plot_hits_table <- function(hits) {

  if (inherits(hits, "GRanges")) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE))
      stop("`hits` is a GRanges; the GenomicRanges package must be installed",
        call. = FALSE)
    df <- as.data.frame(hits)
    df$sequence <- as.character(df$seqnames)
    df$seqnames <- NULL
  } else if (is(hits, "DataFrame") || is.data.frame(hits)) {
    df <- as.data.frame(hits)
  } else {
    stop("`hits` must be a DataFrame, data.frame or GRanges from ",
      "scan_sequences() or scan_sequences_lite()", call. = FALSE)
  }

  if (!"end" %in% colnames(df) && "stop" %in% colnames(df))
    df$end <- df$stop
  for (col in c("motif", "start", "end")) {
    if (!col %in% colnames(df))
      stop("`hits` must have a `", col, "` column",
        if (col == "end") " (or `stop`)", call. = FALSE)
  }
  if (!"sequence" %in% colnames(df))
    stop("`hits` must have a `sequence` column", call. = FALSE)
  if (!is.numeric(df$start) || !is.numeric(df$end))
    stop("`start` and `end` must be numeric", call. = FALSE)
  if (anyNA(df$start) || anyNA(df$end))
    stop("`start` and `end` cannot contain NA values", call. = FALSE)

  start <- pmin(df$start, df$end)
  df$end <- pmax(df$start, df$end)
  df$start <- start
  df$motif <- as.character(df$motif)
  df$sequence <- as.character(df$sequence)
  strand <- if ("strand" %in% colnames(df)) as.character(df$strand)
    else rep("*", nrow(df))
  strand[is.na(strand) | !strand %in% c("+", "-")] <- "*"
  df$strand <- strand
  df

}

## A numeric column scaled from 0 at its smallest value to 1 at its largest,
## with P-values and the like taken as -log10, plus what its legend needs.
plot_hits_values <- function(hits, col, arg) {

  if (!is.character(col) || length(col) != 1L || is.na(col))
    stop("`", arg, "` must be a single column name", call. = FALSE)
  if (!col %in% colnames(hits))
    stop("`hits` has no `", col, "` column to use for `", arg, "`",
      call. = FALSE)
  v <- hits[[col]]
  if (!is.numeric(v))
    stop("`", arg, "` must name a numeric column", call. = FALSE)
  v <- as.numeric(v)
  title <- col
  if (grepl("^(pval|qval|eval)", tolower(col))) {
    v <- -log10(pmax(v, 1e-300))
    title <- paste0("-log10(", col, ")")
  }
  finite <- v[is.finite(v)]
  if (length(finite)) {
    lo <- min(finite)
    hi <- max(finite)
  } else {
    lo <- 0
    hi <- 1
  }
  ## Missing values count as the weakest and infinite ones are clipped.
  v[is.na(v)] <- lo
  v <- pmin(pmax(v, lo), hi)
  frac <- if (hi <= lo) rep(1, length(v)) else (v - lo) / (hi - lo)
  breaks <- pretty(c(lo, hi), n = 3)
  breaks <- breaks[breaks >= lo & breaks <= hi]
  if (!length(breaks)) breaks <- c(lo, hi)
  list(value = v, frac = frac, lo = lo, hi = hi, title = title,
    breaks = breaks, labels = vapply(breaks, format, character(1)))

}

## The tracks to draw: their labels, the position each is drawn from, their
## drawn bounds, the track of each hit, the hits and the x-axis title.
plot_hits_tracks <- function(hits, sequences, ranges) {

  if (!is.null(ranges)) return(plot_hits_range_tracks(hits, ranges))

  lengths <- NULL
  if (!is.null(sequences)) {
    if (is(sequences, "XString"))
      sequences <- methods::as(sequences, paste0(class(sequences), "Set"))
    if (is(sequences, "XStringSet")) {
      lengths <- Biostrings::width(sequences)
      nms <- names(sequences)
      names(lengths) <- if (is.null(nms))
        as.character(seq_along(sequences)) else nms
    } else if (is.numeric(sequences) && !is.null(names(sequences))) {
      lengths <- sequences
    } else {
      stop("`sequences` must be an XStringSet or a named numeric vector of ",
        "sequence lengths", call. = FALSE)
    }
  }

  if (!is.null(lengths)) {
    tracks <- names(lengths)
    if (nrow(hits) && !any(hits$sequence %in% tracks))
      stop("none of the hits are on the sequences given", call. = FALSE)
  } else if ("sequence.i" %in% colnames(hits) && !anyNA(hits$sequence.i)) {
    tracks <- unique(hits$sequence[order(hits$sequence.i)])
  } else {
    tracks <- unique(hits$sequence)
  }
  hit.track <- match(hits$sequence, tracks)
  on <- !is.na(hit.track)
  hits <- hits[on, , drop = FALSE]
  hit.track <- hit.track[on]
  last <- numeric(length(tracks))
  if (length(hit.track)) {
    m <- tapply(hits$end, hit.track, max)
    last[as.integer(names(m))] <- m
  }
  if (!is.null(lengths)) {
    size <- as.numeric(lengths)
    over <- which(last > size)
    if (length(over))
      stop("hits run past the end of ", tracks[over[1]], call. = FALSE)
  } else {
    size <- last
  }

  list(tracks = tracks, offsets = numeric(length(tracks)),
    lo = rep(0.5, length(tracks)), hi = size + 0.5, hit.track = hit.track,
    hits = hits, xlab = "Position")

}

plot_hits_range_tracks <- function(hits, ranges) {

  labs <- NULL
  if (inherits(ranges, "GRanges")) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE))
      stop("`ranges` is a GRanges; the GenomicRanges package must be ",
        "installed", call. = FALSE)
    labs <- names(ranges)
    ranges <- as.data.frame(ranges)
    if (is.null(labs) && "name" %in% colnames(ranges)) labs <- ranges$name
  } else if (is(ranges, "DataFrame") || is.data.frame(ranges)) {
    ranges <- as.data.frame(ranges)
    if ("name" %in% colnames(ranges)) labs <- ranges$name
  } else {
    stop("`ranges` must be a GRanges or a data.frame", call. = FALSE)
  }
  chrom.col <- intersect(c("seqnames", "chrom"), colnames(ranges))
  if (!length(chrom.col) || !all(c("start", "end") %in% colnames(ranges)))
    stop("`ranges` must have `seqnames` (or `chrom`), `start` and `end` ",
      "columns", call. = FALSE)
  if (!nrow(ranges)) stop("`ranges` has no rows to plot", call. = FALSE)
  chrom <- as.character(ranges[[chrom.col[1]]])
  rs <- as.numeric(ranges$start)
  re <- as.numeric(ranges$end)
  if (anyNA(rs) || anyNA(re) || any(re < rs))
    stop("`ranges` must have start <= end", call. = FALSE)
  labs <- if (is.null(labs)) paste0(chrom, ":", rs, "-", re)
    else as.character(labs)

  ## Hits in overlapping ranges are drawn once in each, and hits found more
  ## than once (such as from scanning overlapping ranges) only once.
  pick <- lapply(seq_along(rs), function(i) {
    idx <- which(hits$sequence == chrom[i] & hits$start <= re[i] &
      hits$end >= rs[i])
    idx[!duplicated(hits[idx, c("motif", "start", "end", "strand")])]
  })
  hit.track <- rep(seq_along(pick), lengths(pick))

  if (length(rs) == 1L) {
    offsets <- 0
    xlab <- paste("Position on", chrom)
  } else {
    offsets <- rs - 1
    xlab <- "Position from the start of each range"
  }

  list(tracks = labs, offsets = offsets, lo = rs - 0.5 - offsets,
    hi = re + 0.5 - offsets, hit.track = hit.track,
    hits = hits[unlist(pick), , drop = FALSE], xlab = xlab)

}

## Put intervals in lanes so that none in a lane overlap or come within `gap`
## of each other, each going in the lowest lane free at its start.
plot_hits_lanes <- function(starts, ends, gap) {

  lane <- integer(length(starts))
  last <- numeric(0)
  for (i in order(starts, ends)) {
    free <- which(last + gap <= starts[i])
    if (length(free)) {
      k <- free[1]
    } else {
      k <- length(last) + 1L
    }
    last[k] <- ends[i]
    lane[i] <- k - 1L
  }
  list(lane = lane, n = length(last))

}

## Sizes of the motif legend keys in points: a swatch of the motif's colour,
## followed by its logo when there are logos.
plot_hits_key_dims <- function(text.size, logos) {
  has.logos <- length(logos) > 0L
  widest <- if (has.logos) max(vapply(logos, function(x) ncol(x@motif),
    numeric(1))) else 1
  ## Logos are about a third as wide per position as they are tall.
  logo.width <- 0.75 * widest * text.size
  sw <- 1.2 * text.size
  list(swatch.w = sw, swatch.h = 0.7 * text.size,
    logo.x = sw + 0.4 * text.size, logo.width = logo.width, widest = widest,
    width = if (has.logos) sw + 0.8 * text.size + logo.width else sw,
    height = if (has.logos) 2.2 * text.size else 0.7 * text.size)
}

plot_hits_pt_cm <- function(x) x / 72.27 * 2.54

plot_hits_key_swatch <- function(text.size, logos) {
  d <- plot_hits_key_dims(text.size, logos)
  function(data, params, size) {
    fill <- if (is.null(data$fill) || is.na(data$fill)) "grey30" else
      data$fill
    g <- grid::rectGrob(x = grid::unit(0, "npc"),
      width = grid::unit(d$swatch.w, "pt"),
      height = grid::unit(d$swatch.h, "pt"), just = "left",
      gp = grid::gpar(fill = fill, col = NA))
    attr(g, "width") <- plot_hits_pt_cm(d$width)
    attr(g, "height") <- plot_hits_pt_cm(d$height)
    g
  }
}

plot_hits_key_logo <- function(text.size, logos) {
  d <- plot_hits_key_dims(text.size, logos)
  cache <- new.env(parent = emptyenv())
  function(data, params, size) {
    key <- as.character(data$colour)
    m <- logos[[key]]
    if (is.null(m)) return(grid::nullGrob())
    if (is.null(cache[[key]]))
      cache[[key]] <- motif_to_polygons(m, "ICM", fontDFroboto, NULL, 0.01,
        0.04, 0.01, TRUE, TRUE, NULL, FALSE)
    info <- cache[[key]]
    if (!nrow(info$poly)) return(grid::nullGrob())
    poly <- info$poly[order(info$poly$letter.id), , drop = FALSE]
    ids <- match(poly$letter.id, unique(poly$letter.id))
    ## Letters are scaled to the most information a position can hold, so
    ## that weak positions look weak.
    yr <- max(info$yhi, log2(nrow(m@motif))) - info$ylo
    if (!is.finite(yr) || yr <= 0) yr <- 1
    vp <- grid::viewport(x = grid::unit(d$logo.x, "pt"),
      width = grid::unit(d$logo.width * ncol(m@motif) / d$widest, "pt"),
      height = grid::unit(d$height, "pt"), just = "left")
    g <- grid::polygonGrob(x = (poly$x - info$xlo) / (info$xhi - info$xlo),
      y = (poly$y - info$ylo) / yr, id = ids,
      gp = grid::gpar(fill = poly$fill_colour[!duplicated(ids)], col = NA),
      vp = vp)
    g <- grid::gTree(children = grid::gList(g))
    attr(g, "width") <- plot_hits_pt_cm(d$width)
    attr(g, "height") <- plot_hits_pt_cm(d$height)
    g
  }
}

## Keys of the shade.by and size.by legends, drawn as right-pointing
## triangles that keep the hits' proportions but are no taller than a line
## and a half of text.
plot_hits_key_value <- function(text.size, sized) {
  hit.pts <- HITS_HEIGHT * HITS_LANE_PTS
  largest <- if (sized) min(hit.pts, 1.5 * text.size) else
    min(0.8 * hit.pts, text.size)
  function(data, params, size) {
    h <- if (sized) data$size * largest else largest
    a <- if (is.null(data$alpha) || is.na(data$alpha)) 1 else data$alpha
    g <- grid::polygonGrob(
      x = grid::unit(0.5, "npc") + grid::unit(c(-0.5, 0.5, -0.5) * h, "pt"),
      y = grid::unit(0.5, "npc") + grid::unit(c(-0.5, 0, 0.5) * h, "pt"),
      gp = grid::gpar(fill = grDevices::adjustcolor("grey30", alpha.f = a),
        col = NA))
    attr(g, "width") <- plot_hits_pt_cm(max(1.2 * text.size, largest))
    attr(g, "height") <- plot_hits_pt_cm(max(0.7 * text.size, largest))
    g
  }
}
