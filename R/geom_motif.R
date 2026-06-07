#' Embed motif logos in ggplot2 constructions.
#'
#' These functions expose the logo-drawing machinery of [view_motifs()] and
#' [view_logo()] as composable \pkg{ggplot2} layers, so a sequence logo can be
#' dropped into a plot you are already building (annotating a coverage track,
#' adding an inset to a scatter plot, placing a logo at every tip of a tree,
#' and so on) instead of being returned as a standalone plot.
#'
#' There are two pairs. `geom_motif()` / `annotate_motif()` take a
#' `universalmotif` object (or anything [convert_motifs()] accepts) and honour
#' `use.type` exactly like [view_motifs()]. `geom_logo()` / `annotate_logo()`
#' take a bare row-named numeric matrix with arbitrary values, like
#' [view_logo()]. In each pair the `annotate_*()` function places one logo at
#' fixed coordinates (no data or aesthetic inheritance), while the `geom_*()`
#' function is a true layer whose position comes from aesthetics, so it picks
#' up faceting and can be driven from a data frame.
#'
#' @param mapping Aesthetic mapping created with [ggplot2::aes()]. The layer
#'   understands two positioning modes (see Details): a bounding box
#'   (`xmin`, `xmax`, `ymin`, `ymax`) or an anchor (`x`, `y`). The optional
#'   `motif` aesthetic selects, per group, which entry of a named-list `motif`
#'   (or `logo`) argument to draw, e.g. `aes(x = 0, y = y, motif = label)`.
#' @param data A data frame, if the layer should override the plot data.
#' @param stat,position,na.rm,show.legend,inherit.aes Standard \pkg{ggplot2}
#'   layer arguments. Colours are resolved internally rather than through the
#'   plot's `fill` scale, so the layer never collides with a `fill` mapping in
#'   the host plot and no legend is added; `show.legend` therefore defaults to
#'   `FALSE`.
#' @param motif A single motif (anything [convert_motifs()] accepts) or a
#'   **named list** of motifs. With a named list, map the `motif` aesthetic to
#'   a column whose values are names of that list to draw a different logo per
#'   group.
#' @param logo A single row-named numeric matrix or a **named list** of such
#'   matrices, used the same way as `motif` but drawn with the values taken
#'   as-is (the [view_logo()] path).
#' @param xmin,xmax,ymin,ymax `numeric(1)` Bounding box in data coordinates for
#'   `annotate_motif()` / `annotate_logo()`. The logo is linearly rescaled to
#'   fill the box.
#' @param use.type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#'   The matrix type used for display (motif variants only). See
#'   [view_motifs()].
#' @param fill `character(1)` A single colour for every letter (logo variants),
#'   used when `colour.scheme` is `NULL`.
#' @param colour.scheme `character` A named vector giving a colour for every
#'   individual character. For the motif variants the default is the standard
#'   DNA/RNA/AA palette (or a generated palette for custom alphabets). For the
#'   bare-matrix logo variants the alphabet is unknown, so when `colour.scheme`
#'   is `NULL` every letter is drawn in the single `fill` colour instead.
#' @param fontDF `data.frame` or `DataFrame` Polygon paths for the letters. By
#'   default the bundled [fontDFroboto] glyphs are used.
#' @param width `numeric(1)` In anchored mode, the width of one motif position
#'   in data coordinates (so the logo spans `x` to `x + ncol * width`).
#' @param height `numeric(1)` In anchored mode, the total logo height in data
#'   coordinates, centred on `y`.
#' @param ... Passed to [ggplot2::layer()], and on to the underlying logo
#'   engine: `min.height`, `x.spacer`, `y.spacer`, `sort.positions`,
#'   `sort.positions.decreasing`, `flip.neg` (same meaning as in
#'   [view_motifs()] / [view_logo()]).
#'
#' @return A \pkg{ggplot2} layer to add to a plot with `+`.
#'
#' @details
#' Letters are drawn as filled, borderless polygons whose colour is resolved
#' internally (so the layer never competes for the plot's `fill` scale and adds
#' no legend). Consequently the `colour`, `alpha`, and `linetype` aesthetics
#' have no effect; appearance is controlled only through `colour.scheme` (or
#' `fill` for the logo variants).
#'
#' Two positioning modes are supported and chosen per row from the supplied
#' aesthetics:
#'
#' * **Bounding box**: supply `xmin`, `xmax`, `ymin`, `ymax` (the model used by
#'   [ggplot2::annotation_custom()]). The whole logo is rescaled into the box.
#' * **Anchored**: supply `x` (the left edge) and `y` (the vertical centre);
#'   the box is then `xmin = x`, `xmax = x + ncol * width`,
#'   `ymin = y - height/2`, `ymax = y + height/2`. This is the natural fit for
#'   placing logos at tree tips.
#'
#' To draw several motifs aligned to a common column frame (for example a
#' stack of aligned logos, or aligned tip logos on a tree), align them first
#' with [view_motifs2()] (or [view_motifs()]) and pass the result to
#' `geom_logo()`:
#'
#' \preformatted{
#' aligned <- view_motifs2(motifs, return.raw = TRUE, use.type = "ICM",
#'                         sort.by = "none")
#' ggplot(tipdf, aes(x = 0, y = y, motif = label)) +
#'   geom_logo(logo = aligned, colour.scheme = DNA_COLOURS)
#' }
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' data(examplemotif)
#'
#' ## Annotate a logo onto an existing plot:
#' df <- data.frame(x = 1:20, y = cumsum(rnorm(20)))
#' ggplot(df, aes(x, y)) +
#'   geom_line() +
#'   annotate_motif(examplemotif, xmin = 5, xmax = 11, ymin = 1, ymax = 4)
#'
#' ## A different motif per row, anchored at each y:
#' m1 <- create_motif("TATAAA", name = "a")
#' m2 <- create_motif("GGGCGG", name = "b")
#' tipdf <- data.frame(x = 0, y = c(1, 2), label = c("a", "b"))
#' ggplot(tipdf, aes(x = x, y = y, motif = label)) +
#'   geom_motif(motif = list(a = m1, b = m2), height = 0.8)
#' }
#'
#' @seealso [view_motifs()], [view_logo()], [view_motifs2()], [motif_tree2()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @name geom_motif
#' @export
geom_motif <- function(mapping = NULL, data = NULL, stat = StatMotif,
  position = "identity", ..., motif, use.type = "ICM", colour.scheme = NULL,
  fontDF = NULL, width = 1, height = 0.9, na.rm = FALSE, show.legend = FALSE,
  inherit.aes = TRUE) {

  layer(geom = GeomMotif, stat = stat, data = data, mapping = mapping,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(motif = motif, use.type = use.type,
      colour.scheme = colour.scheme, fontDF = fontDF, width = width,
      height = height, na.rm = na.rm, ...))

}

#' @rdname geom_motif
#' @export
geom_logo <- function(mapping = NULL, data = NULL, stat = StatLogo,
  position = "identity", ..., logo, fill = "black", colour.scheme = NULL,
  fontDF = NULL, width = 1, height = 0.9, na.rm = FALSE, show.legend = FALSE,
  inherit.aes = TRUE) {

  layer(geom = GeomMotif, stat = stat, data = data, mapping = mapping,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(logo = logo, fill = fill, colour.scheme = colour.scheme,
      fontDF = fontDF, width = width, height = height, na.rm = na.rm, ...))

}

#' @rdname geom_motif
#' @export
annotate_motif <- function(motif, xmin, xmax, ymin, ymax, use.type = "ICM",
  colour.scheme = NULL, fontDF = NULL, ...) {

  d <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  layer(geom = GeomMotif, stat = StatMotif, data = d,
    mapping = aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
      ymax = .data$ymax),
    position = "identity", show.legend = FALSE, inherit.aes = FALSE,
    params = list(motif = motif, use.type = use.type,
      colour.scheme = colour.scheme, fontDF = fontDF, na.rm = FALSE, ...))

}

#' @rdname geom_motif
#' @export
annotate_logo <- function(logo, xmin, xmax, ymin, ymax, fill = "black",
  colour.scheme = NULL, fontDF = NULL, ...) {

  d <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  layer(geom = GeomMotif, stat = StatLogo, data = d,
    mapping = aes(xmin = .data$xmin, xmax = .data$xmax, ymin = .data$ymin,
      ymax = .data$ymax),
    position = "identity", show.legend = FALSE, inherit.aes = FALSE,
    params = list(logo = logo, fill = fill, colour.scheme = colour.scheme,
      fontDF = fontDF, na.rm = FALSE, ...))

}

#-----------------------------------------------------------------------------
# ggproto objects

#' @rdname geom_motif
#' @usage NULL
#' @format NULL
#' @export
StatMotif <- ggproto("StatMotif", Stat,
  required_aes = character(0),
  default_aes = aes(xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
    x = NULL, y = NULL, motif = NULL),
  setup_data = function(data, params) {
    logo_validate_positions(data)
    keys <- logo_used_keys(data, params$motif)
    if (!is.null(params$colour.scheme)) {
      for (k in keys) {
        m <- if (identical(k, "._single_")) params$motif else params$motif[[k]]
        mm <- convert_motifs(m)
        if (is.list(mm)) mm <- mm[[1]]
        resolve_logo_colours(rownames(mm@motif), mm@alphabet,
          params$colour.scheme, NULL)
      }
    }
    data
  },
  compute_panel = function(data, scales, motif = NULL, use.type = "ICM",
    colour.scheme = NULL, fontDF = NULL, width = 1, height = 0.9,
    min.height = 0.01, x.spacer = 0.04, y.spacer = 0.01,
    sort.positions = !use.type %in% c("PCM", "PPM"),
    sort.positions.decreasing = TRUE, flip.neg = FALSE, na.rm = FALSE) {

    if (is.null(fontDF)) fontDF <- fontDFroboto
    fontDF <- as.data.frame(fontDF)
    fit.to.height <- if (use.type == "PPM") 1 else NULL

    gen <- function(key) {
      m <- if (identical(key, "._single_")) motif else motif[[key]]
      motif_to_polygons(m, use.type, fontDF, colour.scheme, min.height,
        x.spacer, y.spacer, sort.positions, sort.positions.decreasing,
        fit.to.height, flip.neg)
    }

    compute_logo_panel(data, gen, width, height)

  }
)

#' @rdname geom_motif
#' @usage NULL
#' @format NULL
#' @export
StatLogo <- ggproto("StatLogo", Stat,
  required_aes = character(0),
  default_aes = aes(xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
    x = NULL, y = NULL, motif = NULL),
  setup_data = function(data, params) {
    logo_validate_positions(data)
    keys <- logo_used_keys(data, params$logo)
    for (k in keys) {
      m <- if (identical(k, "._single_")) params$logo else params$logo[[k]]
      mm <- as.matrix(m)
      if (is.null(rownames(mm)))
        stop("`logo` must be a matrix with row names", call. = FALSE)
      if (!is.null(params$colour.scheme))
        resolve_logo_colours(
          unique(unlist(strsplit(rownames(mm), "", fixed = TRUE))),
          NULL, params$colour.scheme, NULL)
    }
    data
  },
  compute_panel = function(data, scales, logo = NULL, fill = "black",
    colour.scheme = NULL, fontDF = NULL, width = 1, height = 0.9,
    min.height = 0.01, x.spacer = 0.04, y.spacer = 0.01,
    sort.positions = FALSE, sort.positions.decreasing = TRUE,
    flip.neg = FALSE, na.rm = FALSE) {

    if (is.null(fontDF)) fontDF <- fontDFroboto
    fontDF <- as.data.frame(fontDF)

    gen <- function(key) {
      m <- if (identical(key, "._single_")) logo else logo[[key]]
      logo_to_polygons(m, fontDF, colour.scheme, fill, min.height, x.spacer,
        y.spacer, sort.positions, sort.positions.decreasing, flip.neg)
    }

    compute_logo_panel(data, gen, width, height)

  }
)

#' @rdname geom_motif
#' @usage NULL
#' @format NULL
#' @export
GeomMotif <- ggproto("GeomMotif", GeomPolygon,
  draw_key = draw_key_blank,
  draw_panel = function(data, panel_params, coord, ...) {
    coords <- coord$transform(data, panel_params)
    ## Keep each polygon's vertices contiguous (order() is stable for the
    ## integer `group`, so within-group vertex order is preserved).
    ord <- order(coords$group)
    coords <- coords[ord, , drop = FALSE]
    ids <- match(coords$group, unique(coords$group))
    first <- !duplicated(coords$group)
    grid::polygonGrob(
      x = coords$x, y = coords$y, id = ids,
      gp = grid::gpar(fill = coords$fill_colour[first], col = NA)
    )
  }
)

#-----------------------------------------------------------------------------
# Internal helpers

## Align a list of motifs to a common, zero-padded column frame and return the
## display matrices, reusing the view_motifs2() PCC aligner. `sort.by = "none"`
## keeps input order; the [RC] name suffix is dropped and names are restored
## from the caller's motif names by position. DNA/RNA only.
##
## The returned list is named by motif name so it can be used as the selector
## list for geom_logo(aes(motif = ...)). Duplicate names would make that lookup
## ambiguous, so by default duplicates are an error; set `dedup = TRUE` to
## disambiguate them with make.unique() instead.
align_motif_mats <- function(motifs, use.type = "ICM", tryRC = TRUE,
  min.overlap = 6, nthreads = 1, dedup = FALSE) {

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  nm <- vapply(motifs, function(x) x@name, character(1))

  if (anyDuplicated(nm)) {
    if (dedup) nm <- make.unique(nm)
    else stop(wmsg("Duplicate motif names (",
        paste(unique(nm[duplicated(nm)]), collapse = ", "),
        "); set `dedup = TRUE` to disambiguate them."), call. = FALSE)
  }

  mats <- view_motifs2(motifs, return.raw = TRUE, use.type = use.type,
    sort.by = "none", tryRC = tryRC, min.overlap = min.overlap,
    nthreads = nthreads)

  names(mats) <- nm
  mats

}

## Validate that every row supplies exactly one positioning mode. Run from
## setup_data() so that the error reaches the user (ggplot2 turns errors
## thrown inside compute_panel() into warnings).
logo_validate_positions <- function(data) {
  box <- c("xmin", "xmax", "ymin", "ymax")
  for (i in seq_len(nrow(data))) {
    hb <- all(box %in% names(data)) &&
      all(!is.na(unlist(data[i, box, drop = TRUE])))
    ha <- all(c("x", "y") %in% names(data)) &&
      !is.na(data$x[i]) && !is.na(data$y[i])
    if (hb && ha)
      stop(wmsg("Supply either a bounding box (xmin/xmax/ymin/ymax) or an ",
          "anchor (x/y), not both."), call. = FALSE)
    if (!hb && !ha)
      stop(wmsg("geom_motif()/geom_logo() needs either xmin/xmax/ymin/ymax ",
          "or x/y aesthetics."), call. = FALSE)
  }
  invisible(NULL)
}

## Resolve, and validate, which entries of a single object / named list are
## needed for the rows in `data`.
logo_used_keys <- function(data, obj) {
  if (!"motif" %in% names(data)) return("._single_")
  keys <- unique(as.character(data$motif))
  if (!is.list(obj) || is.null(names(obj)))
    stop(wmsg("When mapping the `motif` aesthetic, supply a named list of ",
        "motifs (or logos) to select from."), call. = FALSE)
  miss <- setdiff(keys, names(obj))
  if (length(miss))
    stop(wmsg("`motif` value(s) not found in the supplied list: ",
        paste(miss, collapse = ", ")), call. = FALSE)
  keys
}

## Resolve a colour for every individual character.
resolve_logo_colours <- function(chars, alphabet = NULL, colour.scheme = NULL,
  fill = NULL) {

  chars <- unique(chars)

  if (!is.null(colour.scheme)) {
    miss <- setdiff(chars, names(colour.scheme))
    if (length(miss))
      stop(wmsg("`colour.scheme` must be a named vector with every letter ",
          "present; missing: ", paste(miss, collapse = ", ")), call. = FALSE)
    return(colour.scheme)
  }

  def <- if (!is.null(alphabet) && length(alphabet) == 1L)
    switch(as.character(alphabet),
      "DNA" = DNA_COLOURS, "RNA" = RNA_COLOURS, "AA" = AA_COLOURS, NULL)
  else NULL
  if (!is.null(def) && all(chars %in% names(def))) return(def)

  if (!is.null(fill)) return(stats::setNames(rep(fill, length(chars)), chars))

  pal <- grDevices::hcl.colors(max(length(chars), 2L), "Dark 3")
  stats::setNames(pal[seq_along(chars)], chars)

}

## Attach fill colours and record the native coordinate extents used for
## rescaling into a box.
finish_logo_polygons <- function(poly, mat, alphabet, colour.scheme, fill) {

  poly <- poly[!is.na(poly$x), , drop = FALSE]
  ncol_mat <- ncol(mat)

  ## A motif with no letters above `min.height` (e.g. a zero-information
  ## motif drawn as ICM) yields no polygons; return empty extents rather
  ## than letting min()/max() warn on an empty vector.
  if (!nrow(poly)) {
    poly$fill_colour <- character(0)
    return(list(poly = poly, ncol = ncol_mat, xlo = 0.5, xhi = ncol_mat + 0.5,
      ylo = 0, yhi = 1))
  }

  chars <- unique(unlist(strsplit(rownames(mat), "", fixed = TRUE)))
  cols <- resolve_logo_colours(chars, alphabet, colour.scheme, fill)
  poly$fill_colour <- unname(cols[poly$group])

  list(poly = poly, ncol = ncol_mat, xlo = 0.5, xhi = ncol_mat + 0.5,
    ylo = min(poly$y), yhi = max(poly$y))

}

## Native polygons for a universalmotif object (the view_motifs() path).
motif_to_polygons <- function(motif, use.type, fontDF, colour.scheme,
  min.height, x.spacer, y.spacer, sort.positions, sort.positions.decreasing,
  fit.to.height, flip.neg) {

  motif <- convert_motifs(motif)
  if (is.list(motif)) motif <- motif[[1]]
  alphabet <- motif@alphabet
  mat <- convert_type(motif, use.type)@motif

  poly <- prep_single_motif_plot_data(mat, use.type, fontDF, min.height,
    x.spacer, y.spacer, sort.positions, sort.positions.decreasing,
    fit.to.height, 1, flip.neg = flip.neg)

  finish_logo_polygons(poly, mat, alphabet, colour.scheme, NULL)

}

## Native polygons for a bare numeric matrix (the view_logo() path).
logo_to_polygons <- function(x, fontDF, colour.scheme, fill, min.height,
  x.spacer, y.spacer, sort.positions, sort.positions.decreasing, flip.neg) {

  mat <- as.matrix(x)
  if (is.null(rownames(mat)))
    stop("`logo` must be a matrix with row names", call. = FALSE)

  poly <- make_matrix_polygon_data(mat, fontDF, min.height, x.spacer, y.spacer,
    sort.positions, sort.positions.decreasing, NULL, 1, flip.neg = flip.neg)

  finish_logo_polygons(poly, mat, NULL, colour.scheme, fill)

}

## Pick a bounding box for one input row from whichever positioning mode it
## supplies.
logo_row_box <- function(row, info, width, height) {

  full <- function(cols) all(cols %in% names(row)) &&
    all(vapply(cols, function(z) {
      v <- row[[z]]; length(v) == 1L && !is.na(v)
    }, logical(1)))

  has_box <- full(c("xmin", "xmax", "ymin", "ymax"))
  has_anchor <- full(c("x", "y"))

  if (has_box && has_anchor)
    stop(wmsg("Supply either a bounding box (xmin/xmax/ymin/ymax) or an ",
        "anchor (x/y), not both."), call. = FALSE)

  if (has_box) {
    list(xmin = row$xmin, xmax = row$xmax, ymin = row$ymin, ymax = row$ymax)
  } else if (has_anchor) {
    list(xmin = row$x, xmax = row$x + info$ncol * width,
      ymin = row$y - height / 2, ymax = row$y + height / 2)
  } else {
    stop(wmsg("geom_motif()/geom_logo() needs either xmin/xmax/ymin/ymax or ",
        "x/y aesthetics."), call. = FALSE)
  }

}

## Shared per-panel computation: one logo per input row, rescaled into its
## box, with native polygons memoised by motif key.
compute_logo_panel <- function(data, gen, width, height) {

  has_sel <- "motif" %in% names(data)
  cache <- list()
  out <- vector("list", nrow(data))

  ## NOTE (one logo per row): this draws a logo for every row of `data`. With
  ## inherit.aes = TRUE on a multi-row plot, `geom_motif(motif = m)` will place
  ## the same logo at every row's coordinates. That is the expected ggplot2
  ## "one mark per row" semantic, but can surprise users who expect a single
  ## logo; they should pass an explicit single-row `data` (as annotate_*() do).
  for (i in seq_len(nrow(data))) {

    key <- if (has_sel) as.character(data$motif[i]) else "._single_"

    if (is.null(cache[[key]])) {
      info <- tryCatch(gen(key), error = function(e) e)
      if (inherits(info, "error")) {
        if (has_sel)
          stop(wmsg("Could not build a logo for `motif` key \"", key, "\": ",
              conditionMessage(info)), call. = FALSE)
        else stop(info)
      }
      cache[[key]] <- info
    }
    info <- cache[[key]]

    ## Nothing to draw for this row (empty/zero-information logo); leave the
    ## slot NULL so rbind() drops it.
    if (!nrow(info$poly)) next

    box <- logo_row_box(data[i, , drop = FALSE], info, width, height)

    ## NOTE (y baseline not pinned): the native y-range [ylo, yhi] is mapped
    ## linearly onto [ymin, ymax], i.e. the logo is stretched to fill the box.
    ## For ICM/PPM ylo == 0 so the baseline sits at ymin, but for signed
    ## matrices (PWM/CWM) the y = 0 line lands partway up the box rather than at
    ## a fixed height. That is intended box-fill behaviour; revisit if a fixed
    ## zero baseline is ever needed.
    p <- info$poly
    p$x <- box$xmin + (p$x - info$xlo) / (info$xhi - info$xlo) *
      (box$xmax - box$xmin)
    yr <- info$yhi - info$ylo
    if (!is.finite(yr) || yr == 0) yr <- 1
    p$y <- box$ymin + (p$y - info$ylo) / yr * (box$ymax - box$ymin)

    p$group <- paste0(i, ".", p$letter.id)
    p$PANEL <- data$PANEL[i]
    out[[i]] <- p[, c("x", "y", "group", "fill_colour", "PANEL")]

  }

  res <- do.call(rbind, out)
  if (is.null(res) || !nrow(res))
    return(data.frame(x = numeric(0), y = numeric(0), group = integer(0),
      fill_colour = character(0), PANEL = factor()))
  res$group <- match(res$group, unique(res$group))
  res

}
