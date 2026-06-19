#' Generate ggtree-based motif trees using `compare_motifs_lite()`.
#'
#' `motif_tree_lite()` is the leaner counterpart of [motif_tree()]: it builds the
#' distance matrix via [compare_motifs_lite()] (mean Pearson correlation,
#' with built-in p-value / q-value machinery) instead of
#' [compare_motifs()]. The distance for `hclust()` is derived from the
#' mean Pearson correlation matrix as `(1 - score) / 2`, mapping
#' similarity in `[-1, 1]` to a symmetric distance in `[0, 1]`. The
#' tree-construction and rendering steps (`hclust()` followed by
#' [ape::as.phylo()] followed by [ggtree::ggtree()]) match
#' [motif_tree()] exactly, so the visual conventions, the `linecol` /
#' `labels` / `tipsize` slot-mapping arguments, and the `dist`-input
#' short-circuit all behave the same.
#'
#' @param motifs `list` or `dist`. See [convert_motifs()] for accepted
#'   motif formats. Alternatively, a pre-built `dist` object (e.g. the
#'   result of `as.dist((1 - compare_motifs_lite(...)) / 2)`) skips the
#'   comparison step entirely and goes straight to tree construction.
#' @param layout `character(1)`. One of `c('rectangular', 'slanted',
#'   'fan', 'circular', 'radial', 'equal_angle', 'daylight')`. Passed
#'   through to [ggtree::ggtree()]. Defaults to `"rectangular"` so that the
#'   tip logos (see `tiplogo`) are visible; tip logos are only drawn on the
#'   `"rectangular"` and `"slanted"` layouts.
#' @param linecol `character(1)`. Motif slot used to colour the
#'   branches (e.g. `"family"`). Not used when `motifs` is a `dist`.
#' @param labels `character(1)`. Motif slot used to label the tips
#'   (e.g. `"name"`). For `dist` input, only `"name"` is meaningful. Note:
#'   drawing tip labels uses `ggtree::geom_tiplab(align = TRUE)`, which fails
#'   to render under \pkg{ggplot2} >= 4.0 with \pkg{ggtree} <= 3.12 (ggtree's
#'   internal `empty()` calls the removed `ggplot2:::is.waive()`); update
#'   \pkg{ggtree} if you need tip labels. The tip logos (`tiplogo`) are
#'   unaffected.
#' @param tipsize `character(1)`. Motif slot used to size the tips
#'   (e.g. `"icscore"`). Not used when `motifs` is a `dist`.
#' @param legend `logical(1)`. Show legend for line colour and tip size.
#' @param branch.length `character(1)`. If `"none"`, draw a cladogram.
#'   Passed through to [ggtree::ggtree()].
#' @param min.overlap `integer(1)`. Minimum overlap in columns for the
#'   pairwise alignment to count. Forwarded to [compare_motifs_lite()].
#'   Default `6`.
#' @param tryRC `logical(1)`. Also test reverse-complement alignments.
#'   Forwarded to [compare_motifs_lite()]. Default `TRUE`.
#' @param nthreads `numeric(1)`. Threads passed to [compare_motifs_lite()].
#'   `nthreads = 0` uses all available threads.
#' @param progress `logical(1)`. Print progress messages.
#' @param tiplogo `logical(1)`. Draw a sequence logo at each tip using
#'   [geom_logo()]. On by default. Only drawn on the `"rectangular"` and
#'   `"slanted"` layouts (other layouts warn and skip), and not available
#'   for `dist` input. Default `TRUE`.
#' @param tiplogo.align `logical(1)`. Align the tip logos to a common,
#'   zero-padded column frame (via [view_motifs_lite()]) so that homologous
#'   positions line up across tips. DNA/RNA only; falls back to unaligned
#'   per-tip logos otherwise. Default `TRUE`.
#' @param tiplogo.width `numeric(1)`. Width of one motif position, in tree
#'   x-axis units, for the tip logos. Default `1`.
#' @param tiplogo.height `numeric(1)`. Total height of a tip logo, in tree
#'   y-axis units, centred on each tip. Default `0.9`.
#' @param tiplogo.offset `numeric(1)` or `NULL`. Gap between the tips and the
#'   start of the logos, in tree x-axis units. If `NULL`, a small fraction of
#'   the tree width is used. Default `NULL`.
#' @param ... Additional arguments passed through to
#'   [ggtree::ggtree()].
#'
#' @return A `ggplot` (`ggtree`) object.
#'
#' @details
#' The PCC-to-distance conversion `d = (1 - score) / 2` is symmetric
#' (the comparison matrix from [compare_motifs_lite()] is symmetric in
#' matrix mode), so `as.dist()` and `hclust()` accept it directly.
#' Anti-correlated motif pairs land at the maximum distance of `1`;
#' identical motifs land at `0`.
#'
#' For more control over tree construction, run [compare_motifs_lite()]
#' directly with `matrix.out = "score"`, convert to a `dist` object
#' yourself, and pass that `dist` to `motif_tree_lite()`. The function
#' detects the `dist` input and skips the comparison step.
#'
#' @examples
#' \dontrun{
#' library(universalmotif)
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' if (requireNamespace("ggtree", quietly = TRUE)) {
#'   motif_tree_lite(jaspar, linecol = "none", labels = "name",
#'               layout = "rectangular")
#' }
#'
#' ## Equivalent two-step form, useful when you want to inspect the
#' ## distance matrix before drawing the tree:
#' if (requireNamespace("ggtree", quietly = TRUE)) {
#'   score.mat <- compare_motifs_lite(jaspar, matrix.out = "score")
#'   d <- as.dist((1 - score.mat) / 2)
#'   motif_tree_lite(d, labels = "name")
#' }
#' }
#'
#' @seealso [motif_tree()], [compare_motifs_lite()], [merge_similar_lite()],
#'   [geom_logo()], [view_motifs_lite()], [ggtree::ggtree()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @family lite motif functions
#' @export
motif_tree_lite <- function(motifs, layout = "rectangular", linecol = "family",
                        labels = "none", tipsize = "none", legend = TRUE,
                        branch.length = "none", min.overlap = 6,
                        tryRC = TRUE, nthreads = 1, progress = FALSE,
                        tiplogo = TRUE, tiplogo.align = TRUE,
                        tiplogo.width = 1, tiplogo.height = 0.9,
                        tiplogo.offset = NULL, ...) {

  ## --- arg validation ---------------------------------------------------
  if (!layout %in% c("rectangular", "slanted", "fan", "circular", "radial",
                     "equal_angle", "daylight"))
    stop("`layout` must be one of `rectangular`, `slanted`, `fan`, ",
         "`circular`, `radial`, `equal_angle`, or `daylight`", call. = FALSE)
  if (!is.character(linecol) || length(linecol) != 1L)
    stop("`linecol` must be a single string", call. = FALSE)
  if (!is.character(labels) || length(labels) != 1L)
    stop("`labels` must be a single string", call. = FALSE)
  if (!is.character(tipsize) || length(tipsize) != 1L)
    stop("`tipsize` must be a single string", call. = FALSE)
  if (!isTRUEorFALSE(legend))
    stop("`legend` must be a single logical", call. = FALSE)
  if (!is.character(branch.length) || length(branch.length) != 1L)
    stop("`branch.length` must be a single string", call. = FALSE)
  if (!is.numeric(min.overlap) || length(min.overlap) != 1L ||
      min.overlap < 1L)
    stop("`min.overlap` must be a positive integer", call. = FALSE)
  if (!isTRUEorFALSE(tryRC))
    stop("`tryRC` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(progress))
    stop("`progress` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(tiplogo))
    stop("`tiplogo` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(tiplogo.align))
    stop("`tiplogo.align` must be a single logical", call. = FALSE)

  nthreads <- resolve_nthreads(nthreads)

  ## --- distance matrix --------------------------------------------------
  if (is(motifs, "dist")) {

    tree <- ape_fun(ape::as.phylo(hclust(as.dist(motifs))))
    mot_names <- attr(motifs, "Labels")
    if (labels == "name") {
      tree$tip.label <- mot_names
    } else if (labels != "none") {
      warning("Trees from 'dist' objects can only use 'name' labels")
    }

  } else if (is.list(motifs)) {

    motifs <- convert_motifs(motifs)
    if (progress) message("Comparing motifs...")
    score.mat <- compare_motifs_lite(motifs,
                                 matrix.out  = "score",
                                 RC          = tryRC,
                                 min.overlap = min.overlap,
                                 nthreads    = nthreads)
    if (anyNA(score.mat))
      stop(wmsg("Found NA values in comparison matrix; check the inputs ",
                "or relax `min.overlap`"))
    dist.mat <- (1 - score.mat) / 2
    if (progress) message("Constructing phylogeny...")
    tree <- ape_fun(ape::as.phylo(hclust(as.dist(dist.mat))))
    if (labels != "none") {
      mot_names <- sapply(motifs, function(x) x[labels])
      tree$tip.label <- mot_names
    } else {
      mot_names <- vapply(motifs, function(x) x@name, character(1))
    }

  } else {
    stop("Input must be a 'dist' object or a 'list' of motifs", call. = FALSE)
  }

  if (is(motifs, "dist")) {
    if (linecol != "none") warning("'linecol' is not available for 'dist' objects")
    if (tipsize != "none") warning("'tipsize' is not available for 'dist' objects")
  }

  if (progress) message("Building tree...")

  ## Add tip logos (when enabled and the layout supports them) just before
  ## returning, regardless of which branch built the plot.
  finish <- function(p) add_tip_logos(p, motifs, mot_names, layout, tiplogo,
    tiplogo.align, tiplogo.width, tiplogo.height, tiplogo.offset, tryRC,
    min.overlap, nthreads)

  ## --- linecol-coloured branches ---------------------------------------
  if (linecol != "none" && !is(motifs, "dist")) {

    anno_list <- list()
    anno_bycol <- sapply(motifs, function(x) x[linecol])
    anno_unique <- unique(anno_bycol)
    anno_names <- mot_names
    for (i in seq_along(anno_unique)) {
      anno_list <- c(anno_list, list(anno_names[anno_bycol %in% anno_unique[i]]))
    }
    names(anno_list) <- anno_unique

    tree <- ggtree_fun(ggtree::groupOTU(tree, anno_list))

    if (labels != "none") {
      if (layout %in% c("rectangular", "slanted")) {
        p <- ggtree_fun({
          ggtree::ggtree(tree, aes(color = .data$group), layout = layout,
            branch.length = branch.length, ...) +
            ggtree::geom_tiplab(align = TRUE, linesize = 0.5)
        })
      } else {
        p <- ggtree_fun({
          ggtree::ggtree(tree, aes(color = .data$group), layout = layout,
            branch.length = branch.length, ...) +
            ggtree::geom_tiplab2(align = TRUE, linesize = 0.5)
        })
      }
    } else {
      p <- ggtree_fun({
        ggtree::ggtree(tree, aes(color = .data$group), layout = layout,
          branch.length = branch.length, ...)
      })
    }

    if (tipsize != "none") {
      anno_names <- mot_names
      anno_df <- data.frame(label = anno_names,
                            icscore = sapply(motifs, function(x) x[tipsize]))
      if (tipsize %in% c("pval", "qval", "eval")) {
        anno_df$icscore <- -log10(anno_df$icscore)
      }
      p <- ggtree_fun({
        ggtree::`%<+%`(p, anno_df) +
          ggtree::geom_tippoint(aes(size = .data$icscore))
      })
    }

    if (legend) {
      return(finish(p) +
        theme(legend.position = "right", legend.title = element_blank()))
    } else return(finish(p))

  }

  ## --- plain tree (no per-clade colouring) -----------------------------
  if (labels != "none") {

    if (layout %in% c("rectangular", "slanted")) {
      p <- ggtree_fun({
        ggtree::ggtree(tree, layout = layout, branch.length = branch.length, ...) +
          ggtree::geom_tiplab(align = TRUE, linesize = 0.5)
      })
    } else {
      p <- ggtree_fun({
        ggtree::ggtree(tree, layout = layout, branch.length = branch.length, ...) +
          ggtree::geom_tiplab2(align = TRUE, linesize = 0.5)
      })
    }

    if (tipsize != "none" && !is(motifs, "dist")) {
      anno_names <- mot_names
      anno_df <- data.frame(label = anno_names,
                            icscore = sapply(motifs, function(x) x[tipsize]))
      if (tipsize %in% c("pval", "qval", "eval")) {
        anno_df$icscore <- -log10(anno_df$icscore)
      }
      p <- ggtree_fun({
        ggtree::`%<+%`(p, anno_df) +
          ggtree::geom_tippoint(aes(size = .data$icscore))
      })
    }
    return(finish(p))

  } else {

    p <- ggtree_fun(
      ggtree::ggtree(tree, layout = layout, branch.length = branch.length, ...)
    )
    if (tipsize != "none" && !is(motifs, "dist")) {
      anno_names <- mot_names
      anno_df <- data.frame(label = anno_names,
                            icscore = sapply(motifs, function(x) x[tipsize]))
      if (tipsize %in% c("pval", "qval", "eval")) {
        anno_df$icscore <- -log10(anno_df$icscore)
      }
      p <- ggtree_fun({
        ggtree::`%<+%`(p, anno_df) +
          ggtree::geom_tippoint(aes(size = .data$icscore))
      })
    }
    return(finish(p))

  }

}

#-----------------------------------------------------------------------------
# Add per-tip sequence logos to a ggtree plot built by motif_tree_lite(). Returns
# `p` unchanged when tip logos are off, the layout cannot show them, or the
# input was a `dist` object (no matrices to draw).

add_tip_logos <- function(p, motifs, mot_names, layout, tiplogo, tiplogo.align,
  tiplogo.width, tiplogo.height, tiplogo.offset, tryRC, min.overlap,
  nthreads) {

  if (!tiplogo) return(p)

  if (is(motifs, "dist")) {
    warning("tip logos are not available for `dist` input; skipped",
      call. = FALSE)
    return(p)
  }
  if (!layout %in% c("rectangular", "slanted")) {
    warning("tip logos require a `rectangular` or `slanted` layout; skipped",
      call. = FALSE)
    return(p)
  }

  td <- p$data
  td <- td[td$isTip, c("label", "x", "y"), drop = FALSE]
  if (!nrow(td)) return(p)

  alphabet <- motifs[[1]]@alphabet
  cs <- switch(alphabet, "DNA" = DNA_COLOURS, "RNA" = RNA_COLOURS,
    "AA" = AA_COLOURS, NULL)

  if (tiplogo.align && !alphabet %in% c("DNA", "RNA")) {
    warning("tip logo alignment requires DNA/RNA motifs; drawing unaligned ",
      "logos", call. = FALSE)
    tiplogo.align <- FALSE
  }

  if (tiplogo.align) {
    ## dedup = TRUE: the logo names are overwritten with the tip labels just
    ## below, so disambiguating any duplicate motif names here is harmless.
    logos <- align_motif_mats(motifs, use.type = "ICM", tryRC = tryRC,
      min.overlap = min.overlap, nthreads = nthreads, dedup = TRUE)
  } else {
    logos <- lapply(motifs, function(m) convert_type(m, "ICM")@motif)
  }
  names(logos) <- mot_names

  if (is.null(tiplogo.offset)) {
    rng <- max(td$x) - min(td$x)
    tiplogo.offset <- if (rng > 0) rng * 0.05 else 0.5
  }
  td$x0 <- max(td$x) + tiplogo.offset

  ## The tip logos are ICM sequence logos, so within each position the
  ## letters should be stacked tallest-on-top (the standard logo
  ## convention). geom_logo() takes bare matrices and defaults to leaving
  ## the rows in their given order, so ask it to sort the letters by height
  ## here instead.
  p + geom_logo(aes(x = .data$x0, y = .data$y, motif = .data$label),
    data = td, logo = logos, inherit.aes = FALSE, width = tiplogo.width,
    height = tiplogo.height, colour.scheme = cs, sort.positions = TRUE)

}
