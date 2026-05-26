#' Generate ggtree-based motif trees using `compare_motifs2()` (v2).
#'
#' `motif_tree2()` is the leaner counterpart of [motif_tree()]: it builds the
#' distance matrix via [compare_motifs2()] (mean Pearson correlation,
#' with built-in p-value / q-value machinery) instead of
#' [compare_motifs()] (older similarity / distance metrics such as
#' `EUCL` / `WEUCL`). The distance for `hclust()` is derived from the
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
#'   result of `as.dist((1 - compare_motifs2(...)) / 2)`) skips the
#'   comparison step entirely and goes straight to tree construction.
#' @param layout `character(1)`. One of `c('rectangular', 'slanted',
#'   'fan', 'circular', 'radial', 'equal_angle', 'daylight')`. Passed
#'   through to [ggtree::ggtree()].
#' @param linecol `character(1)`. Motif slot used to colour the
#'   branches (e.g. `"family"`). Not used when `motifs` is a `dist`.
#' @param labels `character(1)`. Motif slot used to label the tips
#'   (e.g. `"name"`). For `dist` input, only `"name"` is meaningful.
#' @param tipsize `character(1)`. Motif slot used to size the tips
#'   (e.g. `"icscore"`). Not used when `motifs` is a `dist`.
#' @param legend `logical(1)`. Show legend for line colour and tip size.
#' @param branch.length `character(1)`. If `"none"`, draw a cladogram.
#'   Passed through to [ggtree::ggtree()].
#' @param min.overlap `integer(1)`. Minimum overlap in columns for the
#'   pairwise alignment to count. Forwarded to [compare_motifs2()].
#'   Default `6`.
#' @param tryRC `logical(1)`. Also test reverse-complement alignments.
#'   Forwarded to [compare_motifs2()]. Default `TRUE`.
#' @param nthreads `numeric(1)`. Threads passed to [compare_motifs2()].
#'   `nthreads = 0` uses all available threads.
#' @param progress `logical(1)`. Print progress messages.
#' @param ... Additional arguments passed through to
#'   [ggtree::ggtree()].
#'
#' @return A `ggplot` (`ggtree`) object.
#'
#' @details
#' The PCC-to-distance conversion `d = (1 - score) / 2` is symmetric
#' (the comparison matrix from [compare_motifs2()] is symmetric in
#' matrix mode), so `as.dist()` and `hclust()` accept it directly.
#' Anti-correlated motif pairs land at the maximum distance of `1`;
#' identical motifs land at `0`.
#'
#' For more control over tree construction, run [compare_motifs2()]
#' directly with `matrix.out = "score"`, convert to a `dist` object
#' yourself, and pass that `dist` to `motif_tree2()`. The function
#' detects the `dist` input and skips the comparison step.
#'
#' @examples
#' \dontrun{
#' library(universalmotif)
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' if (requireNamespace("ggtree", quietly = TRUE)) {
#'   motif_tree2(jaspar, linecol = "none", labels = "name",
#'               layout = "rectangular")
#' }
#'
#' ## Equivalent two-step form, useful when you want to inspect the
#' ## distance matrix before drawing the tree:
#' if (requireNamespace("ggtree", quietly = TRUE)) {
#'   score.mat <- compare_motifs2(jaspar, matrix.out = "score")
#'   d <- as.dist((1 - score.mat) / 2)
#'   motif_tree2(d, labels = "name")
#' }
#' }
#'
#' @seealso [motif_tree()], [compare_motifs2()], [merge_similar2()],
#'   [ggtree::ggtree()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
motif_tree2 <- function(motifs, layout = "circular", linecol = "family",
                        labels = "none", tipsize = "none", legend = TRUE,
                        branch.length = "none", min.overlap = 6,
                        tryRC = TRUE, nthreads = 1, progress = FALSE, ...) {

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
    score.mat <- compare_motifs2(motifs,
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
      return(p + theme(legend.position = "right", legend.title = element_blank()))
    } else return(p)

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
    return(p)

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
    return(p)

  }

}
