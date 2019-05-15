#' Generate \pkg{ggplot2} motif trees with \pkg{ggtree}.
#'
#' For more powerful motif tree functions, see the \pkg{motifStack} package.
#' The [motif_tree()] function compares motifs with [compare_motifs()] to create
#' a distance matrix, which is used to generate a phylogeny via \pkg{ape}.
#' This can be plotted with [ggtree::ggtree()]. The purpose of this function
#' is simply to combine the [compare_motifs()] and [ggtree::ggtree()] steps
#' into one. For more control over tree creation, it is recommend to do these
#' steps separately. See the Advanced Usage vignette for such a workthrough.
#'
#' @param motifs `list`, `dist` See [convert_motifs()] for
#'    available formats. Alternatively, the resulting comparison matrix from
#'    [compare_motifs()] (run `as.dist(results)` beforehand; if the comparison was
#'    performed with a similarity metric, make sure to convert to distances first).
#' @param layout `character(1)` One of `c('rectangular', 'slanted', 'fan', 'circular',
#'    'radial', 'equal_angle', 'daylight')`. See [ggtree::ggtree()].
#' @param linecol `character(1)` [universalmotif-class] slot to use to
#'    colour lines (e.g. 'family'). Not available for `dist` input (see examples
#'    for how to add it manually). See [ggtree::ggtree()].
#' @param labels `character(1)` [universalmotif-class] slot to use to label
#'    tips (e.g. 'name'). For `dist` input, only 'name' is available.
#'    See [ggtree::ggtree()].
#' @param tipsize `character(1)` [universalmotif-class] slot to use to
#'    control tip size (e.g. 'icscore'). Not available for `dist` input (see
#'    examples for how to add it manually). See [ggtree::ggtree()].
#' @param legend `logical(1)` Show legend for line colour and tip size.
#'    See [ggtree::ggtree()].
#' @param branch.length `character(1)` If 'none', draw a cladogram.
#'    See [ggtree::ggtree()].
#' @param db.scores `data.frame` See [compare_motifs()].
#' @param method `character(1)` One of `c('PCC', 'MPCC', 'EUCL', 'SW', 'MSW',
#'    'KL', 'MKL')`. See [compare_motifs()].
#' @param use.type `character(1)`c('PPM', 'ICM')`. The latter allows for taking
#'    into account the background
#'    frequencies (only if `relative_entropy = TRUE`). See [compare_motifs()].
#' @param min.overlap `numeric(1)` Minimum overlap required when aligning the
#'    motifs. Setting this to a number higher then the width of the motifs
#'    will not allow any overhangs. Can also be a number less than 1,
#'    representing the minimum fraction that the motifs must overlap. See
#'    [compare_motifs()].
#' @param tryRC `logical(1)` Try the reverse complement of the motifs as well,
#'    report the best score. See [compare_motifs()].
#' @param min.mean.ic `numeric(1)` Minimum information content between the
#'    two motifs for an alignment to be scored. This helps prevent scoring
#'    alignments between low information content regions of two motifs. See
#'    [compare_motifs()].
#' @param min.position.ic `numeric(1)` Minimum information content required between
#'    individual alignment positions for it to be counted in the final alignment
#'    score. It is recommended to use this together with `normalise.scores = TRUE`,
#'    as this will help punish scores resulting from only a fraction of an
#'    alignment.
#' @param relative_entropy `logical(1)` For ICM calculation. See
#'    [convert_type()].
#' @param progress `logical(1)` Deprecated. Does nothing.
#' @param BP `logical(1)` Deprecated. See `nthreads`.
#' @param nthreads `numeric(1)` Run [compare_motifs()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#' @param ... \pkg{ggtree} params. See [ggtree::ggtree()].
#'
#' @return ggplot object.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.tree <- motif_tree(jaspar, linecol = "none", labels = "name",
#'                           layout = "rectangular")
#'
#' \dontrun{
#' ## When inputting a dist object, the linecol and tipsize options are
#' ## not available. To add these manually:
#'
#' library(MotifDb)
#' library(ggtree)
#' library(ggplot2)
#'
#' motifs <- filter_motifs(MotifDb, organism = "Athaliana")[1:50]
#' comparison <- compare_motifs(motifs, method = "MPCC")
#' comparison <- as.dist(1 - comparison)
#' mot.names <- attr(comparison, "Labels")
#' tree <- motif_tree(comparison)
#'
#' annotations <- data.frame(label = mot.names,
#'                           icscore = sapply(motifs, function(x) x["icscore"]),
#'                           family = sapply(motifs, function(x) x["family"]))
#'
#' tree <- tree %<+% annotations +
#'           geom_tippoint(aes(size = icscore)) +
#'           aes(colour = family) +
#'           theme(legend.position = "right",
#'                 legend.title = element_blank())
#' }
#'
#' @references
#'    \insertRef{ggplot2}{universalmotif}
#'
#'    \insertRef{ggtree}{universalmotif}
#'
#' @seealso [motifStack::motifStack()], [compare_motifs()],
#'    [ggtree::ggtree()], [ggplot2::ggplot()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_tree <- function(motifs, layout = "circular", linecol = "family",
                       labels = "none", tipsize = "none", legend = TRUE,
                       branch.length = "none", db.scores, method = "MEUCL",
                       use.type = "PPM", min.overlap = 6,
                       min.position.ic = 0, tryRC = TRUE,
                       min.mean.ic = 0.25, relative_entropy = FALSE,
                       progress = FALSE, BP = FALSE, nthreads = 1, ...) {

  # TODO: allow for user-provided linecol, labels, tipsize instead of just
  #       pulling from motif slots.

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!layout %in% c("rectangular", "slanted", "fan", "circular", "radial",
                     "equal_angle", "daylight")) {
    layout_check <- paste0(" * Incorrect 'layout': expected `rectangular`, ",
                           "`slanted`, `fan`, `circular`, `radial`, `equal_angle`",
                           " or `daylight`; got `", layout, "`")
    layout_check <- wmsg2(layout_check, 4, 2)
    all_checks <- c(all_checks, layout_check)
  }
  if (!method %in% COMPARE_METRICS) {
    method_check <- paste0(" * Incorrect 'method': expected `PCC`, `MPCC`,",
                           " `EUCL`, `MEUCL`, `SW`, `MSW`, `KL` or `MKL`; got `",
                           method, "`")
    method_check <- wmsg2(method_check)
    all_checks <- c(all_checks, method_check)
  }
  if (!use.type %in% c("PPM", "ICM")) {
    use.type_check <- paste0(" * Incorrect 'use.type': expected `PPM` or `ICM`; got `",
                             use.type, "`")
    use.type_check <- wmsg2(use.type_check, 4, 2)
    all_checks <- c(all_checks, use.type_check)
  }
  char_check <- check_fun_params(list(layout = args$layout, linecol = args$linecol,
                                      labels = args$labels, tipesize = args$tipsize,
                                      branch.length = args$branch.length,
                                      method = args$method, use.type = args$use.type),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     min.position.ic = args$min.position.ic),
                                numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(legend = args$legend, tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      progress = args$progress, BP = args$BP),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(all_checks, char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (method %in% c("PCC", "SW"))
    stop("'PCC', 'SW' are not allowed; please use 'MPCC' and 'MSW' instead")

  if (is(motifs, "dist")) {
    tree <- ape::as.phylo(hclust(motifs))
    mot_names <- attr(motifs, "Labels")
    if (labels == "name") {
      tree$tip.label <- mot_names
    } else if (labels != "none") {
      warning("Trees from 'dist' objects can only use 'name' labels")
    }
  } else if (is.list(motifs)) {
    motifs <- convert_motifs(motifs)
    if (progress && !BP)
      cat("Comparing motifs... ")
    else if (progress)
      cat("Comparing motifs...\n")
    tree <- compare_motifs(motifs,
                           use.type = use.type,
                           method = method, tryRC = tryRC,
                           min.overlap = min.overlap,
                           min.mean.ic = min.mean.ic,
                           relative_entropy = relative_entropy,
                           BP = BP, progress = progress)
    if (method == "MPCC") tree <- 1 - tree
    else if (method == "MSW") tree <- 2 - tree
    if (progress) cat("Constructing phylogeny...\n")
    tree <- ape::as.phylo(hclust(as.dist(tree)))
    if (labels != "none") {
      mot_names <- sapply(motifs, function(x) x[labels])
      tree$tip.label <- mot_names
    } else {
      mot_names <- vapply(motifs, function(x) x@name, character(1))
    }
  } else {
    stop("Input must be a 'dist' object or a 'list' of motifs")
  }

  if (is(motifs, "dist")) {
    if (linecol != "none") warning("'linecol' is not available for 'dist' objects")
    if (tipsize != "none") warning("'tipsize' is not available for 'dist' objects")
  }

  if (progress) cat("Building tree...\n")

  if (linecol != "none" && !is(motifs, "dist")) {

    anno_list <- list()
    anno_bycol <- sapply(motifs, function(x) x[linecol])
    anno_unique <- unique(anno_bycol)
    anno_names <- mot_names
    for (i in seq_along(anno_unique)) {
      anno_list <- c(anno_list, list(anno_names[anno_bycol %in% anno_unique[i]]))
    }
    names(anno_list) <- anno_unique

    tree <- ggtree::groupOTU(tree, anno_list)

    if (labels != "none") {

      if (layout %in% c("rectangular", "slanted")) {
        p <- ggtree(tree, aes(color = group), layout = layout,
                    branch.length = branch.length, ...) +
               geom_tiplab(align = TRUE, linesize = 0.5)
      } else {
        p <- ggtree(tree, aes(color = group), layout = layout,
                    branch.length = branch.length, ...) +
               geom_tiplab2(align = TRUE, linesize = 0.5)
      }

    } else {
      p <- ggtree(tree, aes(color = group), layout = layout,
                  branch.length = branch.length, ...)
    }

    if (tipsize != "none") {
      anno_names <- mot_names
      anno_df <- data.frame(label = anno_names,
                            icscore = sapply(motifs, function(x) x[tipsize]))
      if (tipsize %in% c("pval", "qval", "eval")) {
        anno_df$icscore <- -log10(anno_df$icscore)
      }
      p <- p %<+% anno_df + geom_tippoint(aes(size = icscore))
    }

    if (legend) {
      return(p + theme(legend.position = "right", legend.title = element_blank()))
    } else return(p)

  }

  if (labels != "none") {

    if (layout %in% c("rectangular", "slanted")) {
      p <- ggtree(tree, layout = layout, branch.length = branch.length, ...) +
             geom_tiplab(align = TRUE, linesize = 0.5)
    } else {
      p <- ggtree(tree, layout = layout, branch.length = branch.length, ...) +
             geom_tiplab2(align = TRUE, linesize = 0.5)
    }

    if (tipsize != "none" && !is(motifs, "dist")) {
      anno_names <- mot_names
      anno_df <- data.frame(label = anno_names,
                            icscore = sapply(motifs, function(x) x[tipsize]))
      if (tipsize %in% c("pval", "qval", "eval")) {
        anno_df$icscore <- -log10(anno_df$icscore)
      }
      p <- p %<+% anno_df + geom_tippoint(aes(size = icscore))
    }
    return(p)

  } else {

    p <- ggtree(tree, layout = layout, branch.length = branch.length, ...)
    if (tipsize != "none" && !is(motifs, "dist")) {
      anno_names <- mot_names
      anno_df <- data.frame(label = anno_names,
                            icscore = sapply(motifs, function(x) x[tipsize]))
      if (tipsize %in% c("pval", "qval", "eval") && !is(motifs, "dist")) {
        anno_df$icscore <- -log10(anno_df$icscore)
      }
      p <- p %<+% anno_df + geom_tippoint(aes(size = icscore))
    }
    return (p)

  }

}

# grid.arrange(p1, p2 + scale_x_reverse(), nrow = 1) (package=egg)
