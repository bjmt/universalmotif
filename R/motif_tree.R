#' Generate ggplot2 motif trees with ggtree.
#'
#' For more powerful motif tree functions, see the motifStack package.
#'
#' @param motifs \code{list}, \code{dist} See \code{\link{convert_motifs}} for
#'    available formats. 
#' @param layout \code{character(1)} One of \code{c('rectangular', 'slanted', 'fan', 'circular',
#'    'radial', 'equal_angle', 'daylight')}.
#' @param linecol \code{character(1)} \linkS4class{universalmotif} slot to use to
#'    colour lines (e.g. 'family').
#' @param labels \code{character(1)} \linkS4class{universalmotif} slot to use to label
#'    tips (e.g. 'name').
#' @param tipsize \code{character(1)} \linkS4class{universalmotif} slot to use to
#'    control tip size (e.g. 'icscore').
#' @param legend \code{logical(1)} Show legend for line colour and tip size.
#' @param branch.length \code{character(1)} If 'none', draw a cladogram.
#' @param db.scores \code{data.frame} See \code{\link{compare_motifs}}.
#' @param method \code{character(1)} One of \code{'Pearson', 'Euclidean', 'KL'}.
#' @param use.type \code{character(1)} One of \code{'PCM'} (Pearson only),
#'    \code{'PPM'} (any method), \code{'PWM'} (Pearson only), and \code{'ICM'}
#'    (any method). The two allow for taking into account the background
#'    frequencies (for ICM, only if \code{relative_entropy = TRUE}).
#' @param min.overlap \code{numeric(1)} Minimum overlap required when aligning the
#'    motifs. Setting this to a number higher then the width of the motifs
#'    will not allow any overhangs. Can also be a number less than 1,
#'    representing the minimum fraction that the motifs must overlap.
#' @param tryRC \code{logical(1)} Try the reverse complement of the motifs as well,
#'    report the best score.
#' @param min.mean.ic \code{numeric(1)} Minimum information content between the
#'    two motifs for an alignment to be scored. This helps prevent scoring
#'    alignments between low information content regions of two motifs.
#' @param relative_entropy \code{logical(1)} For ICM calculation. See
#'    \code{\link{convert_type}}.
#' @param ... \pkg{ggtree} params.
#'
#' @return \code{ggplot} object.
#'
#' @details
#'    See \code{\link[ggtree]{ggtree}} for more detailed descriptions of
#'    parameters.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.tree <- motif_tree(jaspar, linecol = "none", labels = "name",
#'                           layout = "rectangular")
#'
#' @references
#'    \insertRef{ggplot2}{universalmotif}
#'
#'    \insertRef{ggtree}{universalmotif}
#'
#' @seealso \code{\link[motifStack]{motifStack}}, \code{\link{compare_motifs}},
#'    \code{\link[ggtree]{ggtree}}
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_tree <- function(motifs, layout = "circular", linecol = "family",
                       labels = "none", tipsize = "none", legend = TRUE,
                       branch.length = "none", db.scores, method = "MPCC",
                       use.type = "PPM", min.overlap = 6, tryRC = TRUE,
                       min.mean.ic = 0.5, relative_entropy = FALSE, ...){

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(layout = args$layout, linecol = args$linecol,
                                      labels = args$labels, tipesize = args$tipsize,
                                      branch.length = args$branch.length,
                                      method = args$method, use.type = args$use.type),
                                 numeric(), logical(), "character")
  num_check <- check_fun_params(list(min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic),
                                numeric(), logical(), "numeric")
  logi_check <- check_fun_params(list(legend = args$legend, tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy),
                                 numeric(), logical(), "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (is(motifs, "dist")) {
    tree <- ape::as.phylo(hclust(motifs))
  } else {
    motifs <- convert_motifs(motifs)
    if (!missing(db.scores)) {
      tree <- compare_motifs(motifs, db.scores = db.scores,
                             use.type = use.type,
                             method = method, tryRC = tryRC,
                             min.overlap = min.overlap,
                             min.mean.ic = min.mean.ic,
                             relative_entropy = relative_entropy)
    } else {
      tree <- compare_motifs(motifs,
                             use.type = use.type,
                             method = method, tryRC = tryRC,
                             min.overlap = min.overlap,
                             min.mean.ic = min.mean.ic,
                             relative_entropy = relative_entropy)
    }
    if (method %in% c("MPCC", "PCC", "SW", "MSW")) tree <- 1 / tree  # !!!!
    tree <- ape::as.phylo(hclust(as.dist(tree)))
  }

  if (labels != "none") {
    mot_names <- sapply(motifs, function(x) x[labels])
    tree$tip.label <- mot_names
  } else {
    mot_names <- vapply(motifs, function(x) x["name"], character(1))
  }

  if (linecol != "none") {

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
      anno_df <- data.frame(name = anno_names,
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

    if (tipsize != "none") {
      anno_names <- mot_names
      anno_df <- data.frame(name = anno_names,
                            icscore = sapply(motifs, function(x) x[tipsize]))
      if (tipsize %in% c("pval", "qval", "eval")) {
        anno_df$icscore <- -log10(anno_df$icscore)
      }
      p <- p %<+% anno_df + geom_tippoint(aes(size = icscore))
    }
    return(p)

  } else {

    p <- ggtree(tree, layout = layout, branch.length = branch.length, ...)
    if (tipsize != "none") {
      anno_names <- mot_names
      anno_df <- data.frame(name = anno_names,
                            icscore = sapply(motifs, function(x) x[tipsize]))
      if (tipsize %in% c("pval", "qval", "eval")) {
        anno_df$icscore <- -log10(anno_df$icscore)
      }
      p <- p %<+% anno_df + geom_tippoint(aes(size = icscore))
    }
    return (p)

  }
 
} 

# grid.arrange(p1, p2 + scale_x_reverse(), nrow = 1) (package=egg)
