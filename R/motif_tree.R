#' Generate ggplot2 motif trees with ggtree.
#'
#' For more powerful motif tree functions, see the motifStack package.
#'
#' @param motifs List of \linkS4class{universalmotif} objects or 'dist' object.
#' @param layout Character. One of 'rectangular', 'slanted', 'fan', 'circular',
#'    'radial', 'equal_angle', and 'daylight'.
#' @param linecol Character. \linkS4class{universalmotif} slot to use to
#'    colour lines (e.g. 'family').
#' @param labels Character. \linkS4class{universalmotif} slot to use to label
#'    tips (e.g. 'name').
#' @param tipsize Character. \linkS4class{universalmotif} slot to use to
#'    control tip size (e.g. 'icscore').
#' @param legend Logical. Show legend for line colour and tip size.
#' @param branch.length Character. If 'none', draw a cladogram.
#' @param db.scores data.frame.
#' @param method Character.
#' @param min.overlap Numeric.
#' @param tryRC Logical.
#' @param min.mean.ic Numeric.
#' @param relative_entropy Logical.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#' @param ... ggtree params.
#'
#' @return ggplot2 object.
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
#'    \insertRef{ggtree}{universalmotif}
#'
#'    \insertRef{ggplot2}{universalmotif}

#'
#' @seealso \code{\link[motifStack]{motifStack}}, \code{\link{compare_motifs}},
#'    \code{\link[ggtree]{ggtree}}
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_tree <- function(motifs, layout = "circular", linecol = "family",
                       labels = "none", tipsize = "none", legend = TRUE,
                       branch.length = "none", db.scores, method = "Pearson",
                       min.overlap = 6, tryRC = TRUE, min.mean.ic = 0.5,
                       relative_entropy = TRUE, BPPARAM = SerialParam(), ...){

  if (class(motifs) == "dist") {
    tree <- as.phylo(hclust(motifs))
  } else {
    motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
    tree <- as.phylo(hclust(as.dist(compare_motifs(motifs, db.scores = db.scores,
                                                   method = method, tryRC = tryRC,
                                                   min.overlap = min.overlap,
                                                   min.mean.ic = min.mean.ic,
                                                   relative_entropy = relative_entropy,
                                                   BPPARAM = BPPARAM))))
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

    tree <- groupOTU(tree, anno_list)

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
