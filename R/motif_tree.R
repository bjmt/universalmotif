#' Generate ggplot2 motif trees with ggtree.
#'
#' For more powerful motif tree functions, see motifStack.
#'
#' @param motifs List of motifs or 'dist' object.
#' @param layout Character.
#' @param linecol Character. (family, organism, etc)
#' @param labels Character. (name, altname, etc)
#' @param tipsize Character. (icscore, pval, qval, eval)
#' @param legend Logical.
#' @param branch.length Character.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#' @param ... ggtree params.
#'
#' @return ggplot2 object.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.tree <- motif_tree(jaspar, linecol = "none", labels = "name",
#'                           layout = "rectangular")
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_tree <- function(motifs, layout = "circular", linecol = "family",
                       labels = "none", tipsize = "none", legend = TRUE,
                       branch.length = "none", BPPARAM = bpparam(), ...){

  if (class(motifs) == "dist") {
    tree <- as.phylo(hclust(motifs))
  } else {
    motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
    tree <- as.phylo(hclust(as.dist(compare_motifs(motifs, BPPARAM = BPPARAM))))
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
