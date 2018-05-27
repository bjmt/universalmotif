#' Generate ggplot2 motif logos.
#'
#' Create logos using ggplot2 + ggseqlogo or seqLogo.
#'
#' @param motif universalmotif object.
#' @param engine Character. ggseqlogo, seqLogo, motifStack.
#' @param type Character. ic or prob.
#' @param ... Package-specific params.
#'
#' @return For ggseqlogo: a ggplot2 object.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.logo <- motif_logo(jaspar[[1]], type = "prob")
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_logo <- function(motif, engine = "ggseqlogo", type = "ic", ...) {

  motif <- convert_motifs(motif)

  if (engine == "ggseqlogo") {

    if (type == "ic") type <- "bits"

    ggseqlogo(motif["motif"], method = type, ...) +
      xlab("Position") +
      ylab(c("Probability", "Information content")[c("prob", "bits") %in% type]) +
      theme(axis.ticks.x = element_line(colour = "gray", size = 0.3),
            axis.ticks.y = element_line(colour = "gray", size = 0.3),
            axis.ticks.length = unit(0.3, "cm"),
            axis.line.x = element_line(colour = "gray", size = 0.3),
            axis.line.y = element_line(colour = "gray", size = 0.3),
            axis.text.x = element_text(margin = margin(0.2, unit = "cm"))) +
      scale_y_continuous(limits = c(0, c(1, 2)[c("prob", "bits") %in% type])) 

  } else if (engine == "seqLogo") {

    motif <- convert_motifs(motif, "seqLogo-pwm")
    seqLogo(motif, ic.scale = type, ...)

  } else if (engine == "motifStack") {

    if (type == "ic") type <- TRUE else if (type == "prob") type <- FALSE
    plotMotifLogo(convert_motifs(motif, "motifStack-pfm"), motif["name"],
                  p = motif["bkg"], ic.scale = type, ...)

  } 

}
