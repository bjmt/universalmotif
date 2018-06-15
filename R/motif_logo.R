#' Generate ggplot2 motif logos.
#'
#' Create logos using ggplot2 + ggseqlogo, seqLogo, motifStack or Logolas.
#'
#' @param motif universalmotif object.
#' @param engine Character. ggseqlogo, seqLogo, motifStack, Logolas.
#' @param type Character. For \code{engine = c('ggseqlogo', 'seqLogo',
#'             'motifStack')}: 'ic' or 'prob'. For \code{engine = Logolas}:
#'             'Logo' or 'EDLogo'.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
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
motif_logo <- function(motif, engine = "ggseqlogo", type = "ic",
                       BPPARAM = bpparam(), ...) {

  motif <- convert_motifs(motif, BPPARAM = BPPARAM)

  if (engine == "ggseqlogo") {

    if (type == "ic") type <- "bits"

    motif <- convert_type(motif, "PPM", BPPARAM = BPPARAM)

    ggseqlogo(motif["motif"], method = type, ...) +
      xlab("Position") +
      ylab(c("Probability", "Information content")[c("prob", "bits") %in% type]) +
      # theme(axis.ticks.x = element_line(colour = "gray", size = 0.3),
            # axis.ticks.y = element_line(colour = "gray", size = 0.3),
            # axis.ticks.length = unit(0.3, "cm"),
            # axis.line.x = element_line(colour = "gray", size = 0.3),
            # axis.line.y = element_line(colour = "gray", size = 0.3),
            # axis.text.x = element_text(margin = margin(0.2, unit = "cm"))) +
      scale_y_continuous(limits = c(0, c(1, 2)[c("prob", "bits") %in% type])) 

  } else if (engine == "seqLogo") {

    motif <- convert_motifs(motif, "seqLogo-pwm")
    seqLogo(motif, ic.scale = type, ...)

  } else if (engine == "motifStack") {

    if (type == "ic") type <- TRUE else if (type == "prob") type <- FALSE
    plotMotifLogo(convert_motifs(motif, "motifStack-pfm"), motif["name"],
                  p = motif["bkg"], ic.scale = type, ...)

  } else if (engine == "Logolas") {
  
    motif <- convert_type(motif, "PPM", BPPARAM = BPPARAM)

    if (type %in% c("ic", "prob", "Logo")) type <- "Logo"
    if (type == "EDLogo") type <- "EDLogo"

    logomaker(motif["motif"], type = type, ...)
  
  }

}
