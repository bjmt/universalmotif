######################################################################
## Benjamin Tremblay
##
## Wrappers for functions from other packages
##
######################################################################

#' @export
motifStack <- function(pfms, layout = c("stack", "treeview",
                                       "phylog", "radialPhylog")) {
  motifs <- lapply(pfms, function(x) convert_motifs(x, "motifStack-pfm"))
  motifStack::motifStack(pfms = motifs, layout = layout)
}

#' @export
seqLogo <- function(motif, ic.scale = TRUE, xaxis = TRUE, yaxis = TRUE,
                    xfontsize = 15, yfontsize = 15) {
  motif <- convert_motifs(motif, "seqLogo-pwm")
  seqLogo::seqLogo(pwm = motif, ic.scale = ic.scale, xaxis = xaxis,
                   yaxis = yaxis, xfontsize = xfontsize,
                   yfontsize = yfontsize)
}
