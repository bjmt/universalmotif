######################################################################
## Benjamin Tremblay
##
## Convert motifs between different classes. convert_motifs works on
## one object at a time, so lapply must be used on lists.
##
######################################################################

#' @describeIn convert_motifs Convert from universalmotif.
#' @include utils.R universalmotif-class.R generics.R
setMethod("convert_motifs", signature = "universalmotif",
          definition = function(motif, out_class, pseudoweight = 0.8,
                                background = c("A" = 0.25, "C" = 0.25,
                                               "G" = 0.25, "T" = 0.25)) {
            if (out_class == "PFMatrix") {
              motif <- umot_to_pfmatrix(motif)
            }
            if (out_class == "PWMatrix") {
              if (length(motif@bkg) > 0) {
                background[1] <- motif@bkg[1]
                background[2] <- motif@bkg[2]
                background[3] <- motif@bkg[3]
                background[4] <- motif@bkg[4]
              }
              motif <- umot_to_pfmatrix(motif)
              motif <- TFBSTools::toPWM(motif, type = "log2probratio",
                                        pseudocounts = pseudoweight,
                                        bg = background)
            }
            return(motif)
          })

#' @describeIn convert_motifs Convert a list of motifs.
setMethod("convert_motifs", signature = "list",
          definition = function(motif, out_class, ...) {
            lapply(X = motif, FUN = convert_motifs, out_class, ...)
          })

######################################################################
######################################################################

umot_to_pfmatrix <- function(motif) {

  if (motif@type == "PPM") {
    motif <- convert_type(motif, "PCM")
  }
  if (motif@type == "PWM") {
    motif <- convert_type(motif, "PCM")
  }

  if (!motif@alphabet %in% c("DNA", "RNA")) {
    stop("PFMatrix can only be of type DNA or RNA")
  }
  bg <- c("A" = motif@bkg[1], "C" = motif@bkg[2],
          "G" = motif@bkg[3], "T" = motif@bkg[4])
  motif <- TFBSTools::PFMatrix(ID = motif@name, name = motif@name,
                               strand = motif@strand[1], bg = bg,
                               profileMatrix = motif@motif)
}
