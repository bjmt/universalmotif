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

            motif2 <- NULL

            if (out_class == "TFBSTools-PFMatrix") {
              motif2 <- umot_to_pfmatrix(motif)
            }

            if (out_class == "TFBSTools-PWMatrix") {
              if (length(motif@bkg) > 0) {
                background[1] <- motif@bkg[1]
                background[2] <- motif@bkg[2]
                background[3] <- motif@bkg[3]
                background[4] <- motif@bkg[4]
              }
              motif2 <- umot_to_pfmatrix(motif)
              motif2 <- TFBSTools::toPWM(motif2, type = "log2probratio",
                                        pseudocounts = pseudoweight,
                                        bg = background)
            }

            if (out_class == "TFBSTools-ICMatrix") {
              if (length(motif@bkg) > 0) {
                background[1] <- motif@bkg[1]
                background[2] <- motif@bkg[2]
                background[3] <- motif@bkg[3]
                background[4] <- motif@bkg[4]
              }
              motif2 <- umot_to_pfmatrix(motif)
              motif2 <- TFBSTools::toICM(motif2, pseudocounts = pseudoweight,
                                        bg = background)
            }

            if (out_class == "seqLogo-pwm") {
              if (motif_slots(motif, "type") == "PCM") {
                motif2 <- convert_type(motif, "PPM")
                motif2 <- seqLogo::makePWM(motif_slots(motif2, "motif"))
              }
              if (motif_slots(motif, "type") == "PPM") {
                motif2 <- seqLogo::makePWM(motif_slots(motif, "motif"))
              }
              if (motif_slots(motif, "type") == "PWM") {
                motif2 <- convert_type(motif, "PPM")
                motif2 <- seqLogo::makePWM(motif_slots(motif2, "motif"))
              }
              if (motif_slots(motif, "type") == "ICM") {
                stop("conversion from ICM not currently supported")
              }
            }

            if (out_class == "motifStack-pcm") {
              mstackpcm <- getClass("pcm", where = "motifStack")
              motif2 <- convert_type(motif, "PCM")
              motif2 <- new(mstackpcm, mat = motif_slots(motif2, "motif"),
                            name = motif_slots(motif2, "name"),
                            alphabet = motif_slots(motif2, "alphabet"),
                            background = motif_slots(motif2, "bkg"))
            }

            if (out_class == "motifStack-pfm") {
              mstackpfm <- getClass("pfm", where = "motifStack")
              motif2 <- convert_type(motif, "PPM")
              motif2 <- new(mstackpfm, mat = motif_slots(motif2, "motif"),
                            name = motif_slots(motif, "name"),
                            alphabet = motif_slots(motif, "alphabet"),
                            background = motif_slots(motif, "bkg"))
            }

            if (out_class == "PWMEnrich-PWM") {
              motif2 <- convert_type(motif, "PCM")
              motpwm <- convert_type(motif, "PWM")
              penrpwm <- getClass("PWM", where = "PWMEnrich")
              biomat <- matrix(as.integer(motif_slots(motif2, "motif")),
                               byrow = FALSE, nrow = 4)
              rownames(biomat) <- c("A", "C", "G", "T")
              biopriors <- motif_slots(motif2, "bkg")
              names(biopriors) <- c("A", "C", "G", "T")
              biomat <- PWMEnrich::PFMtoPWM(biomat, type = "log2probratio",
                                            prior.params = biopriors,
                                            pseudo.count = motif_slots(motif2,
                                                                       "pseudoweight"))
              motif2 <- new(penrpwm, name = motif_slots(motif2, "name"),
                            pfm = motif_slots(motif2, "motif"),
                            prior.params = biopriors,
                            pwm = biomat$pwm)
            }

            if (out_class == "Biostrings-PWM") {
              motif2 <- convert_type(motif, "PCM")
              biomat <- matrix(as.integer(motif_slots(motif2, "motif")),
                               byrow = FALSE, nrow = 4)
              rownames(biomat) <- c("A", "C", "G", "T")
              biopriors <- motif_slots(motif2, "bkg")
              names(biopriors) <- c("A", "C", "G", "T")
              motif2 <- Biostrings::PWM(x = biomat, 
                                        type = "log2probratio",
                                        prior.params = biopriors)
            }

            if (is.null(motif2)) stop("unknown 'out_class'")
            return(motif2)
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
  if (motif@type == "ICM") {
    stop("ICM conversion currently not supported")
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
