######################################################################
## Benjamin Tremblay
##
## Convert motifs between different classes. convert_motifs works on
## one object at a time, so lapply must be used on lists.
##
######################################################################

#' @describeIn convert_motifs Convert a list of motifs.
setMethod("convert_motifs", signature = "list",
          definition = function(motif, out_class = "universalmotif",
                                pseudoweight = 0.8,
                                background = c("A" = 0.25, "C" = 0.25,
                                               "G" = 0.25, "T" = 0.25),
                                ...) {
            lapply(X = motif, FUN = convert_motifs, out_class = out_class,
                   pseudoweight = pseudoweight, background = background, ...)
          })

#' @describeIn convert_motifs Convert from \linkS4class{PFMatrix} (TFBSTools).
setMethod("convert_motifs", signature = "PFMatrix",
          definition = function(motif, out_class = "universalmotif", ...) {
            if (all(names(motif@bg) == c("A", "C", "G", "T"))) {
              alphabet <- "DNA"
            } else alphabet <- "RNA"
            motif <- universalmotif(name = motif@name, strand = motif@strand,
                                    bkg = motif@bg, motif = motif@profileMatrix,
                                    alphabet = alphabet,
                                    extrachar = c("ID" = motif@ID,
                                                  "matrixClass" = motif@matrixClass,
                                                  unlist(motif@tags)))
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{PFMatrixList} (TFBSTools).
setMethod("convert_motifs", signature = "PFMatrixList",
          definition = function(motif, out_class = "universalmotif", ...) {
            nummots <- length(motif@listData)
            motiflist <- list()
            motnames <- vector(length = nummots)
            for (i in seq_len(nummots)) {
              motiflist[[i]] <- motif@listData[[i]]
              motnames[i] <- motif@listData[[i]]@name
            }
            motifs <- convert_motifs(motiflist, out_class = out_class)
            names(motifs) <- motnames
            return(motifs)
          })

#' @describeIn convert_motifs Convert from \linkS4class{PWMatrix} (TFBSTools).
setMethod("convert_motifs", signature = "PWMatrix",
          definition = function(motif, out_class = "universalmotif", ...) {
            if (all(names(motif@bg) == c("A", "C", "G", "T"))) {
              alphabet <- "DNA"
            } else alphabet <- "RNA"
            motif <- universalmotif(name = motif@name, strand = motif@strand,
                                    bkg = motif@bg, motif = motif@profileMatrix,
                                    type = "PWM", alphabet = alphabet,
                                    extrachar = c("ID" = motif@ID,
                                                  "matrixClass" = motif@matrixClass,
                                                  unlist(motif@tags)))
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{PFMatrixList} (TFBSTools).
setMethod("convert_motifs", signature = "PWMatrixList",
          definition = function(motif, out_class = "universalmotif", ...) {
            nummots <- length(motif@listData)
            motiflist <- list()
            motnames <- vector(length = nummots)
            for (i in seq_len(nummots)) {
              motiflist[[i]] <- motif@listData[[i]]
              motnames[i] <- motif@listData[[i]]@name
            }
            motifs <- convert_motifs(motiflist, out_class = out_class)
            names(motifs) <- motnames
            return(motifs)
          })

#' @describeIn convert_motifs Convert fom \linkS4class{ICMatrix} (TFBSTools).
setMethod("convert_motifs", signature = "ICMatrix",
          definition = function(motif, out_class = "universalmotif", ...) {
            if (all(names(motif@bg) == c("A", "C", "G", "T"))) {
              alphabet  <- "DNA"
            } else alphabet  <- "RNA"
            motif <- universalmotif(name = motif@name, strand = motif@strand,
                                    bkg = motif@bg, motif = motif@profileMatrix,
                                    type = "ICM", alphabet = alphabet,
                                    extrachar = c("ID" = motif@ID,
                                                  "matrixClass" = motif@matrixClass,
                                                  unlist(motif@tags)))
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{PFMatrixList} (TFBSTools).
setMethod("convert_motifs", signature = "ICMatrixList",
          definition = function(motif, out_class = "universalmotif", ...) {
            nummots <- length(motif@listData)
            motiflist <- list()
            motnames <- vector(length = nummots)
            for (i in seq_len(nummots)) {
              motiflist[[i]] <- motif@listData[[i]]
              motnames[i] <- motif@listData[[i]]@name
            }
            motifs <- convert_motifs(motiflist, out_class = out_class)
            names(motifs) <- motnames
            return(motifs)
          })

#' @describeIn convert_motifs Convert from \linkS4class{pwm} (seqLogo).
setMethod("convert_motifs", signature = "pwm",
          definition = function(motif, out_class = "universalmotif",
                                name, ...) {
            motif <- universalmotif(name = name, motif = motif@pwm,
                                    type = "PPM", alphabet = motif@alphabet)
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{pcm} (motifStack).
setMethod("convert_motifs", signature = "pcm",
          definition = function(motif, out_class = "universalmotif", ...) {
            motif <- universalmotif(name = motif@name, motif = motif@mat,
                                    type = "PCM", alphabet = motif@alphabet,
                                    bkg = motif@background)
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{pfm} (motifStack).
setMethod("convert_motifs", signature = "pfm",
          definition = function(motif, out_class = "universalmotif", ...) {
            motif <- universalmotif(name = motif@name, motif = motif@mat,
                                    type = "PPM", alphabet = motif@alphabet,
                                    bkg = motif@background)
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{PWM} (PWMEnrich).
setMethod("convert_motifs", signature = "PWM",
          definition = function(motif, out_class = "universalmotif", ...) {
            if (all(names(motif@pwm) == c("A", "C", "G", "T"))) {
              alphabet <- "DNA"
            } else alphabet <- "RNA"
            motif <- universalmotif(name = motif@name, motif = motif@pwm,
                                    type = "PWM", alphabet = alphabet,
                                    bkg = motif@prior.params,
                                    extrachar = c("ID" = motif@id))
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{MotifList} (MotifDb).
setMethod("convert_motifs", signature = "MotifList",
          definition = function(motif, out_class = "universalmotif", ...) {
            
            # system.time(convert_motifs(MotifDb)) --> 8 seconds
            # system.time(convert_motifs(MotifDb, out_class = "PWMEnrich-PWM"))
            #     --> 21 seconds

            motnum <- length(motif)
            motnames <- names(motif@listData)
            motifs <- list()
            for (i in seq_len(motnum)) {
              motifs[[motnames[i]]] <- universalmotif(name = motnames[i],
                                                   motif = motif@listData[[i]],
                                                   alphabet = "DNA",
                                                   type = "PPM",
                                                   extrachar = c(

         "providerName" = motif@elementMetadata@listData$providerName[i],
         "providerId" = motif@elementMetadata@listData$providerId[i],
         "dataSource" = motif@elementMetadata@listData$dataSource[i],
         "geneSymbol" = motif@elementMetadata@listData$geneSymbol[i],
         "geneId" = motif@elementMetadata@listData$geneId[i],
         "geneIdType" = motif@elementMetadata@listData$geneIdType[i],
         "proteinId" = motif@elementMetadata@listData$proteinId[i],
         "proteinIdType" = motif@elementMetadata@listData$proteinIdType[i],
         "organism" = motif@elementMetadata@listData$organism[i],
         "bindingSequence" = motif@elementMetadata@listData$bindingSequence[i],
         "tfFamily" = motif@elementMetadata@listData$tfFamily[i],
         "experimentType" = motif@elementMetadata@listData$experimentType[i],
         "pubmedID" = motif@elementMetadata@listData$pubmedID[i],
         "sequenceCount" = motif@elementMetadata@listData$sequenceCount[i]
         
         ))

            }
            if (out_class == "universalmotif") return(motifs)
            motifs <- convert_motifs(motifs, out_class = out_class)
            return(motifs)
          })


## pwm2 is NOT exported by MotIV! Would exporting the class myself allow it
## to interact with MotIV-derived 'pwm2' class objects out in the wild?

#' @describeIn convert_motifs Convert from \linkS4class{pwm2} (MotIV).
setMethod("convert_motifs", signature = "pwm2",
          definition = function(motif, out_class = "universalmotif",
                                name, ...) {
            motif <- universalmotif(name = name, motif = motif@pwm,
                                    type = "PPM", alphabet = motif@alphabet)
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{motif} (rGADEM).
setMethod("convert_motifs", signature = "motif",
          definition = function(motif, out_class = "universalmotif", ...) {
            motif <- universalmotif(name = motif@name, motif = motif@pwm)
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class)
            return(motif)
          })

#' @describeIn convert_motifs Convert from \linkS4class{universalmotif}.
#' @include utils.R universalmotif-class.R generics.R
setMethod("convert_motifs", signature = "universalmotif",
          definition = function(motif, out_class, pseudoweight = 0.8,
                                background = c("A" = 0.25, "C" = 0.25,
                                               "G" = 0.25, "T" = 0.25)) {

            if (out_class == "universalmotif") return(motif)
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

            if (out_class == "MotIV-pwm2") {
              motif2 <- convert_type(motif, "PPM")
              motif2 <- MotIV::makePWM(motif2@motif,
                                       alphabet = motif2@alphabet)
            }

            if (out_class == "rGADEM-motif") {
              motif <- convert_type(motif, "PPM")
              rgadmot <- getClass("motif", where = "rGADEM")
              motif2 <- new("motif", pwm = motif@motif, name = motif@name,
                            consensus = motif@consensus)
            }

            if (is.null(motif2)) stop("unknown 'out_class'")
            return(motif2)
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
