#' Convert motif class.
#'
#' Allows for easy transfer of motif information between different classes
#' defined by other Bioconductor packages. This function is also used by
#' nearly all other functions in this package; so any motifs of a compatible
#' class can be used without needed to convert beforehand.
#'
#' @param motifs Single motif object or list.
#' @param class Desired motif class. Input as 'package-class'. If left empty,
#'              defaults to 'universalmotif-universalmotif'. (See details.)
#' @param BPPARAM See \code{\link[BiocParallel]{SerialParam}}.
#'
#' @return Single motif object or list.
#'
#' @details
#'    The following packge-class combinations can be used as input:
#'    \itemize{
#'       \item MotifDb-MotifList
#'       \item TFBSTools-PFMatrix
#'       \item TFBSTools-PWMatrix
#'       \item TFBSTools-ICMatrix
#'       \item TFBSTools-PFMatrixList
#'       \item TFBSTools-PWMatrixList
#'       \item TFBSTools-ICMatrixList
#'       \item TFBSTools-TFFMFirst
#'       \item seqLogo-pwm
#'       \item motifStack-pcm
#'       \item motifStack-pfm
#'       \item PWMEnrich-PWM
#'       \item motifRG-Motif
#'       \item universalmotif-universalmotif
#'    }
#'
#'    The following package-class combinations can be output:
#'    \itemize{
#'       \item MotIV-pwm2
#'       \item TFBSTools-PFMatrix
#'       \item TFBSTools-PWMatrix
#'       \item TFBSTools-ICMatrix
#'       \item TFBSTools-TFFMFirst
#'       \item seqLogo-pwm
#'       \item motifStack-pcm
#'       \item motifStack-pfm
#'       \item PWMEnrich-PWM
#'       \item Biostrings-PWM (\code{type = 'log2prob'})
#'       \item rGADEM-motif
#'       \item universalmotif-universalmotif
#'    }
#'
#' @examples
#' # convert from universalmotif:
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.motifstack.pfm <- convert_motifs(jaspar, "motifStack-pfm")
#'
#' # convert from another class to universalmotif:
#' library(TFBSTools)
#' data(MA0003.2)
#' motif <- convert_motifs(MA0003.2)
#'
#' # convert from another class to another class
#' motif <- convert_motifs(MA0003.2, "PWMEnrich-PWM")
#'
#' # the 'convert_motifs' function is embedded in the rest of the universalmotif
#' # functions; non-universalmotif class motifs can be used
#' MA0003.2.trimmed <- trim_motifs(MA0003.2)
#' # note: if the motif object going in has information that the
#' # 'universalmotif' class can't hold, it will be lost
#'
#' @references
#'    \insertRef{motiv}{universalmotif}
#'
#'    \insertRef{motifdb}{universalmotif}
#'
#'    \insertRef{tfbstools}{universalmotif}
#'
#'    \insertRef{seqlogo}{universalmotif}
#'
#'    \insertRef{motifstack}{universalmotif}
#'
#'    \insertRef{pwmenrich}{universalmotif}
#'
#'    \insertRef{motifrg}{universalmotif}
#'
#'    \insertRef{rgadem}{universalmotif}
#'
#'    \insertRef{biostrings}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @include universalmotif-class.R
#' @export
setGeneric("convert_motifs", function(motifs,
                                      class = "universalmotif-universalmotif",
                                      BPPARAM = SerialParam())
           standardGeneric("convert_motifs"))

#' @describeIn convert_motifs Convert a list of motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "list"),
          definition = function(motifs, class, BPPARAM) {
            mot_classes <- vapply(motifs, function(x) class(x), character(1))
            mot_classes <- unique(mot_classes)
            if (length(mot_classes) == 1) {
              classin <- strsplit(class, "-")[[1]][2]
              if (mot_classes == classin) return(motifs)
            }
            bplapply(motifs, function(x) convert_motifs(x, class = class,
                                                        BPPARAM = BPPARAM))
          })

#' @describeIn convert_motifs Convert a \linkS4class{universalmotif} object.
#' @export
setMethod("convert_motifs", signature(motifs = "universalmotif"),
          definition = function(motifs, class, BPPARAM) {
            
            out_class <- strsplit(class, "-")[[1]][2]
            out_class_pkg <- strsplit(class, "-")[[1]][1]

            if (out_class_pkg == "universalmotif") {
              return(motifs)
            }

            # MotIV-pwm2
            if (out_class_pkg == "MotIV" && out_class == "pwm2") {
              motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)
              motifs <- makePWM(motifs["motif"],
                                alphabet = motifs["alphabet"])
              return(motifs)
            }

            # TFBSTools- PFMatrix, PWMatrix, and ICMatrix
            if (out_class_pkg == "TFBSTools" && out_class %in% 
                c("PFMatrix", "PWMatrix", "ICMatrix")) {
              motifs <- convert_type(motifs, "PCM")
              bkg <- motifs["bkg"]
              names(bkg) <- DNA_BASES
              extrainfo <- motifs["extrainfo"]
              if (length(motifs["altname"]) == 0) {
                motifs["altname"] <- ""
              }
              strand <- motifs["strand"]
              if (strand %in% c("+-", "-+")) {
                strand <- "*"
              }
              motifs <- PFMatrix(name = motifs["name"],
                                 ID = motifs["altname"],
                                 strand = strand,
                                 bg = bkg,
                                 profileMatrix = motifs["motif"])
              if (out_class == "PFMatrix") {
                motifs <- motifs
              } else if (out_class == "PWMatrix") {
                motifs <- toPWM(motifs, type = "log2probratio",
                                pseudocounts = 0.8,
                                bg = bkg)
              } else if (out_class == "ICMatrix") {
                motifs <- toICM(motifs, pseudocounts = 0.8,
                                bg = bkg)
              }
              if (length(extrainfo) > 0) {
                motifs@tags <- as.list(extrainfo)
              }
              return(motifs)
            }

            # TFBSTools-TFFMFirst
            if (out_class_pkg == "TFBSTools" && out_class == "TFFMFirst") {
              motifs <- convert_type(motifs, "PPM")
              if (!`2` %in% names(motifs@multifreq)) {
                stop("cannot convert without filled multifreq slot")
              }
              if (motifs["alphabet"] != "DNA") stop("alphabet must be DNA")
              bkg <- motifs["bkg"]
              emission <- list(length = ncol(motifs@multifreq$`2`) + 1)
              bkg <- motifs["bkg"]
              names(bkg) <- DNA_BASES
              emission[[1]] <- bkg
              transition <- matrix(rep(0, (ncol(motifs@multifreq$`2`) + 1) *
                                           ncol(motifs@multifreq$`2`)),
                                   nrow = ncol(motifs@multifreq$`2`) + 1,
                                   ncol = ncol(motifs@multifreq$`2`))
              for (i in seq_len(ncol(motifs@multifreq$`2`))) {
                emission[[i + 1]] <- motifs@multifreq$`2`[, i]
                transition[i, i] <- 1
              }
              names(emission) <- rep("state", length(emission))
              transition[nrow(transition), 1] <- 1
              transition[2, 1:2] <- c(0.95, 0.05)
              colnames(transition) <- seq_len(ncol(transition))
              rownames(transition) <- c(0, seq_len(ncol(transition)))
              if (length(motifs["altname"]) < 1) motifs@altname <- "Unknown"
              strand <- motifs["strand"]
              if (strand %in% c("+-", "-+")) strand <- "*"
              family <- motifs["family"]
              if (length(family) < 1) family <- "Unknown"
              motifs <- TFFMFirst(ID = motifs["altname"], name = motifs["name"],
                                  strand = strand, type = "First",
                                  bg = bkg, matrixClass = family,
                                  profileMatrix = motifs["motif"],
                                  tags = as.list(motifs["extrainfo"]),
                                  emission = emission, transition = transition)
              return(motifs)
            }

            # seqLogo-pwm
            if (out_class_pkg == "seqLogo" && out_class == "pwm") {
              # if (!requireNamespace("seqLogo", quietly = TRUE)) {
                # stop("'seqLogo' package not installed")
              # }
              motifs <- convert_type(motifs, "PPM")
              motifs <- seqLogo::makePWM(motifs["motif"])
              return(motifs)
            }

            # motifStack-pcm
            if (out_class_pkg == "motifStack" && out_class == "pcm") {
              pcm_class <- getClass("pcm", where = "motifStack")
              motifs <- convert_type(motifs, "PCM")
              motifs <- new(pcm_class, mat = motifs["motif"],
                            name = motifs["name"],
                            alphabet = motifs["alphabet"],
                            background = motifs["bkg"])
              return(motifs)
            }

            # motifStack-pfm
            if (out_class_pkg == "motifStack" && out_class == "pfm") {
              pfm_class <- getClass("pfm", where = "motifStack")
              motifs <- convert_type(motifs, "PPM")
              motifs <- new(pfm_class, mat = motifs["motif"],
                            name = motifs["name"],
                            alphabet = motifs["alphabet"],
                            background = motifs["bkg"])
              return(motifs)
            }

            # PWMEnrich-PWM
            if (out_class_pkg == "PWMEnrich" && out_class == "PWM") {
              motifs <- convert_type(motifs, "PCM")
              PWM_class <- getClass("PWM", where = "PWMEnrich")
              bio_mat <- matrix(as.integer(motifs["motif"]), byrow = FALSE,
                                nrow = 4)
              rownames(bio_mat) <- DNA_BASES
              bio_priors <- motifs["bkg"]
              names(bio_priors) <- DNA_BASES
              bio_mat <- PFMtoPWM(bio_mat, type = "log2probratio",
                                  prior.params = bio_priors,
                                  pseudo.count = motifs["pseudocount"])
              motifs <- new(PWM_class, name = motifs["name"],
                            pfm = motifs["motif"],
                            prior.params = bio_priors,
                            pwm = bio_mat$pwm)
              return(motifs)
            }

            # Biostrings-PWM
            if (out_class_pkg == "Biostrings" && out_class == "PWM") {
              motifs <- convert_type(motifs, "PCM")
              bio_mat <- matrix(as.integer(motifs["motif"]), byrow = FALSE,
                                nrow = 4)
              rownames(bio_mat) <- DNA_BASES
              bio_priors <- motifs["bkg"]
              names(bio_priors) <- DNA_BASES
              motifs <- PWM(x = bio_mat, type = "log2probratio",
                            prior.params = bio_priors)
              return(motifs)
            }

            # rGADEM-motif
            if (out_class_pkg == "rGADEM" && out_class == "motif") {
              if (!requireNamespace("rGADEM", quietly = FALSE)) {
                stop("'rGADEM' package not installed")
              }
              motifs <- convert_type(motifs, "PPM")
              rGADEM_motif_class <- getClass("motif", where = "rGADEM")
              motifs <- new("motif", pwm = motifs["motif"],
                            name = motifs["name"],
                            consensus = motifs["consensus"])
              return(motifs)
            }

            # custom convert function
            # custom_convert <- NULL
            # tryCatch({custom_convert <- get(class)}, error = function(e) {})
            # if (is.function(custom_convert)) {
              # message("attempting to convert using custom function")
              # motifs <- custom_convert(motifs)
              # return(motifs)
            # }

            stop("unknown 'class'")
          
          })

#' @describeIn convert_motifs Convert \linkS4class{MotifList} motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "MotifList"),
          definition = function(motifs, class, BPPARAM) {
            motifdb_fun <- function(i) {
              x <- motifs
              mot <- universalmotif_cpp(
                             name = x@elementMetadata@listData$providerName[i],
                             altname = x@elementMetadata@listData$geneSymbol[i],
                             family = x@elementMetadata@listData$tfFamily[i],
                             organism = x@elementMetadata@listData$organism[i],
                             motif = x@listData[[i]], alphabet = "DNA",
                             type = "PPM")
              msg <- validObject_universalmotif(mot)
              if (length(msg) > 0) stop(msg)
              mot
            }
            motifs_out <- vector("list", length(motifs))
            motifs_out <- bplapply(seq_len(length(motifs)), motifdb_fun,
                                   BPPARAM = BPPARAM)
            convert_motifs(motifs_out, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert TFFMFirst motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "TFFMFirst"),
          definition = function(motifs, class, BPPARAM) {
            difreq <- motifs@emission[-1]
            difreq <- do.call(c, difreq)
            difreq <- matrix(difreq, nrow = 16) / 4
            rownames(difreq) <- DNA_DI
            colnames(difreq) <- seq_len(ncol(difreq))
            mot <- universalmotif(name = motifs@name, altname = motifs@ID,
                           strand = motifs@strand, bkg = motifs@bg,
                           motif = getPosProb(motifs),
                           multifreq = list(`2` = difreq))
            convert_motifs(motifs_out, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert \linkS4class{PFMatrix} motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "PFMatrix"),
          definition = function(motifs, class, BPPARAM) {
            if (all(names(motifs@bg) %in% DNA_BASES)) {
              alphabet <- "DNA"
            } else alphabet  <- "RNA"
            extrainfo <- motifs@tags
            if (length(extrainfo) > 0) {
              extrainfo <- unlist(extrainfo)
            } else {
              extrainfo <- character()
            }
            motifs <- universalmotif_cpp(name = motifs@name, altname = motifs@ID,
                                     family = motifs@tags$family,
                                     organism = motifs@tags$species,
                                     motif = motifs@profileMatrix,
                                     alphabet = alphabet, type = "PCM",
                                     bkg = motifs@bg,
                                     strand = paste0(motifs@strand,
                                                     collapse = ""),
                                     extrainfo = extrainfo)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert \linkS4class{PWMatrix} motifs.
setMethod("convert_motifs", signature(motifs = "PWMatrix"),
          definition = function(motifs, class, BPPARAM) {
            if (all(names(motifs@bg) %in% DNA_BASES)) {
              alphabet <- "DNA"
            } else alphabet  <- "RNA"
            extrainfo <- motifs@tags
            if (length(extrainfo) > 0) {
              extrainfo <- unlist(extrainfo)
            } else {
              extrainfo <- character()
            }
            motifs <- universalmotif_cpp(name = motifs@name, altname = motifs@ID,
                                     family = motifs@tags$family,
                                     organism = motifs@tags$species,
                                     motif = motifs@profileMatrix,
                                     alphabet = alphabet, type = "PWM",
                                     bkg = motifs@bg,
                                     strand = paste0(motifs@strand,
                                                     collapse = ""),
                                     extrainfo = extrainfo)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert \linkS4class{ICMatrix} motifs.
setMethod("convert_motifs", signature(motifs = "ICMatrix"),
          definition = function(motifs, class, BPPARAM) {
            if (all(names(motifs@bg) %in% DNA_BASES)) {
              alphabet <- "DNA"
            } else alphabet  <- "RNA"
            extrainfo <- motifs@tags
            if (length(extrainfo) > 0) {
              extrainfo <- unlist(extrainfo)
            } else {
              extrainfo <- character()
            }
            motifs <- universalmotif_cpp(name = motifs@name, altname = motifs@ID,
                                     family = motifs@tags$family,
                                     organism = motifs@tags$species,
                                     motif = motifs@profileMatrix,
                                     alphabet = alphabet, type = "ICM",
                                     bkg = motifs@bg,
                                     strand = paste0(motifs@strand,
                                                     collapse = ""),
                                     extrainfo = extrainfo)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert XMatrixList motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "XMatrixList"),
          definition = function(motifs, class, BPPARAM) {
            motif_num <- length(motifs@listData)
            motifs_out <- bplapply(seq_len(motif_num),
                                   function(i) {
                                     motifs@listData[[i]]
                                   }, BPPARAM = BPPARAM)
            motif_names <- unlist(bplapply(seq_len(motif_num),
                                           function(i) {
                                            motifs@listData[[i]]@name
                                           }, BPPARAM = BPPARAM))
            names(motifs_out) <- motif_names
            motifs <- convert_motifs(motifs_out, class = "universalmotif",
                                     BPPARAM = BPPARAM)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert \linkS4class{pwm} motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "pwm"),
          definition = function(motifs, class, BPPARAM) {
            motifs <- universalmotif_cpp(motif = motifs@pwm, type = "PPM",
                                     alphabet = motifs@alphabet)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert \linkS4class{pcm} motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "pcm"),
          definition = function(motifs, class, BPPARAM) {
            motifs <- universalmotif_cpp(name = motifs@name, motif = motifs@mat,
                                     nsites = unique(colSums(motifs@mat))[1],
                                     alphabet = motifs@alphabet,
                                     bkg = motifs@background,
                                     type = "PCM")
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert \linkS4class{pfm} motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "pfm"),
          definition = function(motifs, class, BPPARAM) {
            motifs <- universalmotif_cpp(name = motifs@name, motif = motifs@mat,
                                     alphabet = motifs@alphabet,
                                     bkg = motifs@background,
                                     type = "PPM")
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert \linkS4class{PWM} motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "PWM"),
          definition = function(motifs, class, BPPARAM) {
            if (all(names(motifs@pwm) %in% DNA_BASES)) {
              alphabet <- "DNA"
            } else alphabet <- "RNA"
            motifs <- universalmotif_cpp(name = motifs@name, motif = motifs@pwm,
                                     type = "PWM", alphabet = alphabet,
                                     bkg = motifs@prior.params,
                                     altname = motifs@id)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

#' @describeIn convert_motifs Convert Motif motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "Motif"),
          definition = function(motifs, class, BPPARAM) {
            motifs <- universalmotif_cpp(name = motifs@pattern,
                                     nsites = sum(motifs@count),
                                     alphabet = "DNA",
                                     type = "PCM",
                                     extrainfo = c(score = motifs@score),
                             strand = paste(unique(motifs@match$match.strand),
                                            collapse = ""),
            motif = create_motif(input = DNAStringSet(motifs@match$pattern))@motif)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class, BPPARAM = BPPARAM)
          })

# @describeIn convert_motifs Convert non-\linkS4class{universalmotif} class motifs.
# @export
# setMethod("convert_motifs", signature(motifs = "ANY"),
          # definition = function(motifs, class, BPPARAM) {

            # success <- FALSE

            # ## convert to universalmotif
            # in_class <- class(motifs)[1]
            # in_class_pkg <- attributes(class(motifs))$package

            # if (paste(in_class_pkg, in_class, sep = "-") == class) {
              # return(motifs)
            # }

            # paste(in_class_pkg)
            # paste(in_class)
            # if (!success) stop("unrecognized class")

            # ## convert to desired class
            # motifs <- convert_motifs(motifs, class = class)

            # motifs

          # })
