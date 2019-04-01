#' Convert motif class.
#'
#' Allows for easy transfer of motif information between different classes
#' defined by other Bioconductor packages. This function is also used by
#' nearly all other functions in this package; so any motifs of a compatible
#' class can be used without needed to convert beforehand.
#'
#' @param motifs Single motif object or list. See details.
#' @param class `character(1)` Desired motif class. Input as
#'    'package-class'. If left empty, defaults to
#'    'universalmotif-universalmotif'. (See details.)
#'
#' @return Single motif object or list.
#'
#' @details
#' The following packge-class combinations can be used as input:
#' * MotifDb-MotifList
#' * TFBSTools-PFMatrix
#' * TFBSTools-PWMatrix
#' * TFBSTools-ICMatrix
#' * TFBSTools-PFMatrixList
#' * TFBSTools-PWMatrixList
#' * TFBSTools-ICMatrixList
#' * TFBSTools-TFFMFirst
#' * seqLogo-pwm
#' * motifStack-pcm
#' * motifStack-pfm
#' * PWMEnrich-PWM
#' * motifRG-Motif
#' * universalmotif-universalmotif
#' * matrix
#'
#' The following package-class combinations can be output:
#' * MotIV-pwm2
#' * TFBSTools-PFMatrix
#' * TFBSTools-PWMatrix
#' * TFBSTools-ICMatrix
#' * TFBSTools-TFFMFirst
#' * seqLogo-pwm
#' * motifStack-pcm
#' * motifStack-pfm
#' * PWMEnrich-PWM
#' * Biostrings-PWM (\code{type = 'log2prob'})
#' * rGADEM-motif
#' * universalmotif-universalmotif
#'
#' @examples
#' # convert from universalmotif:
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' if (requireNamespace("motifStack", quietly = TRUE)) {
#'   jaspar.motifstack.pfm <- convert_motifs(jaspar, "motifStack-pfm")
#' }
#'
#' # convert from another class to universalmotif:
#' if (requireNamespace("TFBSTools", quietly = TRUE)) {
#' library(TFBSTools)
#' data(MA0003.2)
#' motif <- convert_motifs(MA0003.2)
#'
#' # convert from another class to another class
#' if (requireNamespace("PWMEnrich", quietly = TRUE)) {
#'   motif <- convert_motifs(MA0003.2, "PWMEnrich-PWM")
#' }
#'
#' # the 'convert_motifs' function is embedded in the rest of the universalmotif
#' # functions; non-universalmotif class motifs can be used
#' MA0003.2.trimmed <- trim_motifs(MA0003.2)
#' # note: if the motif object going in has information that the
#' # 'universalmotif' class can't hold, it will be lost
#' }
#'
#' @references
#'    \insertRef{seqlogo}{universalmotif}
#'
#'    \insertRef{rgadem}{universalmotif}
#'
#'    \insertRef{motiv}{universalmotif}
#'
#'    \insertRef{motifstack}{universalmotif}
#'
#'    \insertRef{biostrings}{universalmotif}
#'
#'    \insertRef{motifdb}{universalmotif}
#'
#'    \insertRef{pwmenrich}{universalmotif}
#'
#'    \insertRef{tfbstools}{universalmotif}
#'
#'    \insertRef{motifrg}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @include universalmotif-class.R
#' @export
setGeneric("convert_motifs", function(motifs,
                                      class = "universalmotif-universalmotif")
           standardGeneric("convert_motifs"))

#' @describeIn convert_motifs Convert a list of motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "list"),
          definition = function(motifs, class) {
            mot_classes <- vapply(motifs, function(x) class(x), character(1))
            mot_classes <- unique(mot_classes)
            if (length(mot_classes) == 1) {
              classin <- strsplit(class, "-")[[1]][2]
              if (mot_classes == classin) return(motifs)
            }
            lapply(motifs, function(x) convert_motifs(x, class = class))
          })

#' @describeIn convert_motifs Convert a \linkS4class{universalmotif} object.
#' @export
setMethod("convert_motifs", signature(motifs = "universalmotif"),
          definition = function(motifs, class) {

            out_class <- strsplit(class, "-")[[1]][2]
            out_class_pkg <- strsplit(class, "-")[[1]][1]

            switch(out_class_pkg,
              "universalmotif" = {
                motifs
              },
              "MotIV" = {
                if (out_class == "pwm2")
                  convert_to_motiv_pwm2(motifs)
                else
                  stop("unknown 'class'")
              },
              "TFBSTools" = {
                if (out_class %in% c("PFMatrix", "PWMatrix", "ICMatrix"))
                  convert_to_tfbstools_matrix(motifs, out_class)
                else if (out_class == "TFFMFirst")
                  convert_to_tfbstools_tffmfirst(motifs)
                else
                  stop("unknown 'class'")
              },
              "seqLogo" = {
                if (out_class == "pwm")
                  convert_to_seqlogo_pwm(motifs)
                else
                  stop("unknown 'class'")
              },
              "motifStack" = {
                switch(out_class,
                       "pcm" = convert_to_motifstack_pcm(motifs),
                       "pfm" = convert_to_motifstack_pfm(motifs),
                               stop("unknown 'class'"))
              },
              "PWMEnrich" = {
                if (out_class == "PWM")
                  convert_to_pwmenrich_pwm(motifs)
                else
                  stop("unknown 'class'")
              },
              "Biostrings" = {
                if (out_class == "PWM")
                  convert_to_biostrings_pwm(motifs)
                else
                  stop("unknown 'class'")
              },
              "rGADEM" = {
                if (out_class == "motif")
                  convert_to_rgadem_motif(motifs)
                else
                  stop("unknown 'class'")
              },
              stop("unknown 'class'")
            )

          })

convert_to_motiv_pwm2 <- function(motifs) {
  motifs <- convert_type(motifs, "PPM")
  if (requireNamespace("MotIV", quietly = TRUE)) {
    motifs <- MotIV::makePWM(motifs["motif"],
                             alphabet = motifs["alphabet"])
  } else {
    stop("package 'MotIV' is not installed")
  }
  motifs
}

convert_to_tfbstools_matrix <- function(motifs, out_class) {
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
  if (requireNamespace("TFBSTools", quietly = TRUE)) {
    motifs <- TFBSTools::PFMatrix(name = motifs["name"],
                       ID = motifs["altname"],
                       strand = strand, bg = bkg,
                       profileMatrix = motifs["motif"])
    switch(out_class,
      "PFMatrix" = {
        motifs <- motifs
      },
      "PWMatrix" = {
        motifs <- TFBSTools::toPWM(motifs, type = "log2probratio",
                        pseudocounts = 0.8, bg = bkg)
      },
      "ICMatrix" = {
        motifs <- TFBSTools::toICM(motifs, pseudocounts = 0.8, bg = bkg)
      }
    )
  } else {
    stop("package 'TFBSTools' is not installed") 
  }
  if (length(extrainfo) > 0) {
    motifs@tags <- as.list(extrainfo)
  }
  motifs
}

convert_to_tfbstools_tffmfirst <- function(motifs) {
  motifs <- convert_type(motifs, "PPM")
  if (!`2` %in% names(motifs@multifreq)) {
    stop("cannot convert without filled multifreq slot")
  }
  if (motifs["alphabet"] != "DNA") stop("alphabet must be DNA")
  bkg <- motifs["bkg"]
  emission <- list(length = ncol(motifs@multifreq[["2"]]) + 1)
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
  if (requireNamespace("TFBSTools", quietly = TRUE)) {
    motifs <- TFBSTools::TFFMFirst(ID = motifs["altname"],
                        name = motifs["name"],
                        strand = strand, type = "First",
                        bg = bkg, matrixClass = family,
                        profileMatrix = motifs["motif"],
                        tags = as.list(motifs["extrainfo"]),
                        emission = emission, transition = transition)
  } else {
    stop("package 'TFBSTools' is not installed") 
  }
  motifs
}

convert_to_seqlogo_pwm <- function(motifs) {
  motifs <- convert_type(motifs, "PPM")
  if (requireNamespace("seqLogo", quietly = TRUE)) {
    motifs <- seqLogo::makePWM(motifs["motif"])
  } else {
    stop("'seqLogo' package not installed")
  }
  motifs
}

convert_to_motifstack_pcm <- function(motifs) {
  if (requireNamespace("motifStack", quietly = TRUE)) {
    pcm_class <- getClass("pcm", where = "motifStack")
    motifs <- convert_type(motifs, "PCM")
    motifs <- new(pcm_class, mat = motifs["motif"],
                  name = motifs["name"],
                  alphabet = motifs["alphabet"],
                  background = motifs["bkg"])
  } else {
    stop("'motifStack' package not installed")
  }
  colnames(motifs@mat) <- seq_len(ncol(motifs@mat))
  names(motifs@background) <- rownames(motifs@mat)
  motifs
}

convert_to_motifstack_pfm <- function(motifs) {
  if (requireNamespace("motifStack", quietly = TRUE)) {
    pfm_class <- getClass("pfm", where = "motifStack")
    motifs <- convert_type(motifs, "PPM")
    motifs <- new(pfm_class, mat = motifs["motif"],
                  name = motifs["name"],
                  alphabet = motifs["alphabet"],
                  background = motifs["bkg"])
  } else {
    stop("'motifStack' package not installed")
  }
  colnames(motifs@mat) <- seq_len(ncol(motifs@mat))
  names(motifs@background) <- rownames(motifs@mat)
  motifs
}

convert_to_pwmenrich_pwm <- function(motifs) {
  if (requireNamespace("PWMEnrich", quietly = TRUE)) {
    motifs <- convert_type(motifs, "PCM")
    PWM_class <- getClass("PWM", where = "PWMEnrich")
    bio_mat <- matrix(as.integer(motifs["motif"]), byrow = FALSE,
                      nrow = 4)
    rownames(bio_mat) <- DNA_BASES
    bio_priors <- motifs["bkg"]
    names(bio_priors) <- DNA_BASES
    bio_mat <- PWMEnrich::PFMtoPWM(bio_mat, type = "log2probratio",
                        prior.params = bio_priors,
                        pseudo.count = motifs["pseudocount"])
    motifs <- new(PWM_class, name = motifs["name"],
                  pfm = motifs["motif"],
                  prior.params = bio_priors,
                  pwm = bio_mat$pwm)
  } else {
    stop("package 'PWMEnrich' not installed")
  }
  motifs
}

convert_to_biostrings_pwm <- function(motifs) {
  motifs <- convert_type(motifs, "PCM")
  bio_mat <- matrix(as.integer(motifs["motif"]), byrow = FALSE,
                    nrow = 4)
  rownames(bio_mat) <- DNA_BASES
  bio_priors <- motifs["bkg"]
  names(bio_priors) <- DNA_BASES
  motifs <- PWM(x = bio_mat, type = "log2probratio",
                prior.params = bio_priors)
  motifs
}

convert_to_rgadem_motif <- function(motifs) {
  if (requireNamespace("rGADEM", quietly = FALSE)) {
    motifs <- convert_type(motifs, "PPM")
    rGADEM_motif_class <- getClass("motif", where = "rGADEM")
    motifs <- new("motif", pwm = motifs["motif"],
                  name = motifs["name"],
                  consensus = motifs["consensus"])
  } else {
    stop("'rGADEM' package not installed")
  }
  motifs
}

#' @describeIn convert_motifs Convert MotifList motifs. (\pkg{MotifDb})
#' @export
setMethod("convert_motifs", signature(motifs = "MotifList"),
          definition = function(motifs, class) {
            motifdb_fun <- function(i) {
              x <- motifs
              mot <- universalmotif_cpp(
                             altname = x@elementMetadata@listData$providerName[i],
                             name = x@elementMetadata@listData$geneSymbol[i],
                             family = x@elementMetadata@listData$tfFamily[i],
                             organism = x@elementMetadata@listData$organism[i],
                             motif = x@listData[[i]], alphabet = "DNA",
                             type = "PPM")
              msg <- validObject_universalmotif(mot)
              if (length(msg) > 0) stop(msg)
              mot
            }
            motifs_out <- vector("list", length(motifs))
            motifs_out <- lapply(seq_len(length(motifs)), motifdb_fun)
            convert_motifs(motifs_out, class = class)
          })

#' @describeIn convert_motifs Convert TFFMFirst motifs. (\pkg{TFBSTools})
#' @export
setMethod("convert_motifs", signature(motifs = "TFFMFirst"),
          definition = function(motifs, class) {
            if (requireNamespace("TFBSTools", quietly = TRUE)) {
              difreq <- motifs@emission[-1]
              difreq <- do.call(c, difreq)
              difreq <- matrix(difreq, nrow = 16) / 4
              rownames(difreq) <- DNA_DI
              colnames(difreq) <- seq_len(ncol(difreq))
              mot <- universalmotif(name = motifs@name, altname = motifs@ID,
                             strand = motifs@strand, bkg = motifs@bg,
                             motif = TFBSTools::getPosProb(motifs),
                             multifreq = list(`2` = difreq))
            } else {
              stop("package 'TFBSTools' is not installed") 
            }
            convert_motifs(mot, class = class)
          })

#' @describeIn convert_motifs Convert PFMatrix motifs. (\pkg{TFBSTools})
#' @export
setMethod("convert_motifs", signature(motifs = "PFMatrix"),
          definition = function(motifs, class) {
            if (all(names(motifs@bg) %in% DNA_BASES)) {
              alphabet <- "DNA"
            } else alphabet  <- "RNA"
            extrainfo <- motifs@tags
            if (length(extrainfo) > 0) {
              extrainfo <- unlist(extrainfo)
            } else {
              extrainfo <- character()
            }
            nsites <- sum(motifs@profileMatrix[, 1])
            motifs <- universalmotif_cpp(name = motifs@name, altname = motifs@ID,
                                     family = motifs@tags$family,
                                     nsites = nsites,
                                     organism = paste(motifs@tags$species, collapse = "/"),
                                     motif = motifs@profileMatrix,
                                     alphabet = alphabet, type = "PCM",
                                     bkg = motifs@bg,
                                     strand = paste0(motifs@strand,
                                                     collapse = ""),
                                     extrainfo = extrainfo)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Convert PWMatrix motifs. (\pkg{TFBSTools})
setMethod("convert_motifs", signature(motifs = "PWMatrix"),
          definition = function(motifs, class) {
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
                                     organism = paste(motifs@tags$species, collapse = "/"),
                                     motif = motifs@profileMatrix,
                                     alphabet = alphabet, type = "PWM",
                                     bkg = motifs@bg,
                                     strand = paste0(motifs@strand,
                                                     collapse = ""),
                                     extrainfo = extrainfo)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Convert ICMatrix motifs. (\pkg{TFBSTools})
setMethod("convert_motifs", signature(motifs = "ICMatrix"),
          definition = function(motifs, class) {
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
                                     organism = paste(motifs@tags$species, collapse = "/"),
                                     motif = motifs@profileMatrix,
                                     alphabet = alphabet, type = "ICM",
                                     bkg = motifs@bg,
                                     strand = paste0(motifs@strand,
                                                     collapse = ""),
                                     extrainfo = extrainfo)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Convert XMatrixList motifs. (\pkg{TFBSTools})
#' @export
setMethod("convert_motifs", signature(motifs = "XMatrixList"),
          definition = function(motifs, class) {
            motif_num <- length(motifs@listData)
            motifs_out <- lapply(seq_len(motif_num),
                                   function(i) {
                                     motifs@listData[[i]]
                                   })
            motif_names <- unlist(lapply(seq_len(motif_num),
                                           function(i) {
                                            motifs@listData[[i]]@name
                                           }))
            names(motifs_out) <- motif_names
            motifs <- convert_motifs(motifs_out, class = "universalmotif")
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Convert pwm motifs. (\pkg{seqLogo})
#' @export
setMethod("convert_motifs", signature(motifs = "pwm"),
          definition = function(motifs, class) {
            motifs <- universalmotif_cpp(motif = motifs@pwm, type = "PPM",
                                     alphabet = motifs@alphabet)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Convert pcm motifs. (\pkg{motifStack})
#' @export
setMethod("convert_motifs", signature(motifs = "pcm"),
          definition = function(motifs, class) {
            motifs <- universalmotif_cpp(name = motifs@name, motif = motifs@mat,
                                     nsites = unique(colSums(motifs@mat))[1],
                                     alphabet = motifs@alphabet,
                                     bkg = motifs@background,
                                     type = "PCM")
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Convert pfm motifs. (\pkg{motifStack})
#' @export
setMethod("convert_motifs", signature(motifs = "pfm"),
          definition = function(motifs, class) {
            motifs <- universalmotif_cpp(name = motifs@name, motif = motifs@mat,
                                     alphabet = motifs@alphabet,
                                     bkg = motifs@background,
                                     type = "PPM")
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Convert PWM motifs. (\pkg{PWMEnrich})
#' @export
setMethod("convert_motifs", signature(motifs = "PWM"),
          definition = function(motifs, class) {
            if (all(names(motifs@pwm) %in% DNA_BASES)) {
              alphabet <- "DNA"
            } else alphabet <- "RNA"
            motifs <- universalmotif_cpp(name = motifs@name, motif = motifs@pwm,
                                     type = "PWM", alphabet = alphabet,
                                     bkg = motifs@prior.params,
                                     altname = motifs@id)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Convert Motif motifs. (\pkg{motifRG})
#' @export
setMethod("convert_motifs", signature(motifs = "Motif"),
          definition = function(motifs, class) {
            motifs <- universalmotif_cpp(name = motifs@pattern,
                                     nsites = sum(motifs@count),
                                     alphabet = "DNA",
                                     type = "PCM",
                                     extrainfo = c(score = motifs@score),
                             strand = paste(unique(motifs@match$match.strand),
                                            collapse = ""),
            motif <- create_motif(input = DNAStringSet(motifs@match$pattern))@motif)
            msg <- validObject_universalmotif(motifs)
            if (length(msg) > 0) stop(msg)
            convert_motifs(motifs, class = class)
          })

#' @describeIn convert_motifs Create motif from matrices.
#' @export
setMethod("convert_motifs", signature(motifs = "matrix"),
          definition = function(motifs, class) {
            motifs <- create_motif(motifs)
            convert_motifs(motifs, class = class)
          })

# @describeIn convert_motifs Convert non-\linkS4class{universalmotif} class motifs.
# @export
# setMethod("convert_motifs", signature(motifs = "ANY"),
          # definition = function(motifs, class) {

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
