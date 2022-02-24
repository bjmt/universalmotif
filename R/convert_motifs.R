#' Convert motif class.
#'
#' Allows for easy transfer of motif information between different classes as
#' defined by other Bioconductor packages. This function is also used by
#' nearly all other functions in this package, so any motifs of a compatible
#' class can be used without needing to be converted beforehand.
#'
#' @param motifs Single motif object or list. See details.
#' @param class `character(1)` Desired motif class. Input as
#'    'package-class'. If left empty, defaults to
#'    'universalmotif-universalmotif'. (See details.)
#'
#' @return Single motif object or list.
#'
#' @details
#' ## Input
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
#' ## Output
#' The following package-class combinations can be output:
#' * MotifDb-MotifList
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
#' * universalmotif-universalmotif (the default, no need to specify this)
#'
#' Note: MotifDb-MotifList output was a later addition to [convert_motifs()].
#' As a result, to stay consistent with previous behaviour most functions
#' will always convert MotifDb-MotifList objects to a list of `universalmotif`
#' motifs, even if other formats would be simply returned as is (e.g. for
#' other formats, [filter_motifs()] will return the input format; for
#' MotifDb-MotifList, a list of `universalmotif` objects will be returned).
#'
#' @examples
#' # Convert from universalmotif:
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' if (requireNamespace("motifStack", quietly = TRUE)) {
#'   jaspar.motifstack.pfm <- convert_motifs(jaspar, "motifStack-pfm")
#' }
#'
#' # Convert from another class to universalmotif:
#' if (requireNamespace("TFBSTools", quietly = TRUE)) {
#' library(TFBSTools)
#' data(MA0003.2)
#' motif <- convert_motifs(MA0003.2)
#'
#' # Convert from another class to another class
#' if (requireNamespace("PWMEnrich", quietly = TRUE)) {
#'   motif <- convert_motifs(MA0003.2, "PWMEnrich-PWM")
#' }
#'
#' # The 'convert_motifs' function is embedded in the rest of the universalmotif
#' # functions: non-universalmotif class motifs can be used
#' MA0003.2.trimmed <- trim_motifs(MA0003.2)
#' # Note: if the motif object going in has information that the
#' # 'universalmotif' class can't hold, it will be lost
#' }
#'
#' @references
#'
#' Bembom O (2018). *seqLogo: Sequence logos for DNA sequence
#' alignments*. R package version 1.46.0.
#'
#' Droit A, Gottardo R, Robertson G, Li L (2014). *rGADEM: de novo
#' motif discovery*. R package version 2.28.0.
#'
#' Mercier E, Gottardo R (2014). *MotIV: Motif Identification and
#' Validation*. R package version 1.36.0.
#'
#' Ou J, Wolfe SA, Brodsky MH, Zhu LJ (2018). “motifStack for the
#' analysis of transcription factor binding site evolution.” *Nature
#' Methods*, **15**, 8-9. doi: 10.1038/nmeth.4555.
#'
#' Shannon P, Richards M (2018). *MotifDb: An Annotated Collection of
#' Protein-DNA Binding Sequence Motifs*. R package version 1.22.0.
#'
#' Stojnic R, Diez D (2015). *PWMEnrich: PWM enrichment analysis*. R
#' package version 4.16.0.
#'
#' Tan G, Lenhard B (2016). “TFBSTools: an R/Bioconductor package for
#' transcription factor binding site analysis.” *Bioinformatics*,
#' **32**, 1555-1556. doi: 10.1093/bioinformatics/btw024.
#'
#' Yao Z (2012). *motifRG: A package for discriminative motif
#' discovery, designed for high throughput sequencing dataset*. R
#' package version 1.24.0.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @include universalmotif-class.R
#' @export
setGeneric("convert_motifs", function(motifs,
                                      class = "universalmotif-universalmotif")
           standardGeneric("convert_motifs"))

#' @describeIn convert_motifs Generate an error to remind users to run
#'   [to_list()] instead of using the column from [to_df()] directly.
#' @export
setMethod("convert_motifs", signature(motifs = "AsIs"),
  definition = function(motifs, class) {
    stop(wmsg(
      "If you are providing the `motif` column from the `data.frame` output of ",
      "to_df(), then please use to_list() instead. If this error ",
      "message does not apply to you, then remove the `AsIs` class attribute ",
      "and try again (`class(mylist) <- NULL`)."
    ))
})

#' @describeIn convert_motifs Convert a list of motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "list"),
          definition = function(motifs, class) {

  motifs <- unlist(motifs)
  if (!length(motifs)) stop("Input is an empty list")
  mot_classes <- unique(vapply(motifs, function(x) class(x), character(1)))
  if (length(mot_classes) == 1) {
    classin <- strsplit(class, "-", fixed = TRUE)[[1]][2]
    if (mot_classes == classin) return(motifs)
  }
  if (class == "MotifDb-MotifList") {
    motifs <- lapply(motifs, function(x) convert_motifs(x))
    motifs <- convert_to_motifdb_motiflist(motifs)
  } else {
    motifs <- lapply(motifs, function(x) convert_motifs(x, class = class))
  }

  motifs

})

#' @describeIn convert_motifs Convert a \linkS4class{universalmotif} object.
#' @export
setMethod("convert_motifs", signature(motifs = "universalmotif"),
          definition = function(motifs, class) {

  out_class <- strsplit(class, "-", fixed = TRUE)[[1]][2]
  out_class_pkg <- strsplit(class, "-", fixed = TRUE)[[1]][1]

  switch(out_class_pkg,
    "universalmotif" = {
      validObject_universalmotif(motifs)
      motifs
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
    "MotifDb" = {
      if (out_class == "MotifList") 
        convert_to_motifdb_motiflist(motifs)
      else
        stop("unknown 'class'")
    },
    stop("unknown 'class'")
  )

})

convert_to_tfbstools_matrix <- function(motifs, out_class) {
  motifs <- convert_type(motifs, "PCM")
  bkg <- motifs["bkg"][DNA_BASES]
  # names(bkg) <- DNA_BASES
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
                        pseudocounts = 1, bg = bkg)
      },
      "ICMatrix" = {
        motifs <- TFBSTools::toICM(motifs, pseudocounts = 1, bg = bkg)
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
  if (!"2" %in% names(motifs@multifreq)) {
    stop("cannot convert without filled multifreq slot")
  }
  if (motifs["alphabet"] != "DNA") stop("alphabet must be DNA")
  bkg <- motifs["bkg"][DNA_BASES]
  emission <- list(length = ncol(motifs@multifreq[["2"]]) + 1)
  emission[[1]] <- bkg
  transition <- matrix(rep(0, (ncol(motifs@multifreq$"2") + 1) *
                               ncol(motifs@multifreq$"2")),
                       nrow = ncol(motifs@multifreq$"2") + 1,
                       ncol = ncol(motifs@multifreq$"2"))
  for (i in seq_len(ncol(motifs@multifreq$"2"))) {
    emission[[i + 1]] <- motifs@multifreq$"2"[, i]
    transition[i, i] <- 1
  }
  names(emission) <- rep("state", length(emission))
  transition[nrow(transition), 1] <- 1
  transition[2, 1:2] <- c(0.95, 0.05)
  colnames(transition) <- seq_len(ncol(transition))
  rownames(transition) <- c(0, seq_len(ncol(transition)))
  if (length(motifs@altname) < 1) motifs@altname <- "Unknown"
  strand <- motifs@strand
  if (strand %in% c("+-", "-+")) strand <- "*"
  family <- motifs@family
  if (length(family) < 1) family <- "Unknown"
  if (requireNamespace("TFBSTools", quietly = TRUE)) {
    motifs <- TFBSTools::TFFMFirst(ID = motifs@altname,
                        name = motifs@name,
                        strand = strand, type = "First",
                        bg = bkg, matrixClass = family,
                        profileMatrix = motifs@motif,
                        tags = as.list(motifs@extrainfo),
                        emission = emission, transition = transition)
  } else {
    stop("package 'TFBSTools' is not installed") 
  }
  motifs
}

convert_to_seqlogo_pwm <- function(motifs) {
  motifs <- convert_type(motifs, "PPM")
  if (requireNamespace("seqLogo", quietly = TRUE)) {
    motifs <- seqLogo::makePWM(motifs@motif)
  } else {
    stop("'seqLogo' package not installed")
  }
  motifs
}

convert_to_motifdb_motiflist <- function(motifs) {
  if (!is.list(motifs)) motifs <- list(motifs)
  if (requireNamespace("MotifDb", quietly = TRUE)) {
    motiflist_class <- getClass("MotifList", where = "MotifDb")
    motifs <- convert_type(motifs, "PPM")
    m <- lapply(motifs, function(x) x["motif"])
    for (i in seq_along(m)) {
      colnames(m[[i]]) <- seq_len(ncol(m[[i]]))
    }
    m_names <- vapply(motifs, function(x) x["name"], character(1))
    m_altnames <- vapply(motifs, function(x) {
        x <- x["altname"]
        if (!length(x)) "<ID:unknown>"
        else x
      }, character(1))
    m_family <- vapply(motifs, function(x) {
        x <- x["family"]
        if (!length(x)) "<family:unknown>"
        else x
      }, character(1))
    m_organisms <- vapply(motifs, function(x) {
        x <- x["organism"]
        if (!length(x)) "<organism:unknown>"
        else x
      }, character(1))
    m_nsites <- as.numeric(sapply(motifs, function(x) {
        x <- x["nsites"]
        if (!length(x)) NA_real_
        else x
      }))
    if (any(duplicated(m_names))) {
      stop("MotifDb-MotifList cannot be created from motifs with duplicated 'name' slots",
        call. = FALSE)
    }
    names(m) <- m_names
    extras <- c("dataSource", "geneSymbol", "geneId", "proteinId", "proteinType",
      "bindingSequence", "bindingDomain", "experimentType", "pubmedID")
    meta <- DataFrame(
      providerName = m_names,
      providerId = m_altnames,
      dataSource = "<dataSource:unknown>",
      geneSymbol = NA_character_,
      geneId = NA_character_,
      proteinId = NA_character_,
      proteinIdType = NA_character_,
      organism = m_organisms,
      sequenceCount = m_nsites,
      bindingSequence = NA_character_,
      bindingDomain = NA_character_,
      tfFamily = m_family,
      experimentType = NA_character_,
      pubmedID = NA_character_
    )
    for (i in seq_along(m)) {
      m_ex_i <- motifs[[i]]["extrainfo"]
      for (j in seq_along(extras)) {
        m_ex_i_j <- m_ex_i[extras[j]]
        if (!is.na(m_ex_i_j)) {
          meta[[extras[j]]][i] <- m_ex_i_j
        }
      }
    }
    assoctab <- data.frame(
      motif = m_names,
      tf.gene = NA_character_,
      tf.ensg = NA_character_
    )
    motifs <- new(motiflist_class, listData = m,
      elementMetadata = meta,
      manuallyCuratedGeneMotifAssociationTable = assoctab)
  } else {
    stop("'MotifDb' package not installed")
  }
  motifs
}

convert_to_motifstack_pcm <- function(motifs) {
  if (requireNamespace("motifStack", quietly = TRUE)) {
    pcm_class <- getClass("pcm", where = "motifStack")
    motifs <- convert_type(motifs, "PCM")
    motifs <- new(pcm_class, mat = motifs@motif,
                  name = motifs@name,
                  alphabet = motifs@alphabet,
                  background = motifs@bkg[DNA_BASES])
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
    motifs <- new(pfm_class, mat = motifs@motif,
                  name = motifs@name,
                  alphabet = motifs@alphabet,
                  background = motifs@bkg[DNA_BASES])
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
    bio_mat <- matrix(as.integer(motifs@motif), byrow = FALSE,
                      nrow = 4)
    rownames(bio_mat) <- DNA_BASES
    bio_priors <- motifs@bkg[DNA_BASES]
    bio_mat <- PWMEnrich::PFMtoPWM(bio_mat, type = "log2probratio",
                        prior.params = bio_priors,
                        pseudo.count = motifs@pseudocount)
    motifs <- new(PWM_class, name = motifs["name"],
                  pfm = motifs@motif,
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
  bio_priors <- motifs@bkg[DNA_BASES]
  motifs <- PWM(x = bio_mat, type = "log2probratio",
                prior.params = bio_priors)
  motifs
}

convert_to_rgadem_motif <- function(motifs) {
  if (requireNamespace("rGADEM", quietly = FALSE)) {
    motifs <- convert_type(motifs, "PPM")
    rGADEM_motif_class <- getClass("motif", where = "rGADEM")
    motifs <- new(rGADEM_motif_class, pwm = motifs@motif,
                  name = motifs@name,
                  consensus = motifs@consensus)
  } else {
    stop("'rGADEM' package not installed")
  }
  motifs
}

#' @describeIn convert_motifs Convert MotifList motifs. (\pkg{MotifDb})
#' @export
setMethod("convert_motifs", signature(motifs = "MotifList"),
          definition = function(motifs, class) {

  x <- motifs

  altname <- x@elementMetadata@listData$providerName
  name <- x@elementMetadata@listData$geneSymbol
  family <- x@elementMetadata@listData$tfFamily
  organism <- x@elementMetadata@listData$organism
  motif <- x@listData
  nsites <- as.numeric(x@elementMetadata@listData$sequenceCount)
  dataSource <- x@elementMetadata@listData$dataSource

  motifdb_fun <- function(i) {

    mot <- universalmotif_cpp(altname = altname[i], name = name[i],
                              family = family[i], organism = organism[i],
                              motif = motif[[i]], alphabet = "DNA",
                              type = "PPM", nsites = nsites[i],
                              extrainfo = c(dataSource = dataSource[i]))

    validObject_universalmotif(mot)
    mot

  }

  motifs_out <- lapply(seq_len(length(x)), motifdb_fun)
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
    mot <- universalmotif_cpp(name = motifs@name, altname = motifs@ID,
                   strand = motifs@strand, bkg = motifs@bg,
                   motif = TFBSTools::getPosProb(motifs))
    mot@multifreq <- list("2" = difreq)
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
                               family = paste0(motifs@tags$family, collapse = " / "),
                               nsites = nsites,
                               organism = paste0(motifs@tags$species,
                                                 collapse = "/"),
                               motif = motifs@profileMatrix,
                               alphabet = alphabet, type = "PCM",
                               bkg = motifs@bg,
                               strand = collapse_cpp(motifs@strand),
                               extrainfo = extrainfo)
  validObject_universalmotif(motifs)
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
                               family = paste0(motifs@tags$family, collapse = " / "),
                               organism = paste0(motifs@tags$species,
                                                 collapse = "/"),
                               motif = motifs@profileMatrix,
                               alphabet = alphabet, type = "PWM",
                               bkg = motifs@bg,
                               strand = collapse_cpp(motifs@strand),
                               extrainfo = extrainfo)
  validObject_universalmotif(motifs)
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
                           family = paste0(motifs@tags$family, collapse = " / "),
                           organism = paste0(motifs@tags$species,
                                             collapse = "/"),
                           motif = motifs@profileMatrix,
                           alphabet = alphabet, type = "ICM",
                           bkg = motifs@bg,
                           strand = collapse_cpp(motifs@strand),
                           extrainfo = extrainfo)
  validObject_universalmotif(motifs)
  convert_motifs(motifs, class = class)

})

#' @describeIn convert_motifs Convert XMatrixList motifs. (\pkg{TFBSTools})
#' @export
setMethod("convert_motifs", signature(motifs = "XMatrixList"),
          definition = function(motifs, class) {

  motif_num <- length(motifs@listData)
  motifs_out <- lapply(seq_len(motif_num),
                       function(i) motifs@listData[[i]])
  motif_names <- unlist(lapply(seq_len(motif_num),
                               function(i) motifs@listData[[i]]@name))
  motifs <- lapply(motifs_out, convert_motifs)
  convert_motifs(motifs, class = class)

})

#' @describeIn convert_motifs Convert pwm motifs. (\pkg{seqLogo})
#' @export
setMethod("convert_motifs", signature(motifs = "pwm"),
          definition = function(motifs, class) {
  motifs <- universalmotif_cpp(motif = motifs@pwm, type = "PPM",
                           alphabet = motifs@alphabet)
  validObject_universalmotif(motifs)
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
  validObject_universalmotif(motifs)
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
  validObject_universalmotif(motifs)
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
  validObject_universalmotif(motifs)
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
  validObject_universalmotif(motifs)
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
