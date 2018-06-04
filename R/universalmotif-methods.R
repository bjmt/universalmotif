#' @param x universalmotif object.
#' @param i Slot.
#' @include universalmotif-class.R
#' @rdname universalmotif-class
#' @aliases [,universalmotif-method
setMethod("[", "universalmotif", function(x, i) {

  if (missing(i)) {
    i <- c("name", "altname", "family", "organism", "motif", "alphabet", "type",
               "icscore", "nsites", "pseudoweight", "bkg", "bkgsites",
               "consensus", "strand", "pval", "qval", "eval", "extrainfo")
  }

  if (all(i == "motif")) return(x@motif)
  
  return_list <- lapply(i, function(y) slot(x, y))
  names(return_list) <- i
  if ("motif" %in% names(return_list)) {
    return_list$motif <- x@motif
  }

  if (length(return_list) <= 1) return(unlist(return_list))
  return(return_list)

})

#' @param value Object to replace slot with.
#' @rdname universalmotif-class
#' @aliases [<-,universalmotif-method
setMethod("[<-", "universalmotif", function(x, i, value) {
  slot(x, i) <- value
  if (validObject(x)) return(x)
})

#' @param .Object universalmotif object.
#' @param name Character. Motif name.
#' @param altname Character. Alternate motif name.
#' @param family Character. Transcription factor family.
#' @param organism Character. Species of origin.
#' @param motif Matrix.
#' @param alphabet Character. 'DNA', 'RNA', 'AA', 'custom', or a combined 
#'                 string representing the letters.
#' @param type Character. 'PCM', 'PPM', 'PWM', or 'ICM'.
#' @param icscore Numeric. Total information content. Automatically generated.
#' @param nsites Numeric. Number of sites the motif was constructed from.
#' @param pseudoweight Numeric. Correction to be applied to prevent \code{-Inf}
#'                     from apearing in PWM matrices.
#' @param bkg Numeric. Must sum to 1 and be equal in length to the alphabet
#'            length.
#' @param bkgsites Numeric. Total number of sites used to find the motif.
#' @param consensus Character. Consensus string. Automatically generated for 
#'                  'DNA', 'RNA', and 'AA' alphabets.
#' @param strand Character. Whether the motif is specific to a certain strand.
#' @param pval Numeric. P-value associated with motif.
#' @param qval Numeric. Adjusted P-value associated with motif.
#' @param eval Numeric. E-value associated with motif.
#' @param extrainfo Character. Any other extra information, represented as
#'                  a named character vector.
#' @name universalmotif
#' @rdname universalmotif-class
#' @aliases initialize,universalmotif-method
setMethod("initialize", signature = "universalmotif",
          definition = function(.Object, name, altname, family, 
                                organism, motif,
                                alphabet = "DNA", type, icscore, nsites,
                                pseudoweight = 0.8, bkg, bkgsites,
                                consensus, strand = "+-", pval,
                                qval, eval, extrainfo) {
            
            if (missing(name) || length(name) == 0 || is.na(name)) {
              name <- "new motif"
            }
            .Object@name <- name

            if (missing(altname) || length(altname) == 0 ||
                is.na(altname)) {
              altname <- character(0)
            }
            .Object@altname <- altname

            if(missing(family) || length(family) == 0 ||
               is.na(family)) {
              family <- character(0)
            }
            .Object@family <- family

            if (missing(organism) || length(organism) == 0 ||
                is.na(organism)) {
              organism <- character(0)
            }
            .Object@organism <- organism

            if (missing(motif)) stop("missing motif matrix")

            .Object@motif <- motif
            # if (!missing(type)) {
              # if (type == "PCM") {
                # .Object@motif <- round(motif)
              # }
            # }

            .Object@alphabet <- alphabet
            
            if (missing(alphabet) || length(alphabet) == 0 ||
                is.na(alphabet)) {
              if (nrow(motif) == 4) alphabet <- "DNA" else {
                if (nrow(motif) == 20) {
                  alphabet <- "AA"
                } else alphabet <- "custom"
              }
            }

            if (alphabet == "DNA") {
              rownames(.Object@motif) <- DNA_BASES
            } else if (alphabet == "RNA") {
              rownames(.Object@motif)  <- RNA_BASES
            } else if (alphabet == "AA") {
              rownames(.Object@motif) <- AA_STANDARD
            } else if (length(strsplit(alphabet, "")[[1]]) == nrow(motif) &&
                       alphabet != "custom") {
              rownames(.Object@motif) <- strsplit(alphabet, "")[[1]]
            } 

            if (missing(type) || length(type) == 0 || is.na(type)) {
              if (all(motif >= 1 || motif == 0)) {
                type <- "PCM" 
              } else if ((all(colSums(motif) > 0.99 ||
                              colSums(motif) < 1.01)) &&
                         all(motif >= 0)) {
                type <- "PPM"
              } else if (all(colSums(motif) >= 0)) {
                type <- "ICM"
              } else type <- "PWM"
            }
            .Object@type <- type

            if (missing(nsites) || length(nsites) == 0  || is.na(nsites)) {
              if (type == "PCM") {
                nsites <- sum(motif[, 1])
              } else nsites <- numeric(0)
            }
            .Object@nsites <- nsites

            if (missing(bkg) || length(bkg) == 0 || is.na(bkg)) {
              if (alphabet %in% c("DNA", "RNA")) {
                bkg <- rep(0.25, 4)
              } else if (alphabet == "AA") {
                bkg <- rep(0.05, 20)
              } else {
                bkg <- rep(1 / nrow(motif), nrow(motif))
              }
            }
            .Object@bkg <- bkg

            if (missing(icscore) || length(icscore) == 0 || is.na(icscore)) {
              icscores <- apply(motif, 2, position_icscore,
                                bkg = bkg, type = type,
                                pseudoweight = pseudoweight, nsites = nsites)
              icscore <- sum(icscores)
            }
            .Object@icscore <- icscore

            .Object@pseudoweight <- pseudoweight

            if (missing(bkgsites) || length(bkgsites) == 0 ||
                is.na(bkgsites)) {
              bkgsites <- numeric(0)
            }
            .Object@bkgsites <- bkgsites

            if (missing(consensus) || length(consensus) == 0 ||
              is.na(consensus)) {
              if (.Object@alphabet %in% c("DNA", "RNA")) {
                consensus <- apply(motif, 2, get_consensus, alphabet = alphabet,
                                   type = type, pseudoweight = pseudoweight)
                .Object@consensus <- paste(consensus, collapse = "")
              } else if (.Object@alphabet == "AA") {
                consensus <- apply(motif, 2, get_consensusAA, type = type,
                                   pseudoweight = pseudoweight)
                .Object@consensus <- consensus
              } else .Object@consensus <- character(0)
            }

            if (length(.Object@consensus) > 0) {
              colnames(.Object@motif) <- consensus
            }

            .Object@strand <- strand

            if (missing(pval) || length(pval) == 0 || is.na(pval)) {
              pval <- numeric(0)
            }
            .Object@pval <- pval

            if (missing(qval) || length(qval) == 0 || is.na(qval)) {
              qval <- numeric(0)
            }
            .Object@qval <- qval

            if (missing(eval) || length(eval) == 0 || is.na(eval)) {
              eval <- numeric(0)
            }
            .Object@eval <- eval

            if (missing(extrainfo) || length(extrainfo) == 0 || 
                is.na(extrainfo)) {
              extrainfo <- character(0)
            }
            .Object@extrainfo <- extrainfo

            validObject(.Object)
            .Object

          })

#' @param object \linkS4class{universalmotif} object.
#' @rdname universalmotif-class
#' @aliases show,universalmotif-method
setMethod("show", signature = "universalmotif",
          definition = function(object) {
            cat("\n       Motif name:   ", object@name, "\n", sep = "")
            if (length(object@altname) > 0) {
              cat("   Alternate name:   ", object@altname, "\n", sep = "")
            }
            if (length(object@family)) {
              cat("           Family:   ", object@family, "\n", sep = "")
            }
            if (length(object@organism)) {
              cat("         Organism:   ", object@organism, "\n", sep = "")
            }
            cat("         Alphabet:   ", object@alphabet, "\n", sep = "")
            cat("             Type:   ", object@type, "\n", sep = "")
            cat("          Strands:   ", object@strand, "\n", sep = "")
            cat("         Total IC:   ", object@icscore, "\n", sep = "")
            if (length(object@consensus) > 0) {
              cat("        Consensus:   ", object@consensus, "\n", sep = "")
            }
            if (length(object@nsites) > 0) {
              cat("     Target sites:   ", object@nsites, "\n", sep = "")
            }
            if (length(object@bkgsites) > 0) {
              cat(" Background sites:   ", object@bkgsites, "\n", sep = "")
            }
            if (length(object@pval) > 0) {
              cat("          P-value:   ", object@pval, "\n", sep = "")
            }
            if (length(object@qval) > 0) {
              cat("          Q-value:   ", object@qval, "\n", sep = "")
            }
            if (length(object@eval) > 0) {
              cat("          E-value:   ", object@eval, "\n", sep = "")
            }
            if (length(object@extrainfo) > 0 ) {
              cat("       Extra info:   ")
              for (i in seq_len(length(object@extrainfo))) {
                if (!is.null(names(object@extrainfo))) {
                  to_show <- paste0(names(object@extrainfo[i]), ": ",
                                    object@extrainfo[i])
                } else to_show <- object@extrainfo[i]
                if (i == 1) {cat(to_show, "\n"); next}
                cat("                     ", to_show,
                    "\n", sep = "")
              }
            }
            cat("\n")
            print(object@motif, digits = 3)
            invisible(NULL)
          })
