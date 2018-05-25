#' @describeIn universalmotif Accessor function.
#' @include universalmotif-class.R
setMethod("[", "universalmotif", function(x, i, drop = "missing") {

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

#' @describeIn universalmotif Replacement method.
setMethod("[<-", "universalmotif", function(x, i, value) {
  slot(x, i) <- value
  if (validObject(x)) return(x)
})

#' @describeIn universalmotif Initialize method.
setMethod("initialize", signature = "universalmotif",
          definition = function(.Object, name, altname, family, 
                                organism, motif,
                                alphabet = "DNA", type, icscore, nsites,
                                pseudoweight = 0.8, bkg, bkgsites,
                                consensus, strand = "+-", pval,
                                qval, eval, extrainfo) {
            
            if (missing(name)) name <- "new motif"
            .Object@name <- name

            if (missing(altname)) altname <- character(0)
            .Object@altname <- altname

            if(missing(family)) family <- character(0)
            .Object@family <- family

            if (missing(organism)) organism <- character(0)
            .Object@organism <- organism

            if (missing(motif)) stop("missing motif matrix")
            .Object@motif <- motif

            .Object@alphabet <- alphabet
            
            if (missing(alphabet)) {
              if (nrow(motif) == 4) alphabet <- "DNA" else {
                if (nrow(motif) == 20) {
                  alphabet <- "AA"
                } else alphabet <- "custom"
              }
            }

            if (alphabet == "DNA") {
              rownames(.Object@motif) <- c("A", "C", "G", "T")
            } else if (alphabet == "RNA") {
              rownames(.Object@motif)  <- c("A", "C", "G", "U")
            } else if (alphabet == "AA") {
              rownames(.Object@motif) <- c("A", "C", "D", "E", "F", "G",
                                          "H", "I", "K", "L", "M", "N", "P",
                                          "Q", "R", "S", "T", "V", "W", "Y")
            }

            if (missing(type)) {
              if (all(colSums(motif) >= 2)) type <- "PCM" else {
                if (all(motif >= 0)) type <- "PPM" else {
                  type <- "PWM"
                  message("assumed 'type' as being PWM")
                }
              }
            }
            .Object@type <- type

            if (missing(nsites) || length(nsites) == 0  || is.na(nsites)) {
              if (type == "PCM") {
                nsites <- sum(motif[, 1])
              } else nsites <- numeric(0)
            }
            .Object@nsites <- nsites

            if (missing(bkg)) {
              if (alphabet %in% c("DNA", "RNA")) {
                bkg <- rep(0.25, 4)
              } else if (alphabet == "AA") {
                bkg <- rep(0.05, 20)
              } else if (alphabet == "custom") {
                stop("motif 'bkg' required for 'custom' alphabet")
              }
            }
            .Object@bkg <- bkg

            if (missing(icscore)) {
              icscores <- apply(motif, 2, position_icscore,
                                bkg = bkg, type = type,
                                pseudoweight = pseudoweight, nsites = nsites)
              icscore <- sum(icscores)
            }
            .Object@icscore <- icscore

            .Object@pseudoweight <- pseudoweight

            if (missing(bkgsites)) bkgsites <- numeric(0)
            .Object@bkgsites <- bkgsites

            if (missing(consensus)) {
              if (.Object@alphabet %in% c("DNA", "RNA")) {
                # consensus <- consensusString(.Object@motif, threshold = 0.25,
                                             # ambiguityMap = IUPAC_CODE_MAP)
                # .Object@consensus <- consensus
                # consensus <- strsplit(consensus, "")[[1]]
                consensus <- apply(motif, 2, get_consensus, alphabet = alphabet,
                                   type = type, pseudoweight = pseudoweight)
                .Object@consensus <- paste(consensus, collapse = "")
              } else if (.Object@alphabet == "AA") {
                consensus <- apply(motif, 2, get_consensusAA)
                .Object@consensus <- consensus
              } else .Object@consensus <- character(0)
            }

            if (length(.Object@consensus) > 0) colnames(.Object@motif) <- consensus

            .Object@strand <- strand

            if (missing(pval)) pval <- numeric(0)
            .Object@pval <- pval

            if (missing(qval)) qval <- numeric(0)
            .Object@qval <- qval

            if (missing(eval)) eval <- numeric(0)
            .Object@eval <- eval

            if (missing(extrainfo)) extrainfo <- character(0)
            .Object@extrainfo <- extrainfo

            .Object

          })

#' @describeIn universalmotif Show method.
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
            cat("             Type:   ", object@type, "\n", sep = "")
            cat("          Strands:   ", object@strand, "\n", sep = "")
            cat("         Total IC:   ", object@icscore, "\n", sep = "")
            cat("        Consensus:   ", object@consensus, "\n", sep = "")
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
