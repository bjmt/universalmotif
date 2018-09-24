#' @param x [universalmotif-class] Motif.
#' @param i `character` Slot.
#' @include universalmotif-class.R
#' @rdname universalmotif-class
#' @aliases [,universalmotif-method
setMethod("[", "universalmotif", function(x, i) {

  if (missing(i)) {
    i <- c("name", "altname", "family", "organism", "motif", "alphabet", "type",
               "icscore", "nsites", "pseudocount", "bkg", "bkgsites",
               "consensus", "strand", "pval", "qval", "eval",
               "multifreq", "extrainfo")
  }

  if (all(i == "motif")) return(x@motif)
  if (all(i == "multifreq")) return(x@multifreq)
  
  return_list <- lapply(i, function(y) slot(x, y))
  names(return_list) <- i
  if ("motif" %in% names(return_list)) {
    return_list$motif <- x@motif
  }
  if ("multifreq" %in% names(return_list)) {
    return_list$multifreq <- x@multifreq
  }
 
  if (length(return_list) <= 1) {
    return_list <- unlist(return_list)
    if (i == "extrainfo") {
      names(return_list) <- gsub("extrainfo.", "", names(return_list))
      return(return_list)
    }
  }
  return(return_list)

})

#' @param value Object to replace slot with.
#' @rdname universalmotif-class
#' @aliases [<-,universalmotif-method
setMethod("[<-", "universalmotif", function(x, i, value) {
  if (i == "icscore") stop("'icscore' is generated automatically")
  if (i == "multifreq") stop("please use 'add_multifreq'")
  if (i == "consensus" && x@alphabet %in% c("DNA", "RNA", "AA")) {
    stop("consensus string for ", x@alphabet, " motifs is generated automatically")
  }
  slot(x, i) <- value
  # msg <- validObject_universalmotif(x)
  # if (length(msg) > 0) stop(msg) else x
  if (validObject(x)) x
})

#' @param .Object [universalmotif-class] Final motif.
#' @param name `character` Motif name.
#' @param altname `character` Alternate motif name.
#' @param family `character` Transcription factor family.
#' @param organism `character` Species of origin.
#' @param motif `matrix` Each column represents a position in the motif.
#' @param alphabet `character(1)` One of `c('DNA', 'RNA', 'AA', 'custom')`,
#'    or a combined string representing the letters.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#' @param icscore `numeric` Total information content. Automatically generated.
#' @param nsites `numeric` Number of sites the motif was constructed from.
#' @param pseudocount `pseudocount(1)` Correction to be applied to prevent `-Inf`
#'   from apearing in PWM matrices.
#' @param bkg `numeric` Must sum to 1 and be equal in length to the alphabet
#'            length. If missing, assumes a uniform background.
#' @param bkgsites `numeric` Total number of sites used to find the motif.
#' @param consensus `character` Consensus string. Automatically generated for 
#'    'DNA', 'RNA', and 'AA' alphabets.
#' @param strand `character(1)` Whether the motif is specific to a certain strand.
#' @param pval `numeric` P-value associated with motif.
#' @param qval `numeric` Adjusted P-value associated with motif.
#' @param eval `numeric` E-value associated with motif.
#' @param multifreq `list` See [add_multifreq()].
#' @param extrainfo `character` Any other extra information, represented as
#'    a named character vector.
#' @name universalmotif
#' @rdname universalmotif-class
#' @aliases initialize,universalmotif-method
setMethod("initialize", signature = "universalmotif",
          definition = function(.Object, name, altname, family, 
                                organism, motif,
                                alphabet = "DNA", type, icscore, nsites,
                                pseudocount = 0.8, bkg, bkgsites,
                                consensus, strand = "+-", pval,
                                qval, eval, multifreq, extrainfo) {
            
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

            if (!alphabet %in% c("DNA", "RNA", "AA", "custom")) {
              alphabet <- strsplit(alphabet, "")[[1]]
              alphabet <- paste(sort(alphabet), collapse = "")
            }
            .Object@alphabet <- alphabet

            if (missing(alphabet) || length(alphabet) == 0 ||
                is.na(alphabet)) {
              if (nrow(motif) == 4) alphabet <- "DNA" else {
                if (nrow(motif) == 20) {
                  alphabet <- "AA"
                } else alphabet <- "custom"
              }
            }

            if (missing(type) || length(type) == 0 || is.na(type)) {
              if (all(motif >= 1 | motif == 0)) {
                type <- "PCM" 
              } else if ((all(colSums(motif) > 0.99 &
                              colSums(motif) < 1.01)) &&
                         all(motif >= 0)) {
                type <- "PPM"
              } else if (all(motif >= 0)) {
                type <- "ICM"
              } else type <- "PWM"
            }
            .Object@type <- type

            if (missing(nsites) || length(nsites) == 0  || is.na(nsites)) {
              if (type == "PCM") {
                nsites <- sum(motif[, 1])
              } else nsites <- numeric(0)
            } else if (type == "PCM" && any(colSums(motif) != nsites)) {
              for (i in seq_len(ncol(motif))) {
                motif[, i] <- motif[, i] / sum(motif[, i])
              }
              motif <- apply(motif, 2, ppm_to_pcmC, nsites = nsites)
              .Object@motif <- motif
            }
            .Object@nsites <- nsites

            if (is.null(rownames(.Object@motif))) {
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
            } else {
              .Object@motif <- .Object@motif[order(rownames(.Object@motif)), ]
            }
            if (missing(bkg) || length(bkg) == 0 || is.na(bkg)) {
              if (alphabet %in% c("DNA", "RNA")) {
                bkg <- rep(0.25, 4)
              } else if (alphabet == "AA") {
                bkg <- rep(0.05, 20)
              } else {
                bkg <- rep(1 / nrow(motif), nrow(motif))
              }
            }
            if (length(bkg) != nrow(.Object@motif)) {
              stop("bkg vector length should be ", nrow(.Object@motif))
            }
            .Object@bkg <- bkg

            if (missing(icscore) || length(icscore) == 0 || is.na(icscore)) {
              icscores <- apply(motif, 2, position_icscoreC,
                                bkg = bkg, type = type,
                                pseudocount = pseudocount,
                                nsites = ifelse(length(nsites) == 0, 100, nsites))
              icscore <- sum(icscores)
            }
            .Object@icscore <- icscore

            .Object@pseudocount <- pseudocount

            if (missing(bkgsites) || length(bkgsites) == 0 ||
                is.na(bkgsites)) {
              bkgsites <- numeric(0)
            }
            .Object@bkgsites <- bkgsites

            if (missing(consensus) || length(consensus) == 0 ||
              is.na(consensus)) {
              if (.Object@alphabet %in% c("DNA", "RNA")) {
                consensus <- apply(motif, 2, get_consensusC, alphabet = alphabet,
                                   type = type, pseudocount = pseudocount)
                .Object@consensus <- paste(consensus, collapse = "")
              } else if (.Object@alphabet == "AA") {
                consensus <- apply(motif, 2, get_consensusAAC, type = type,
                                   pseudocount = pseudocount)
                .Object@consensus <- consensus
              } else .Object@consensus <- character(0)
            }

            if (length(.Object@consensus) > 0) {
              names(consensus) <- NULL
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

            if (!missing(multifreq)) .Object@multifreq <- multifreq

            if (missing(extrainfo) || length(extrainfo) == 0 || 
                is.na(extrainfo)) {
              extrainfo <- character(0)
            }
            .Object@extrainfo <- extrainfo

            validObject(.Object)
            .Object

          })

#' @param object [universalmotif-class] Motif.
#' @rdname universalmotif-class
#' @aliases show,universalmotif-method
setMethod("show", signature = "universalmotif",
          definition = function(object) {
            name <- object@name
            if (nchar(name) > 40) {
              name <- paste0(substr(name, 1, 40), "...")
            }
            cat("\n       Motif name:   ", name, "\n", sep = "")
            if (length(object@altname) > 0) {
              altname <- object@altname
              if (nchar(altname) > 40) {
                altname <- paste0(substr(name, 1, 40), "...")
              }
              cat("   Alternate name:   ", altname, "\n", sep = "")
            }
            if (length(object@family)) {
              family <- object@family
              if (nchar(family) > 40) {
                family <- paste0(substr(family, 1, 40), "...")
              }
              cat("           Family:   ", family, "\n", sep = "")
            }
            if (length(object@organism)) {
              organism <- object@organism
              if (nchar(organism) > 40) {
                organism <- paste0(substr(organism, 1, 40), "...")
              } 
              cat("         Organism:   ", organism, "\n", sep = "")
            }
            alphabet <- object@alphabet
            if (nchar(alphabet) > 40) {
              alphabet <- paste0(substr(alphabet, 1, 40), "...")
            }
            cat("         Alphabet:   ", alphabet, "\n", sep = "")
            cat("             Type:   ", object@type, "\n", sep = "")
            if (object@alphabet %in% c("DNA", "RNA")) {
              cat("          Strands:   ", object@strand, "\n", sep = "")
            }
            cat("         Total IC:   ", object@icscore, "\n", sep = "")
            if (length(object@consensus) > 0) {
              consensus <- object@consensus
              if (nchar(consensus) > 40) {
                consensus <- paste0(substr(consensus, 1, 40), "...")
              }
              cat("        Consensus:   ", consensus, "\n", sep = "")
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
            if (length(object@multifreq) > 0) {
              toprint <- paste(names(object@multifreq), collapse = ", ")
              cat("   k-letter freqs:  ", toprint, "\n")
            }
            if (length(object@extrainfo) > 0 ) {
              extrainfo <- object@extrainfo
              cat("       Extra info:   ")
              for (i in seq_len(length(extrainfo))) {
                if (!is.null(names(extrainfo))) {
                  to_show <- paste0(names(extrainfo[i]), ": ",
                                    extrainfo[i])
                } else to_show <- extrainfo[i]
                if (nchar(to_show) > 40) {
                  to_show <- paste0(substr(to_show, 1, 40), "...")
                }
                if (i == 1) {cat(to_show, "\n"); next}
                cat("                     ", to_show,
                    "\n", sep = "")
              }
            }
            cat("\n")
            print(round(object@motif * 100) / 100)
            invisible(NULL)
          })

# as.matrix
# as.character
# names

#' @rdname universalmotif-class
#' @aliases as.data.frame,universalmotif-method
setMethod("as.data.frame", signature(x = "universalmotif"),
          definition = function(x) {
            name <- x@name
            altname <- x@altname
            if (length(altname) == 0) altname <- as.character(NA)
            family <- x@family
            if (length(family) == 0) family <- as.character(NA)
            organism <- x@organism
            if (length(organism) == 0) organism <- as.character(NA)
            alphabet <- x@alphabet
            icscore <- x@icscore
            nsites <- x@nsites
            if (length(nsites) == 0) nsites <- as.numeric(NA)
            bkgsites <- x@bkgsites
            if (length(bkgsites) == 0) bkgsites <- as.numeric(NA)
            consensus <- x@consensus
            if (length(consensus) == 0) consensus <- as.character(NA)
            strand <- x@strand
            pval <- x@pval
            if (length(pval) == 0) pval <- as.numeric(NA)
            qval <- x@qval
            if (length(qval) == 0) qval <- as.numeric(NA)
            eval <- x@eval
            if (length(eval) == 0) eval <- as.numeric(NA)
            df <- data.frame(name = name, altname = altname, family = family,
                             organism = organism, alphabet = alphabet,
                             icscore = icscore, nsites = nsites,
                             bkgsites = bkgsites,
                             consensus = consensus, strand = strand,
                             pval = pval, qval = qval, eval = eval,
                             stringsAsFactors = FALSE)
            rownames(df) <- NULL
            df
          })

#' @param select `numeric` Columns to keep.
#' @rdname universalmotif-class
#' @aliases subset,universalmotif-method
setMethod("subset", signature(x = "universalmotif"),
          definition = function(x, select) {
            mot <- x["motif"][, select]
            motif <- universalmotif_cpp(motif = mot, name = x@name,
                                        altname = x@altname, family = x@family,
                                        organism = x@organism,
                                        alphabet = x@alphabet, nsites = x@nsites,
                                        pseudocount = x@pseudocount, bkg = x@bkg,
                                        bkgsites = x@bkgsites, strand = x@strand,
                                        pval = x@pval, qval = x@qval, eval = x@eval,
                                        extrainfo = x@extrainfo)
            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif
          })

#' @rdname universalmotif-class
#' @aliases normalize,universalmotif-method
setMethod("normalize", signature(object = "universalmotif"),
          definition = function(object) {
            type <- object["type"]
            pseudo <- object["pseudocount"]
            if (pseudo == 0) pseudo <- 1
            object <- convert_type(object, "PCM")
            convert_type(object, type, pseudocount = pseudo)
          })

#' @rdname universalmotif-class
#' @aliases rowMeans,universalmotif-method
setMethod("rowMeans", signature(x = "universalmotif"),
          definition = function(x) rowMeans(x["motif"]))

#' @rdname universalmotif-class
#' @aliases colMeans,universalmotif-method
setMethod("colMeans", signature(x = "universalmotif"),
          definition = function(x) colMeans(x["motif"]))

#' @rdname universalmotif-class
#' @aliases colSums,universalmotif-method
setMethod("colSums", signature(x = "universalmotif"),
          definition = function(x) colSums(x["motif"]))

#' @rdname universalmotif-class
#' @aliases rowSums,universalmotif-method
setMethod("rowSums", signature(x = "universalmotif"),
          definition = function(x) rowSums(x["motif"]))

#' @rdname universalmotif-class
#' @aliases nrow,universalmotif-method
setMethod("nrow", signature = "universalmotif",
          definition = function(x) nrow(x["motif"]))

#' @rdname universalmotif-class
#' @aliases ncol,universalmotif-method
setMethod("ncol", signature = "universalmotif",
          definition = function(x) ncol(x["motif"]))

#' @rdname universalmotif-class
#' @aliases colnames,universalmotif-method
setMethod("colnames", signature(x = "universalmotif"),
          definition = function(x) colnames(x["motif"]))

#' @rdname universalmotif-class
#' @aliases rownames,universalmotif-method
setMethod("rownames", signature(x = "universalmotif"),
          definition = function(x) rownames(x["motif"]))

#' @param ... [universalmotif-class] Motifs.
#' @param deparse.level Unused.
#' @rdname universalmotif-class
#' @aliases cbind,universalmotif-method
setMethod("cbind", signature = "universalmotif",
          definition = function(..., deparse.level = 0) {

            mots <- list(...)
            if (length(mots) == 1) return(mots[[1]])
            mot.alphabet <- unique(vapply(mots, function(x) x["alphabet"], character(1)))
            if (length(mot.alphabet) != 1)
              stop("all motifs must have the same alphabet")
            mots <- convert_type(mots, "PPM")
            mot.name <- vapply(mots, function(x) x["name"], character(1))
            mot.altname <- do.call(c, lapply(mots, function(x) x["altname"]))
            mot.family <- do.call(c, lapply(mots, function(x) x["family"]))
            mot.organism <- do.call(c, lapply(mots, function(x) x["organism"]))
            mot.motif <- lapply(mots, function(x) x["motif"])
            mot.nsites <- do.call(c, lapply(mots, function(x) x["nsites"]))
            mot.bkg <- lapply(mots, function(x) x["bkg"])
            mot.bkgsites <- lapply(mots, function(x) x["bkgsites"])
            mot.bkgsites <- do.call(c, mot.bkgsites)
            mot.strand <- unique(vapply(mots, function(x) x["strand"], character(1)))
            mot.pval <- do.call(c, lapply(mots, function(x) x["pval"]))
            mot.qval <- do.call(c, lapply(mots, function(x) x["qval"]))
            mot.eval <- do.call(c, lapply(mots, function(x) x["eval"]))
            mot.extrainfo <- lapply(mots, function(x) x["extrainfo"])
            mot.extrainfo <- do.call(c, mot.extrainfo)

            mot.motif <- do.call(cbind, mot.motif)
            mot.name <- paste(mot.name, collapse = "/")
            if (length(mot.altname) > 0) mot.altname <- paste(mot.altname, collapse = "/")
            if (length(mot.family) > 0) mot.family <- paste(mot.family, collapse = "/")
            if (length(mot.organism) > 0) mot.organism <- paste(mot.organism, collapse = "/")
            if (length(mot.nsites) > 1) mot.nsites <- max(mot.nsites)
            mot.bkg <- colMeans(do.call(rbind, mot.bkg))
            if (length(mot.bkgsites) > 1) mot.bkgsites <- max(mot.bkgsites)
            if (length(mot.strand) > 1) mot.strand <- "+-"
            if (length(mot.pval) > 1) mot.pval <- min(mot.pval)
            if (length(mot.qval) > 1) mot.qval <- min(mot.qval)
            if (length(mot.eval) > 1) mot.eval <- min(mot.eval)

            motif <- universalmotif_cpp(motif = mot.motif, name = mot.name,
                                        altname = mot.altname, family = mot.family,
                                        organism = mot.organism, nsites = mot.nsites,
                                        alphabet = mot.alphabet, bkg = mot.bkg,
                                        bkgsites = mot.bkgsites, strand = mot.strand,
                                        pval = mot.pval, qval = mot.qval, eval = mot.eval,
                                        extrainfo = mot.extrainfo)
            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif

          })
