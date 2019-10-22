setMethod("initialize", signature = "universalmotif_gapped",
          definition = function(.Object) {
  .Object@isgapped <- FALSE
  .Object
})

#' @param x [universalmotif-class] Motif.
#' @param i `character` Slot.
#' @include universalmotif-class.R
#' @rdname universalmotif-class
#' @aliases [,universalmotif-method
setMethod("[", "universalmotif", function(x, i) {

  validObject_universalmotif(x)

  if (missing(i)) return(universalmotif_to_list(x))

  if (all(i == "motif")) return(x@motif)
  if (all(i == "multifreq")) return(x@multifreq)
  if (all(i == "bkg")) return(x@bkg)

  return_list <- lapply(i, function(y) slot(x, y))

  if (length(return_list) <= 1) {

    return_list <- unlist(return_list)

  } else {

    names(return_list) <- i

    if ("motif" %in% names(return_list))
      return_list$motif <- x@motif
    if ("multifreq" %in% names(return_list))
      return_list$multifreq <- x@multifreq
    if ("bkg" %in% names(return_list))
      return_list$bkg <- x@bkg

  }

  return_list

})

#' @param value Object to replace slot with.
#' @rdname universalmotif-class
#' @aliases [<-,universalmotif-method
setMethod("[<-", "universalmotif", function(x, i, value) {

  if (i %in% c("icscore", "multifreq", "consensus", "motif"))
    stop(wmsg("this slot is unmodifiable with [<-"))

  slot(x, i) <- value

  validObject_universalmotif(x)
  x

})

#' @param .Object [universalmotif-class] Final motif.
#' @param name `character(1)` Motif name.
#' @param altname `character(1)` Alternate motif name.
#' @param family `character(1)` Transcription factor family.
#' @param organism `character(1)` Species of origin.
#' @param motif `matrix` Each column represents a position in the motif.
#' @param alphabet `character(1)` One of `c('DNA', 'RNA', 'AA')`,
#'    or a combined string representing the letters.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#' @param icscore `numeric(1)` Total information content. Automatically generated.
#' @param nsites `numeric(1)` Number of sites the motif was constructed from.
#' @param pseudocount `numeric(1)` Correction to be applied to prevent `-Inf`
#'   from appearing in PWM matrices.
#' @param bkg `numeric` A vector of probabilities, each between 0 and 1. If
#'    higher order backgrounds are provided, then the elements of the vector
#'    must be named.
#' @param bkgsites `numeric(1)` Total number of sites used to find the motif.
#' @param consensus `character(1)` Consensus string. Automatically generated for 
#'    'DNA', 'RNA', and 'AA' alphabets.
#' @param strand `character(1)` Whether the motif is specific to a certain strand.
#' @param pval `numeric(1)` P-value associated with motif.
#' @param qval `numeric(1)` Adjusted P-value associated with motif.
#' @param eval `numeric(1)` E-value associated with motif.
#' @param multifreq `list` See [add_multifreq()].
#' @param extrainfo `character` Any other extra information, represented as
#'    a named character vector.
#' @param gapinfo `universalmotif_gapped(1)` Gapped motif information.
#' @name universalmotif
#' @rdname universalmotif-class
#' @aliases initialize,universalmotif-method
setMethod("initialize", signature = "universalmotif",
          definition = function(.Object, name, altname, family, organism, motif,
                                alphabet = "DNA", type, icscore, nsites,
                                pseudocount = 1, bkg, bkgsites, consensus,
                                strand = "+-", pval, qval, eval, multifreq,
                                extrainfo, gapinfo) {

  message("Please use create_motif() instead.")

  if (!missing(name)) .Object@name <- name
  if (!missing(altname)) .Object@altname
  if (!missing(family)) .Object@family <- family
  if (!missing(organism)) .Object@organism <- organism
  if (!missing(motif)) .Object@motif <- motif
  .Object@alphabet <- alphabet
  if (!missing(type)) .Object@type <- type
  if (!missing(icscore)) .Object@icscore <- icscore
  if (!missing(nsites)) .Object@nsites <- nsites
  .Object@pseudocount <- pseudocount
  if (!missing(bkg)) .Object@bkg <- bkg
  if (!missing(bkgsites)) .Object@bkgsites <- bkgsites
  if (!missing(consensus)) .Object@consensus <- consensus
  .Object@strand <- strand
  if (!missing(pval)) .Object@pval <- pval
  if (!missing(qval)) .Object@qval <- qval
  if (!missing(eval)) .Object@eval <- eval
  if (!missing(multifreq)) .Object@multifreq <- multifreq
  if (!missing(extrainfo)) .Object@extrainfo <- extrainfo
  if (!missing(gapinfo)) .Object@gapinfo <- gapinfo

  validObject(.Object)
  .Object

})

#' @param object [universalmotif-class] Motif.
#' @rdname universalmotif-class
#' @aliases show,universalmotif-method
setMethod("show", signature = "universalmotif",
          definition = function(object) {

  obj.check <- tryCatch(validObject_universalmotif(object, FALSE),
                        error = function(e) return("error"))
  if (length(obj.check) > 0)
    warning(wmsg("Something is wrong with the universalmotif object and it ",
                 "may display incorrectly. Run validObject(motif)",
                 " to diagnose the problem."),
            immediate. = TRUE, call. = FALSE)

  name <- object@name
  if (nchar(name) > 40) name <- collapse_cpp(c(substr(name, 1, 40), "..."))

  cat(collapse_cpp(c("\n       Motif name:   ", name, "\n")))

  if (length(object@altname) > 0) {
    altname <- object@altname
    if (nchar(altname) > 40)
      altname <- collapse_cpp(c(substr(altname, 1, 40), "..."))
    cat(collapse_cpp(c("   Alternate name:   ", altname, "\n")))
  }

  if (length(object@family) > 0) {
    family <- object@family
    if (nchar(family) > 40)
      family <- collapse_cpp(c(substr(family, 1, 40), "..."))
    cat(collapse_cpp(c("           Family:   ", family, "\n")))
  }

  if (length(object@organism) > 0) {
    organism <- object@organism
    if (nchar(organism) > 40)
      organism <- collapse_cpp(c(substr(organism, 1, 40), "..."))
    cat(collapse_cpp(c("         Organism:   ", organism, "\n")))
  }

  alphabet <- object@alphabet
  if (nchar(alphabet) > 40)
    alphabet <- collapse_cpp(c(substr(alphabet, 1, 40), "..."))

  cat(collapse_cpp(c("         Alphabet:   ", alphabet, "\n")))
  cat(collapse_cpp(c("             Type:   ", object@type, "\n")))

  if (object@alphabet %in% c("DNA", "RNA"))
    cat(collapse_cpp(c("          Strands:   ", object@strand, "\n")))

  cat(collapse_cpp(c("         Total IC:   ", round(object@icscore, 2), "\n")))

  if (length(object@consensus) > 0) {
    consensus <- object@consensus
    if (object@gapinfo@isgapped) {
      consensus <- safeExplode(consensus)
      gaploc1 <- c(1, object@gapinfo@gaploc + 1)
      gaploc2 <- c(object@gapinfo@gaploc, length(consensus))

      consensus <- mapply(function(x, y) collapse_cpp(consensus[x:y]),
                          gaploc1, gaploc2, SIMPLIFY = FALSE)
      consensus[-length(consensus)] <- lapply(consensus[-length(consensus)],
                                              function(x) paste0(x, ".."))
      consensus <- collapse_cpp(unlist(consensus))
    }
    if (nchar(consensus) > 40)
      consensus <- collapse_cpp(c(substr(consensus, 1, 40), "..."))
    cat(collapse_cpp(c("        Consensus:   ", consensus, "\n")))
  }

  if (length(object@nsites) > 0)
    cat(collapse_cpp(c("     Target sites:   ", object@nsites, "\n")))

  if (length(object@bkgsites) > 0)
    cat(collapse_cpp(c(" Background sites:   ", object@bkgsites, "\n")))

  if (length(object@pval) > 0)
    cat(collapse_cpp(c("          P-value:   ", object@pval, "\n")))

  if (length(object@qval) > 0)
    cat(collapse_cpp(c("          Q-value:   ", object@qval, "\n")))

  if (length(object@eval) > 0)
    cat(collapse_cpp(c("          E-value:   ", object@eval, "\n")))

  if (length(object@multifreq) > 0) {
    toprint <- paste0(names(object@multifreq), collapse = ", ")
    cat(collapse_cpp(c("   k-letter freqs:   ", toprint, "\n")))
  }

  if (object@gapinfo@isgapped) {
    g <- list(object@gapinfo@mingap, object@gapinfo@maxgap)
    g <- mapply(function(x, y) paste0(x, "-", y), g[[1]], g[[2]])
    g <- paste0(g, collapse = ", ")
    if (nchar(g) > 40) g <- collapse_cpp(c(substr(g, 1, 40), "..."))
    l <- object@gapinfo@gaploc
    l <- list(l, l + 1)
    l <- mapply(function(x, y) paste0(x, "-", y), l[[1]], l[[2]])
    l <- paste0(l, collapse = ", ")
    if (nchar(l) > 40) l <- collapse_cpp(c(substr(l, 1, 40), "..."))
    cat(collapse_cpp(c("    Gap locations:   ", l, "\n")))
    cat(collapse_cpp(c("        Gap sizes:   ", g, "\n")))
  }

  if (length(object@extrainfo) > 0 ) {

    extrainfo <- object@extrainfo
    cat("       Extra info:   ")

    if (length(extrainfo) > 3) extrainfo <- c(extrainfo[1:3], c("..." = "..."))

    for (i in seq_along(extrainfo)) {

      if (!is.null(names(extrainfo)[i]) && names(extrainfo)[i] != "...")
        to_show <- collapse_cpp(c("[", names(extrainfo[i]), "] ", extrainfo[i]))
      else
        to_show <- extrainfo[i]

      if (nchar(to_show) > 40)
        to_show <- collapse_cpp(c(substr(to_show, 1, 40), "..."))

      if (i == 1) {
        cat(to_show, "\n")
        next
      }

      cat(collapse_cpp(c("                     ", to_show, "\n")))

    }

  }

  cat("\n")

  motif <- round(object@motif * 100) / 100

  if (!object@gapinfo@isgapped) {

    print(motif)

  } else {

    gapmat <- matrix(NaN, nrow = nrow(motif), ncol = 1,
                     dimnames = list(rownames(motif), "000"))

    motif <- split_gapped(motif, object@gapinfo@gaploc)

    motif[-length(motif)] <- lapply(motif[-length(motif)], function(x) cbind(x, gapmat))
    motif <- do.call(cbind, motif)

    o <- capture.output(print(motif))
    o <- gsub("000", "  ", o)
    o <- gsub("NaN", "..", o)
    cat(o, sep = "\n")

  }

  invisible(object)

})

# as.matrix
# as.character
# names

#' @export
BiocGenerics::cbind

#' @export
BiocGenerics::subset

#' @export
BiocGenerics::rownames

#' @export
BiocGenerics::colnames

#' @export
BiocGenerics::ncol

#' @export
BiocGenerics::nrow

#' @export
BiocGenerics::rowSums

#' @export
BiocGenerics::colSums

#' @export
BiocGenerics::colMeans

#' @export
BiocGenerics::rowMeans

#' @export
BiocGenerics::normalize

#' @export
BiocGenerics::as.data.frame

#' @rdname universalmotif-class
#' @aliases as.data.frame,universalmotif-method
setMethod("as.data.frame", signature(x = "universalmotif"),
          definition = function(x) {

  data.frame(

    name = x@name,
    altname = ifelse(length(x@altname) == 0, NA_character_, x@altname),
    family = ifelse(length(x@family) == 0, NA_character_, x@family),
    organism = ifelse(length(x@organism) == 0, NA_character_, x@organism),
    alphabet = x@alphabet,
    icscore = x@icscore,
    nsites = ifelse(length(x@nsites) == 0, NA_real_, x@nsites),
    bkgsites = ifelse(length(x@bkgsites) == 0, NA_real_, x@bkgsites),
    consensus = ifelse(length(x@consensus) == 0, NA_character_, x@consensus),
    strand = x@strand,
    pval = ifelse(length(x@pval) == 0, NA_real_, x@pval),
    qval = ifelse(length(x@qval) == 0, NA_real_, x@qval),
    eval = ifelse(length(x@eval) == 0, NA_real_, x@eval),

    check.names = FALSE,
    fix.empty.names = FALSE,
    stringsAsFactors = FALSE

  )

})

#' @param select `numeric` Columns to keep.
#' @rdname universalmotif-class
#' @aliases subset,universalmotif-method
setMethod("subset", signature(x = "universalmotif"),
          definition = function(x, select) {

  mot <- x@motif[, select]
  motif <- universalmotif_cpp(motif = mot, name = x@name,
                              altname = x@altname, family = x@family,
                              organism = x@organism,
                              alphabet = x@alphabet, nsites = x@nsites,
                              pseudocount = x@pseudocount, bkg = x@bkg,
                              bkgsites = x@bkgsites, strand = x@strand,
                              pval = x@pval, qval = x@qval, eval = x@eval,
                              extrainfo = x@extrainfo)

  validObject_universalmotif(motif)
  motif

})

#' @rdname universalmotif-class
#' @aliases normalize,universalmotif-method
setMethod("normalize", signature(object = "universalmotif"),
          definition = function(object) {
  type <- object@type
  pseudo <- object@pseudocount
  if (pseudo == 0) pseudo <- 1
  object <- convert_type(object, "PCM")
  convert_type(object, type, pseudocount = pseudo)
})

#' @rdname universalmotif-class
#' @aliases rowMeans,universalmotif-method
setMethod("rowMeans", signature(x = "universalmotif"),
          definition = function(x) rowMeans(x@motif))

#' @rdname universalmotif-class
#' @aliases colMeans,universalmotif-method
setMethod("colMeans", signature(x = "universalmotif"),
          definition = function(x) colMeans(x@motif))

#' @rdname universalmotif-class
#' @aliases colSums,universalmotif-method
setMethod("colSums", signature(x = "universalmotif"),
          definition = function(x) colSums(x@motif))

#' @rdname universalmotif-class
#' @aliases rowSums,universalmotif-method
setMethod("rowSums", signature(x = "universalmotif"),
          definition = function(x) rowSums(x@motif))

#' @rdname universalmotif-class
#' @aliases nrow,universalmotif-method
setMethod("nrow", signature = "universalmotif",
          definition = function(x) nrow(x@motif))

#' @rdname universalmotif-class
#' @aliases ncol,universalmotif-method
setMethod("ncol", signature = "universalmotif",
          definition = function(x) ncol(x@motif))

#' @rdname universalmotif-class
#' @aliases colnames,universalmotif-method
setMethod("colnames", signature(x = "universalmotif"),
          definition = function(x) colnames(x@motif))

#' @rdname universalmotif-class
#' @aliases rownames,universalmotif-method
setMethod("rownames", signature(x = "universalmotif"),
          definition = function(x) rownames(x@motif))

#' @export
c.universalmotif <- function(...) list(...)

#' @param ... [universalmotif-class] Motifs.
#' @param deparse.level Unused.
#' @rdname universalmotif-class
#' @aliases cbind,universalmotif-method
setMethod("cbind", signature = "universalmotif",
          definition = function(..., deparse.level = 0) {

  mots <- list(...)
  for (i in seq_along(mots)) validObject_universalmotif(mots[[i]])

  if (length(mots) == 1) return(mots[[1]])

  mot.alphabet <- unique(vapply(mots, function(x) x@alphabet, character(1)))
  if (length(mot.alphabet) != 1) stop("all motifs must have the same alphabet")
  mots <- convert_type(mots, "PPM")
  mot.name <- unique(vapply(mots, function(x) x@name, character(1)))
  mot.altname <- unique(do.call(c, lapply(mots, function(x) x@altname)))
  mot.family <- unique(do.call(c, lapply(mots, function(x) x@family)))
  mot.organism <- unique(do.call(c, lapply(mots, function(x) x@organism)))
  mot.motif <- lapply(mots, function(x) x@motif)
  mot.nsites <- do.call(c, lapply(mots, function(x) x@nsites))
  mot.bkg <- lapply(mots, function(x) x@bkg)
  mot.bkgsites <- lapply(mots, function(x) x@bkgsites)
  mot.bkgsites <- do.call(c, mot.bkgsites)
  mot.strand <- unique(vapply(mots, function(x) x@strand, character(1)))
  mot.pval <- do.call(c, lapply(mots, function(x) x@pval))
  mot.qval <- do.call(c, lapply(mots, function(x) x@qval))
  mot.eval <- do.call(c, lapply(mots, function(x) x@eval))
  mot.extrainfo <- lapply(mots, function(x) x@extrainfo)
  mot.extrainfo <- unique(do.call(c, mot.extrainfo))

  mot.motif <- do.call(cbind, mot.motif)
  mot.name <- paste0(mot.name, collapse = "/")

  if (length(mot.altname) > 0) mot.altname <- paste0(mot.altname, collapse = "/")
  if (length(mot.family) > 0) mot.family <- paste0(mot.family, collapse = "/")
  if (length(mot.organism) > 0) mot.organism <- paste0(mot.organism, collapse = "/")
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

  validObject_universalmotif(motif)
  motif

})
