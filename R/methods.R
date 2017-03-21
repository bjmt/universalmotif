######################################################################
## Benjamin Tremblay
##
## assorted methods with no home
##
######################################################################

#' @include universalmotif-methods.R
#' @describeIn convert_type Convert type on a list of motifs
setMethod("convert_type", "list", function(motif, out_type,
                                                     pseudoweight = NULL) {
  lapply(motif, convert_type, out_type, pseudoweight)
})

#' @describeIn convert_motifs Convert a list of motifs.
setMethod("convert_motifs", signature = "list",
          definition = function(motif, out_class, ...) {
            lapply(X = motif, FUN = convert_motifs, out_class, ...)
          })

#' @describeIn create_motif Create motif from a consensus string.
setMethod("create_motif", signature(consensus = "character",
                                    matrix = "missing",
                                    sequences = "missing"),
          definition = function(consensus, name, out_type,
                                out_class, pseudoweight, alphabet,
                                background, nsites) {
            consensus <- strsplit(consensus, split = "")[[1]]
            motif <- vapply(consensus, consensus_to_ppm, numeric(4))
            if (missing(name)) stop("please provide a name", call. = FALSE)
            motif <- universalmotif(name = name, motif = motif,
                                    pseudoweight = pseudoweight,
                                    alphabet =  alphabet,
                                    bkg = background,
                                    nsites = nsites,
                                    type = out_type)
            if (out_type != "PPM") {
              motif <- convert_type(motif, out_type = out_type,
                                    pseudoweight = pseudoweight,
                                    background = background)
            }
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class,
                                    pseudoweight = pseudoweight,
                                    background = background)
            return(motif)
          })

#' @describeIn create_motif Create motif from a matrix object.
setMethod("create_motif", signature(consensus = "missing",
                                    matrix = "matrix",
                                    sequences = "missing"),
          definition = function(matrix, name, out_type,
                                out_class, pseudoweight, alphabet,
                                background, nsites) {
            if (missing(name)) stop("please provide a name", call. = FALSE)
            motif <- universalmotif(name = name, motif = matrix,
                                    pseudoweight = pseudoweight,
                                    alphabet = alphabet,
                                    bkg = background,
                                    nsites = nsites)
            if (out_type != motif_slots(motif, "type")) {
              motif <- convert_type(motif, out_type = out_type,
                                    pseudoweight = pseudoweight,
                                    background = background)
            }
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class,
                                    pseudoweight = pseudoweight,
                                    background = background)
            return(motif)
          })

#' @describeIn create_motif Create motif from an XStringSet.
setMethod("create_motif", signature(consensus = "missing",
                                    matrix = "missing",
                                    sequences = "XStringSet"),
          definition = function(sequences, name, out_type,
                                out_class, pseudoweight, alphabet,
                                background, nsites) {
            if (missing(name)) stop("please provide a name", call. = FALSE)
            if (length(unique(Biostrings::width(sequences))) != 1) {
              stop("all sequences must be of the same width", call. = FALSE)
            }
            sequences <- as.character(sequences)
            sequences <- lapply(sequences, function(x) strsplit(x, split = "")[[1]])
            sequences <- reversevapply(sequences, function(x)
                                       paste(x, collapse = ""), character(1))
            pfm <- sapply(sequences, string_to_pfm, alphabet = alphabet)
            motif <- universalmotif(name = name, motif = pfm,
                                    pseudoweight = pseudoweight,
                                    alphabet = alphabet,
                                    bkg = background,
                                    nsites = nsites)
            if (out_type != "PCM") {
              motif <- convert_type(motif, out_type = out_type,
                                    pseudoweight = pseudoweight,
                                    background = background)
            }
            if (out_class == "universalmotif") return(motif)
            motif <- convert_motifs(motif, out_class = out_class,
                                    pseudoweight = pseudoweight,
                                    background = background)
            return(motif)
          })
