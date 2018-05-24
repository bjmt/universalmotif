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

#' @describeIn create_motif Create motif from a consensus string.
setMethod("create_motif", signature(consensus = "character",
                                    matrix = "missing",
                                    sequences = "missing"),
          definition = function(consensus, name, out_type,
                                out_class, pseudoweight, alphabet,
                                background, nsites) {
            consensus <- strsplit(consensus, split = "")[[1]]
            motif <- vapply(consensus, consensus_to_ppm, numeric(4))
            motif <- universalmotif(name = name, motif = motif,
                                    pseudoweight = pseudoweight,
                                    alphabet = alphabet,
                                    bkg = background,
                                    nsites = nsites)
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
            if (length(unique(Biostrings::width(sequences))) != 1) {
              stop("all sequences must be of the same width", call. = FALSE)
            }
            sequences <- as.character(sequences)
            sequences <- lapply(sequences, function(x)
                                strsplit(x, split = "")[[1]])
            sequences <- withinlistvapply(sequences, function(x)
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

#' @describeIn filter_motifs Filter a list of motifs.
setMethod("filter_motifs", signature = "list",
          definition = function(motifs, ...) {
            motifs <- lapply(motifs, filter_motifs, ...)
            motifs <- motifs[vapply(motifs, function(x) !is.null(x), logical(1))]
            return(motifs)
          })

#' @describeIn trim_motifs Trim a list of motifs.
setMethod("trim_motifs", signature = "list",
          definition = function(motifs, ic_cutoff) {
            lapply(motifs, trim_motifs, ic_cutoff = ic_cutoff)
          })

#' @describeIn trim_motifs Trim a list of motifs which are not
#'    \linkS4class{universalmotif}.
setMethod("trim_motifs", signature = "ANY",
          definition = function(motifs, ic_cutoff) {
            theclass <- class(motifs)
            theclasses <- c("PFMatrix" = "TFBSTools-PFMatrix",
                            "PWMatrix" = "TFBSTools-PWMatrix",
                            "ICMatrix" = "TFBSTools-ICMatrix",
                            "pwm" = "seqLogo-pwm",
                            "pcm" = "motifStack-pcm",
                            "pfm" = "motifStack-pfm",
                            "PWM" = "PWMEnrich-PWM",
                            "motif" = "rGADEM-motif")
            motifs <- convert_motifs(motifs)
            motifs <- trim_motifs(motifs = motifs, ic_cutoff = ic_cutoff)
            motifs <- convert_motifs(motifs, 
                                     theclasses[names(theclasses) == theclass]
                                     )
            return(motifs)
          })
