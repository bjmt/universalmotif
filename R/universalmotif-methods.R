######################################################################
## Benjamin Tremblay
##
## Misc methods for universalmotif-class
##
######################################################################

#' @describeIn universalmotif Show method.
#' @include utils.R universalmotif-class.R generics.R
setMethod("show", signature = "universalmotif",
          definition = function(object) {
            cat("\n       Motif name:   ", object@name, "\n",
                "             Type:   ", object@type, "\n", 
                "          Strands:   ", paste(object@strand, collapse = " "), "\n",
                "         Total IC:   ", object@icscore, "\n", 
                "        Consensus:   ", object@consensus, "\n", sep = "")
            if (length(object@nsites) > 0) {
              cat("     Target sites:   ", object@nsites, "\n", sep = "")
            }
            if (length(object@bkgsites) > 0) {
              cat(" Background sites:   ", object@bkgsites, "\n", sep = "")
            }
            if (length(object@pval) > 0) {
              cat("          P-value:   ", object@pval, "\n", sep = "")
            }
            if (length(object@eval) > 0) {
              cat("          E-value:   ", object@eval, "\n", sep = "")
            }
            if (length(object@extrachar) > 0 || length(object@extranum) > 0) {
              cat("       Extra info:   ")
              if (length(object@extrachar) > 0) {
                cat(paste(names(object@extrachar), collapse = ", "))
                if (length(object@extranum) > 0) cat(", ")
              }
              if (length(object@extranum) > 0) {
                cat(paste(names(object@extranum), collapse = ", "))
              }
              cat("\n")
            }
            cat("\n")
            print(object@motif)
            invisible(NULL)
          })

#' @describeIn universalmotif Initialize method.
setMethod("initialize", signature = "universalmotif",
          definition = function(.Object, name, motif,
                      alphabet = "DNA", #letters = character(0),
                      type = character(0), icscore = numeric(0),
                      nsites = numeric(0), pseudoweight = 0.8,
                      bkg = numeric(0), consensus = character(0),
                      strand = c("+", "-"), extrachar = character(0),
                      extranum = numeric(0), pval = numeric(0),
                      eval = numeric(0), bkgsites = numeric(0)) {

            # required fields for construction:
                # - name
                # - motif matrix
                # - type

            if (missing(name)) stop("motif must have a name")
            .Object@name <- name

            if (missing(motif)) stop("missing motif matrix")

            if (!missing(alphabet)) {
              if (!alphabet %in% c("DNA", "RNA", "AA", "custom")) {
                stop("acceptable alphabets are: 'DNA', 'RNA', 'AA', or 'custom'")
              }
            }
            .Object@alphabet <- alphabet

            if (alphabet == "DNA") letters <- c("A", "C", "G", "T")
            if (alphabet == "RNA") letters <- c("A", "C", "G", "U")
            if (alphabet == "AA") letters <-  c("A", "C", "D", "E", "F",
                                                "G", "H", "I", "K", "L",
                                                "M", "N", "P", "Q", "R",
                                                "S", "T", "V", "W", "Y")
            # .Object@letters <- letters

            if (nrow(motif) != length(letters)) {
              stop("matrix must have same number of rows as letters")
            } else {
              if (!all(letters %in% rownames(motif))) rownames(motif) <- letters
            }
            colnames(motif) <- NULL
            .Object@motif <- motif

            if (missing(type) || !type %in% c("PCM", "PPM", "PWM")) {
              stop("type must be provided as 'PCM', 'PPM', or 'PWM'")
            }
            .Object@type <- type

            .Object@nsites <- nsites

            .Object@pseudoweight <- pseudoweight

            if (missing(bkg) &&
                alphabet %in% c("DNA", "RNA")) bkg <- c(0.25, 0.25, 0.25, 0.25) 
            .Object@bkg <- bkg

            if (alphabet == "DNA") consensus <- apply(motif, 2, get_consensus,
                                                      alphabet = "DNA",
                                                      type = type,
                                                      pseudoweight = pseudoweight)
            if (alphabet == "RNA") consensus <- apply(motif, 2, get_consensus,
                                                      alphabet = "RNA",
                                                      type = type,
                                                      pseudoweight = pseudoweight)
            .Object@consensus <- consensus
            
            icscores <- apply(motif, 2, position_icscore,
                                 bkg = bkg, type = type,
                                 pseudoweight = pseudoweight) 
            # .Object@icscores <- icscores
            .Object@icscore <- sum(icscores)


            if (!all(strand == c("+", "-") && !any(strand %in% c("+", "-")))) {
              strand <- c("+", "-")
            }
            .Object@strand <- strand

            .Object@extrachar <- extrachar

            .Object@extranum <- extranum

            .Object@pval <- pval

            .Object@bkgsites <- bkgsites

            .Object@eval <- eval
            
            .Object

          })

#' @describeIn universalmotif Accessor function.
setMethod("motif_slots", "universalmotif", function(object) {
          toreturn1 <- NULL
          toreturn2 <- NULL
          if (length(object@extranum) > 0) toreturn1 <- object@extranum
          if (length(object@extrachar) > 0) toreturn2 <- object@extrachar
          if (is.null(toreturn1) && is.null(toreturn2)) return(invisible(NULL))
          return(list(extranum = toreturn1, extrachar = toreturn2))
         })

#' @describeIn universalmotif Convert type between PCM, PPM and PWM.
setMethod("convert_type", "universalmotif", function(motif, out_type) {

            if (motif@type == "PCM") {
              if (out_type == "PPM") {
                possums <- colSums(motif@motif)
                for (i in seq_len(ncol(motif@motif))) {
                  motif@motif[, i] <- pcm_to_ppm(motif@motif[, i],
                                       possum = possums[i],
                                       pseudoweight = motif@pseudoweight)
                }
                motif@type <- "PPM"
                return(motif)
              }
              if (out_type == "PWM") {
                motif@motif <- apply(motif@motif, 2, pcm_to_ppm)
                motif@motif <- apply(motif@motif, 2, ppm_to_pwm,
                                     background = motif@bkg,
                                     pseudoweight = motif@pseudoweight,
                                     nsites = motif@nsites)
                motif@type <- "PWM"
              }
            }

            if (motif@type == "PPM") {
              if (out_type == "PCM") {
                motif@motif <- apply(motif@motif, 2, ppm_to_pcm,
                                     nsites = ifelse(length(motif@nsites) > 0,
                                                     motif@nsites, 100))
                motif@type <- "PCM"
                return(motif)
              }
              if (out_type == "PWM") {
                motif@motif <- apply(motif@motif, 2, ppm_to_pwm,
                                     background = motif@bkg,
                                     pseudoweight = motif@pseudoweight,
                                     nsites = motif@nsites)

                motif@type <- "PWM"
              }
            }

            if (motif@type == "PWM") {
              if (out_type == "PCM") {
                motif@motif <- apply(motif@motif, 2, pwm_to_ppm)
                motif@motif <- apply(motif@motif, 2, ppm_to_pcm)
                motif@type <- "PCM"
              }
              if (out_type == "PPM") {
                motif@motif <- apply(motif@motif, 2, pwm_to_ppm)
                motif@type <- "PPM"
              }
            }

            return(motif)

         })
