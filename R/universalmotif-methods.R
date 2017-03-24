######################################################################
## Benjamin Tremblay
##
## Misc methods for universalmotif-class
##
######################################################################

#' @describeIn universalmotif Show method.
#' @include utils.R universalmotif-class.R AllGenerics.R
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
                # - motif matrix

            if (missing(name)) name <- "placeholder-name"
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

            if (missing(type) || !type %in% c("PCM", "PPM", "PWM", "ICM")) {
              if (all(colSums(motif) >= 2)) type <- "PCM" else {
                if (all(motif >= 0)) type <- "PPM" else {
                  type <- "PWM"
                  message("assumed 'type' as being PWM")
                }
              }
            }
            .Object@type <- type

            if (type == "PCM" && length(nsites) == 0) {
              nsites <- unique(colSums(motif))
              if (length(nsites) != 1) stop("all colSums must be equal")
            }
            .Object@nsites <- nsites

            .Object@pseudoweight <- pseudoweight

            if (missing(bkg) &&
                alphabet %in% c("DNA", "RNA")) bkg <- c(0.25, 0.25, 0.25, 0.25) 
            .Object@bkg <- bkg

            consensus <- apply(motif, 2, get_consensus, alphabet = alphabet,
                               type = type, pseudoweight = pseudoweight)

            consensus <- paste(consensus, collapse = "")
            .Object@consensus <- consensus
            
            icscores <- apply(motif, 2, position_icscore,
                                 bkg = bkg, type = type,
                                 pseudoweight = pseudoweight, nsites = nsites) 
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
setMethod("motif_slots", "universalmotif", function(object, slots) {

          if (missing(slots)) {
            slots <- c("name", "motif", "alphabet", "type", "icscore",
                       "nsites", "pseudoweight", "bkg", "consensus",
                       "strand", "pval", "extrachar", "extranum",
                       "eval", "bkgsites")
          }
          
          if (all(slots == "motif")) return(object@motif)

          returnlist <- lapply(slots, function(x) slot(object, x))
          names(returnlist) <- slots
          if ("motif" %in% names(returnlist)) {
            returnlist$motif <- object@motif
          }

          returnlist <- returnlist[vapply(returnlist, function(x) length(x) > 0,
                                          logical(1))]
          if (length(returnlist) == 1) returnlist <- unlist(returnlist)

          return(returnlist)

         })

#' @describeIn universalmotif Replacement method.
setMethod("motif_slots<-", "universalmotif", function(object, slot, value) {
           slot(object, slot) <- value
           if (validObject(object)) return(object)
         })

#' @describeIn convert_type Convert type for \linkS4class{universalmotif}.
setMethod("convert_type", "universalmotif", function(motif, out_type,
                                                     pseudoweight = NULL,
                                                     background = NULL,
                                                     IC_floor = FALSE) {

            # the functions referred to within are stored in utils.R
            # (purely for historical reasons; may be moved at some point)

            # NOTE: the universalmotif ICM implementation is different from
            # the TFBSTools::toICM implementation!
            # (as far I can tell, the universalmotif PWM and TFBSTools::toPWM
            # implementations are identical)

            # CAREFUL WITH IFELSE!! It will recycle 'yes' and 'no' arguments!

            if (motif_slots(motif, "type") == out_type) return(motif)

            if (motif@type == "PCM") {
              if (out_type == "PPM") {
                possums <- colSums(motif@motif)
                for (i in seq_len(ncol(motif@motif))) {
                  motif@motif[, i] <- pcm_to_ppm(motif@motif[, i],
                                       possum = possums[i],
                                       pseudoweight = ifelse(is.null(pseudoweight),
                                                             motif@pseudoweight,
                                                             pseudoweight))
                }
                motif@type <- "PPM"
                return(motif)
              }
              if (out_type == "PWM") {
                motif@motif <- apply(motif@motif, 2, pcm_to_ppm,
                                     pseudoweight = ifelse(is.null(pseudoweight),
                                                           motif@pseudoweight,
                                                           pseudoweight))
                motif@motif <- apply(motif@motif, 2, ppm_to_pwm,
                                     background = if (!is.null(background)) 
                                                         background else motif@bkg,
                                     pseudoweight = ifelse(is.null(pseudoweight),
                                                           motif@pseudoweight,
                                                           pseudoweight),
                                     nsites = motif@nsites)
                motif@type <- "PWM"
                return(motif)
              }
              if (out_type == "ICM") {
                motif@motif <- apply(motif@motif, 2, pcm_to_ppm,
                                     pseudoweight = ifelse(is.null(pseudoweight),
                                                           motif@pseudoweight,
                                                           pseudoweight))
                motif@motif <- apply(motif@motif, 2, ppm_to_icm,
                                     bkg = background, IC_floor = IC_floor)
                motif@type <- "ICM"
                return(motif)
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
                                     background = if (!is.null(background))
                                                         background else motif@bkg,
                                     pseudoweight = ifelse(is.null(pseudoweight),
                                                           motif@pseudoweight,
                                                           pseudoweight),
                                     nsites = motif@nsites)

                motif@type <- "PWM"
                return(motif)
              }
              if (out_type == "ICM") {
                motif@motif <- apply(motif@motif, 2, ppm_to_pcm,
                                     nsites = motif@nsites)
                motif@motif <- apply(motif@motif, 2, pcm_to_ppm,
                                     pseudoweight = ifelse(is.null(pseudoweight),
                                                           motif@pseudoweight,
                                                           pseudoweight))
                motif@motif <- apply(motif@motif, 2, ppm_to_icm,
                                     bkg = background, IC_floor = IC_floor)
                motif@type <- "ICM"
                return(motif)
              }
            }

            if (motif@type == "PWM") {
              if (out_type == "PCM") {
                motif@motif <- apply(motif@motif, 2, pwm_to_ppm,
                                     background = if (!is.null(background))
                                                         background else motif@bkg)
                motif@motif <- apply(motif@motif, 2, ppm_to_pcm,
                                     nsites = motif@nsites)
                motif@type <- "PCM"
                return(motif)
              }
              if (out_type == "PPM") {
                motif@motif <- apply(motif@motif, 2, pwm_to_ppm,
                                     background = if (!is.null(background))
                                                         background else motif@bkg)
                motif@type <- "PPM"
                return(motif)
              }
              if (out_type == "ICM") {
                motif@motif <- apply(motif@motif, 2, pwm_to_ppm,
                                     background = if (!is.null(background))
                                                         background else motif@bkg)
                motif@motif <- apply(motif@motif, 2, ppm_to_icm,
                                     bkg = background, IC_floor = IC_floor)
                motif@type <- "ICM"
                return(motif)
              }
            }

            if (motif@type == "ICM") {
              stop("motifs cannot be converted from ICM", call. = FALSE)
            }

            stop("unrecognized motif type", call. = FALSE)

         })

#' @describeIn filter_motifs Filter \linkS4class{universalmotif}.
setMethod("filter_motifs", signature = "universalmotif",
          definition = function(motifs, min_width, alphabet, type, icscore,
                                nsites, strand, pval, eval, bkgsites,
                                extrachar, extranum) {

            if (!missing(min_width)) {
              if (ncol(motifs@motif) > min_width) return(NULL)
            }

            if (!missing(alphabet)) {
              if (motifs@alphabet != alphabet) return(NULL)
            }

            if (!missing(type)) {
              if (motifs@type != type) return(NULL)
            }

            if (!missing(icscore)) {
              if (motifs@icscore > icscore) return(NULL)
            }

            if (!missing(nsites)) {
              if (motifs@nsites > nsites) return(NULL)
            }

            if (!missing(strand)) {
              if (motifs@strand != strand) return(NULL)
            }

            if (!missing(pval)) {
              if (motifs@pval < pval) return(NULL)
            }

            if (!missing(eval)) {
              if (motifs@eval < eval) return(NULL)
            }
            
            if (!missing(bkgsites)) {
              if (motifs@bkgsites > bkgsites) return(NULL)
            }
            
            if (!missing(extrachar)) {
              charnames <- names(extrachar)
              if (!any(mapply(function(x, y, z) {
                                checkdel <- x == z[which(names(z) == y)]
                                if (length(checkdel) == 0) return(FALSE)
                                return(checkdel)
                              },
                              x = extrachar, y = charnames,
                              MoreArgs = list(z = motifs@extrachar)))) {
                 return(NULL)
              }
            }
                             

            # right now the extranum filter is pretty useless, but I don't
            # think there's a quick solution

            if (!missing(extranum)) {
              numnames <- names(numchar)
              if (!any(mapply(function(x, y, z) {
                                checkdel <- x == z[which(names(z) == y)]
                                if (length(checkdel) == 0) FALSE else checkdel
                              },
                              x = extranum, y = numnames,
                              MoreArgs = list(z = motifs@extranum)))) {
                 return(NULL)
              }
            }

            return(motifs)

          })

#' @describeIn trim_motifs Trim \linkS4class{universalmotif}.
setMethod("trim_motifs", signature = "universalmotif",
          definition = function(motifs, ic_cutoff = 0.25) {
            thescores <- apply(motifs@motif, 2, position_icscore,
                               bkg = motifs@bkg, type = motifs@type,
                               pseudoweight = motifs@pseudoweight,
                               nsites = motifs@nsites)
            tocut <- vector(length = length(thescores))
            for (i in seq_along(thescores)) {
              if (thescores[i] <= ic_cutoff) tocut[i] <- i else break
            }
            for (i in rev(seq_along(thescores))) {
              if (thescores[i] <= ic_cutoff) tocut[i] <- i else break
            }
            tocut <- tocut[tocut != 0]
            motifs@motif <- motifs@motif[, -tocut]
            motifs@consensus <- apply(motifs@motif, 2, get_consensus,
                                      alphabet = motifs@alphabet,
                                      type = motifs@type,
                                      pseudoweight = motifs@pseudoweight)
            motifs@icscore <- apply(motifs@motif, 2, position_icscore,
                                    bkg = motifs@bkg, type = motifs@type,
                                    pseudoweight = motifs@pseudoweight,
                                    nsites = motifs@nsites)
            motifs@icscore <- sum(motifs@icscore)
            if (ncol(motifs@motif) == 0) return(NULL) else return(motifs)
          })
