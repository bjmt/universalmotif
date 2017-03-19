######################################################################
## Benjamin Tremblay
##
## Misc methods for universalmotif-class
##
######################################################################

setMethod("show", signature = "universalmotif",
          definition = function(object) {
            cat(" Motif name:   ", object@name, "\n",
                "       Type:   ", object@type, "\n", 
                "    Strands:   ", paste(object@strand, collapse = " "), "\n",
                # "         IC:   ", paste(object@icscores, collapse = ", "), "\n",
                "   Total IC:   ", object@icscore, "\n", 
                "  Consensus:   ", object@consensus, "\n", sep = "")
            if (length(object@extra) > 0) {
              cat(" Extra info:   ", paste(names(object@extra), collapse = " "),
                  "\n", sep = "")
            }
            cat("\n")
            print(object@motif)
            invisible(NULL)
          })

setMethod("initialize", signature = "universalmotif",
          definition = function(.Object, name, motif,
                      alphabet = "DNA", #letters = character(0),
                      type = character(0), icscore = numeric(0),
                      nsites = numeric(0), pseudoweight = 0.8,
                      bkg = numeric(0), consensus = character(0),
                      strand = c("+", "-"), extra = character(0)) {

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

            .Object@extra <- extra
            
            .Object
          })
