#' @describeIn create_motif Create motif from a consensus string.
#' @include universalmotif-class.R
setMethod("create_motif", signature(consensus = "character",
                                    matrix = "missing",
                                    sequences = "missing"),
          definition = function(consensus, name, pseudoweight, alphabet,
                                bkg, nsites) {
            consensus <- strsplit(consensus, split = "")[[1]]
            if (missing(alphabet)) {
              if (any(consensus %in% c("E", "F", "I", "P", "Q", "X", "Z"))) {
                motif <- vapply(consensus, consensus_to_ppmAA, numeric(20))
                alphabet <- "AA"
              } else {
                motif <- vapply(consensus, consensus_to_ppm, numeric(4))
                if (any(consensus == "U")) {
                  alphabet <- "RNA" 
                } else if (any(consensus %in% c("A", "C", "G", "T"))){
                  alphabet <- "DNA"
                }
              }
            }
            motif <- universalmotif(name = name, motif = motif,
                                    pseudoweight = pseudoweight,
                                    alphabet = alphabet, bkg = bkg,
                                    nsites = nsites)
            motif
          })

#' @describeIn create_motif Create motif from a matrix.
setMethod("create_motif", signature(consensus = "missing",
                                    matrix = "matrix",
                                    sequences = "missing"),
          definition = function(matrix, name, pseudoweight,
                                alphabet, bkg, nsites) {
            if (missing(alphabet)) {
              if (all(rownames(matrix) %in% c("A", "C", "D", "E", "F", "G", "H",
                                              "I", "K", "L", "M", "N", "P", "Q",
                                              "R", "S", "T", "V", "W", "Y"))) {
                alphabet <- "AA"
              } else if (all(rownames(matrix) %in% c("A", "C", "G", "T"))) {
                alphabet  <- "DNA"
              } else if (all(rownames(matrix) %in% c("A", "C", "G", "U"))) {
                alphabet <- "RNA"
              } else if (nrow(matrix) == 20) {
                alphabet <- "AA" 
              } else if (nrow(matrix == 4)) {
                alphabet <- "DNA" 
              } else alphabet <- "custom"
            }
            motif <- universalmotif(name = name, motif = matrix,
                                    pseudoweight = pseudoweight,
                                    alphabet = alphabet, bkg = bkg,
                                    nsites = nsites)
            motif
          })

#' @describeIn create_motif Create motif from an XStringSet object.
setMethod("create_motif", signature(consensus = "missing",
                                    matrix = "missing",
                                    sequences = "XStringSet"),
          definition = function(sequences, name, pseudoweight, alphabet,
                                bkg, nsites) {

            if (length(unique(width(sequences))) != 1) {
              stop("all sequences must be the same width")
            }
            alphabet <- sequences@elementType
            if (alphabet == "DNAString") {
              sequences <- DNAMultipleAlignment(sequences)
              sequences <- consensusMatrix(sequences)
              motif <- universalmotif(name = name,
                                      motif = sequences[1:4, ],
                                      pseudoweight = pseudoweight,
                                      alphabet = "DNA",
                                      bkg = bkg, nsites = nsites)
            } else if (alphabet == "RNAString") {
              sequences <- RNAMultipleAlignment(sequences)
              sequences <- consensusMatrix(sequences)
              motif <- universalmotif(name = name,
                                      motif = sequences[1:4, ],
                                      pseudoweight = pseudoweight,
                                      alphabet = "RNA",
                                      bkg = bkg, nsites = nsites)
            } else if (alphabet == "AAString") {
              stop("cannot handle AAString")
              # sequences <- AAMultipleAlignment(sequences)
              # sequences <- consensusMatrix(sequences)
              # motif <- universalmotif(name = name,
                                      # motif = sequences[1:20, ],
                                      # pseudoweight = pseudoweight,
                                      # alphabet = "AA",
                                      # bkg = bkg, nsites = nsites)
            }

            motif

          })

#' @describeIn convert_motifs Convert a list of motifs.
setMethod("convert_motifs", signature(motifs = "list"),
          definition = function(motifs, class) {
            lapply(motifs, function(x) convert_motifs(x, class = class))
          })

#' @describeIn convert_motifs Convert a universalmotif object.
setMethod("convert_motifs", signature(motifs = "universalmotif"),
          definition = function(motifs, class) {
            
            out_class <- strsplit(class, "-")[[1]][2]
            out_class_pkg <- strsplit(class, "-")[[1]][1]

            if (out_class == "universalmotif") {
              return(motifs)
            }

            # MotIV-pwm2
            if (out_class_pkg == "MotIV" && out_class == "pwm2") {
              motifs <- convert_type(motifs, "PPM")
              motifs <- makePWM(motifs["motif"],
                                alphabet = motifs["alphabet"])
              return(motifs)
            }

            # TFBSTools- PFMatrix, PWMatrix, and ICMatrix
            if (out_class_pkg == "TFBSTools") {
              motifs <- convert_type(motifs, "PCM")
              bkg <- motifs["bkg"]
              names(bkg) <- c("A", "C", "G", "T")
              if (length(motifs["altname"]) == 0) {
                motifs["altname"] <- ""
              }
              strand <- motifs["strand"]
              if (strand %in% c("+-", "-+")) {
                strand <- "*"
              }
              motifs <- PFMatrix(name = motifs["name"],
                                 ID = motifs["altname"],
                                 strand = strand,
                                 bg = bkg,
                                 profileMatrix = motifs["motif"])
              if (out_class == "PFMatrix") {
                motifs <- motifs
              } else if (out_class == "PWMatrix") {
                motifs <- toPWM(motifs, type = "log2probratio",
                                pseudocounts = 0.8,
                                bg = bkg)
              } else if (out_class == "ICMatrix") {
                motifs <- toICM(motifs, pseudocounts = 0.8,
                                bg = bkg)
              }
              return(motifs)
            }

            # seqLogo-pwm
            if (out_class_pkg == "seqLogo" && out_class == "pwm") {
              # if (!requireNamespace("seqLogo", quietly = TRUE)) {
                # stop("'seqLogo' package not installed")
              # }
              motifs <- convert_type(motifs, "PPM")
              motifs <- seqLogo::makePWM(motifs["motif"])
              return(motifs)
            }

            # motifStack-pcm
            if (out_class_pkg == "motifStack" && out_class == "pcm") {
              if (!requireNamespace("motifStack", quietly = TRUE)) {
                stop("'motifStack' package not installed")
              }
              pcm_class <- getClass("pcm", where = "motifStack")
              motifs <- convert_type(motifs, "PCM")
              motifs <- new(pcm_class, mat = motifs["motif"],
                            name = motifs["name"],
                            alphabet = motifs["alphabet"],
                            background = motifs["bkg"])
              return(motifs)
            }

            # motifStack-pfm
            if (out_class_pkg == "motifStack" && out_class == "pfm") {
              if (!requireNamespace("motifStack", quietly = TRUE)) {
                stop("'motifStack' package not installed")
              }
              pfm_class <- getClass("pfm", where = "motifStack")
              motifs <- convert_type(motifs, "PPM")
              motifs <- new(pfm_class, mat = motifs["motif"],
                            name = motifs["name"],
                            alphabet = motifs["alphabet"],
                            background = motifs["bkg"])
              return(motifs)
            }

            # PWMEnrich-PWM
            if (out_class_pkg == "PWMEnrich" && out_class == "PWM") {
              if (!requireNamespace("PWMEnrich", quietly = TRUE)) {
                stop("'PWMEnrich' package not installed")
              }
              motifs <- convert_type(motifs, "PCM")
              PWM_class <- getClass("PWM", where = "PWMEnrich")
              bio_mat <- matrix(as.integer(motifs["motif"]), byrow = FALSE,
                                nrow = 4)
              rownames(bio_mat) <- c("A", "C", "G", "T")
              bio_priors <- motifs["bkg"]
              names(bio_priors) <- c("A", "C", "G", "T")
              bio_mat <- PWMEnrich::PFMtoPWM(bio_mat, type = "log2probratio",
                                             prior.params = bio_priors,
                                             pseudo.count = motifs["pseudoweight"])
              motifs <- new(PWM_class, name = motifs["name"],
                            pfm = motifs["motif"],
                            prior.params = bio_priors,
                            pwm = bio_mat$pwm)
              return(motifs)
            }

            # Biostrings-PWM
            if (out_class_pkg == "Biostrings" && out_class == "PWM") {
              # if (!requireNamespace("Biostrings", quietly = TRUE)) {
                # stop("'Biostrings' package not installed")
              # }
              motifs <- convert_type(motifs, "PCM")
              bio_mat <- matrix(as.integer(motifs["motif"]), byrow = FALSE,
                                nrow = 4)
              rownames(bio_mat) <- c("A", "C", "G", "T")
              bio_priors <- motifs["bkg"]
              names(bio_priors) <- c("A", "C", "G", "T")
              motifs <- PWM(x = bio_mat, type = "log2probratio",
                            prior.params = bio_priors)
              return(motifs)
            }

            # rGADEM-motif
            if (out_class_pkg == "rGADEM" && out_class == "motif") {
              if (!requireNamespace("rGADEM", quietly = FALSE)) {
                stop("'rGADEM' package not installed")
              }
              motifs <- convert_type(motifs, "PPM")
              rGADEM_motif_class <- getClass("motif", where = "rGADEM")
              motifs <- new("motif", pwm = motifs["motif"],
                            name = motifs["name"],
                            consensus = motifs["consensus"])
              return(motifs)
            }

            stop("unknown 'class'")
          
          })

#' @describeIn convert_motifs Convert non-universalmotif class motifs.
setMethod("convert_motifs", signature(motifs = "ANY"),
          definition = function(motifs, class) {
          
            success <- FALSE

            ## convert to universalmotif
            in_class <- class(motifs)[1]
            in_class_pkg <- attributes(class(motifs))$package

            if (paste(in_class_pkg, in_class, sep = "-") == class) {
              return(motifs)
            }

            # MotifDb-MotifList
            if (in_class_pkg == "MotifDb" && in_class == "MotifList") {
              motifs_out <- list()
              motifdb_fun <- function(x) {
                universalmotif(name = x@elementMetadata@listData$providerName,
                               altname = x@elementMetadata@listData$geneSymbol,
                               family = x@elementMetadata@listData$tfFamily,
                               organism = x@elementMetadata@listData$organism,
                               motif = x@listData[[1]], alphabet = "DNA",
                               type = "PPM")
              }
              for (i in seq_len(length(motifs))) {
                motifs_out[[i]] <- motifdb_fun(motifs[i])
              }
              motifs <- motifs_out
              success <- TRUE
            }

            # TFBSTools-PFMatrix
            if (in_class_pkg == "TFBSTools" && in_class == "PFMatrix") {
              if (all(names(motifs@bg) %in% c("A", "C", "G", "T"))) {
                alphabet <- "DNA"
              } else alphabet  <- "RNA"
              motifs <- universalmotif(name = motifs@name, altname = motifs@ID,
                                       family = motifs@tags$family,
                                       organism = motifs@tags$species,
                                       motif = motifs@profileMatrix,
                                       alphabet = alphabet, type = "PCM",
                                       bkg = motifs@bg,
                                       strand = paste0(motifs@strand, collapse = ""))
              success <- TRUE
            }

            # TFBSTools-PWMatrix
            if (in_class_pkg == "TFBSTools" && in_class == "PWMatrix") {
              if (all(names(motifs@bg) %in% c("A", "C", "G", "T"))) {
                alphabet <- "DNA"
              } else alphabet  <- "RNA"
              motifs <- universalmotif(name = motifs@name, altname = motifs@ID,
                                       family = motifs@tags$family,
                                       organism = motifs@tags$species,
                                       motif = motifs@profileMatrix,
                                       alphabet = alphabet, type = "PWM",
                                       bkg = motifs@bg,
                                       strand = paste0(motifs@strand, collapse = ""))
              success <- TRUE
            }

            # TFBSTools-ICMatrix
            if (in_class_pkg == "TFBSTools" && in_class == "ICMatrix") {
              if (all(names(motifs@bg) %in% c("A", "C", "G", "T"))) {
                alphabet <- "DNA"
              } else alphabet  <- "RNA"
              motifs <- universalmotif(name = motifs@name, altname = motifs@ID,
                                       family = motifs@tags$family,
                                       organism = motifs@tags$species,
                                       motif = motifs@profileMatrix,
                                       alphabet = alphabet, type = "ICM",
                                       bkg = motifs@bg,
                                       strand = paste0(motifs@strand, collapse = ""))
              success <- TRUE
            }

            # TFBSTools-PFMatrixList and -PWMatrixList and -ICMatrixList
            if (in_class_pkg == "TFBSTools" && in_class %in% c("PFMatrixList",
                                                              "PWMatrixList",
                                                              "ICMatrixList")) {
              motif_num <- length(motifs@listData)
              motifs_out <- list()
              motif_names <- vector(length = motif_num)
              for (i in seq_len(motif_num)) {
                motifs_out[[i]] <- motifs@listData[[i]]
                motif_names[i] <- motifs@listData[[i]]@name
              }
              names(motifs_out) <- motif_names
              motifs <- convert_motifs(motifs_out, class = "universalmotif")
              success <- TRUE
            }

            # seqLogo-pwm
            if (in_class_pkg == "seqLogo" && in_class == "pwm") {
              motifs <- universalmotif(motif = motifs@pwm, type = "PPM",
                                       alphabet = motifs@alphabet)
              success <- TRUE
            }

            # motifStack-pcm
            if (in_class_pkg == "motifStack" && in_class == "pcm") {
              motifs <- universalmotif(name = motifs@name, motif = motifs@mat,
                                       alphabet = motifs@alphabet,
                                       bkg = motifs@background,
                                       type = "PCM")
              success <- TRUE
            }

            # motifStack-pfm
            if (in_class_pkg == "motifStack" && in_class == "pfm") {
              motifs <- universalmotif(name = motifs@name, motif = motifs@mat,
                                       alphabet = motifs@alphabet,
                                       bkg = motifs@background,
                                       type = "PPM")
              success <- TRUE
            }

            # PWMEnrich-PWM
            if (in_class_pkg == "PWMEnrich" && in_class == "PWM") {
              if (all(names(motifs@pwm) %in% c("A", "C", "G", "T"))) {
                alphabet <- "DNA"
              } else alphabet <- "RNA"
              motifs <- universalmotif(name = motifs@name, motif = motifs@pwm,
                                       type = "PWM", alphabet = alphabet,
                                       bkg = motifs@prior.params,
                                       altname = motifs@id)
              success <- TRUE
            }

            if (!success) stop("unrecognized class")

            ## convert to desired class
            motifs <- convert_motifs(motifs, class = class)

            motifs

          })
