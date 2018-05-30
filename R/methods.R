#' @describeIn create_motif Create motif from a consensus string.
#' @include universalmotif-class.R
#' @export
setMethod("create_motif", signature(input = "character"),
          definition = function(input, alphabet, type, name, pseudoweight, 
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo) {
            if (missing(alphabet)) alphabet <- "missing"
            consensus <- input
            consensus.all <- consensus
            if (length(consensus) > 1) {
              consensus <- consensus[1]
            }
            consensus <- strsplit(consensus, split = "")[[1]]
            if (alphabet %in% c("DNA", "RNA") && length(consensus.all) == 1) {
              motif <- vapply(consensus, consensus_to_ppm, numeric(4))
            } else if (alphabet == "AA" && length(consensus.all) == 0) {
              motif <- vapply(consensus, consensus_to_ppmAA, numeric(20))
            } else if (!missing(alphabet)) {
              motif <- consensusMatrix(paste(consensus, collapse = ""))
            }
            if (!alphabet %in% c("DNA", "RNA", "AA", "custom", "missing") &&
                length(consensus.all) == 1) {
              alph.deparsed <- strsplit(alphabet, "")[[1]]
              if (any(!consensus %in% alph.deparsed)) {
                stop("consensus string does not match provided alphabet")
              }
              motif2 <- vector("list", length(alph.deparsed))
              mot_len <- length(consensus)
              for (i in alph.deparsed) {
                motif2[[i]] <- motif[rownames(motif) == i, ]
                if (length(motif2[[i]]) == 0) motif2[[i]] <- rep(0, mot_len)
              }
              motif <- matrix(unlist(motif2), ncol = mot_len, byrow = TRUE)
            }
            if (alphabet == "missing") {
              if (any(consensus %in% c("E", "F", "I", "P", "Q", "X", "Z")) &&
                  !any(consensus %in% c("O", "U"))) {
                motif <- vapply(consensus, consensus_to_ppmAA, numeric(20))
                alphabet <- "AA"
              } else if (any(consensus == "U") &&
                         !any(consensus %in% c("E", "F", "I", "J", "L", "O",
                                               "P", "Q", "T", "X", "Z"))) {
                alphabet <- "RNA" 
                motif <- vapply(consensus, consensus_to_ppm, numeric(4))
              } else if (any(consensus %in% DNA_ALPHABET[-c(16:18)]) &&
                         !any(consensus %in% c("E", "F", "I", "J", "L", "O",
                                               "P", "Q", "X", "Z", "U"))) {
                alphabet <- "DNA"
                motif <- vapply(consensus, consensus_to_ppm, numeric(4))
              } else if (length(consensus.all) == 1) {
                stop("cannot create a motif using a single consensus string without an alphabet")
              }
            }
            margs <- list(name = name, pseudoweight = pseudoweight)
            if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))
            if (!missing(nsites)) {
              margs <- c(margs, list(nsites = nsites))
            } else {
              margs <- c(margs, list(nsites = length(consensus.all)))
            }
            if (!missing(altname)) margs <- c(margs, list(altname = altname))
            if (!missing(family)) margs <- c(margs, list(family = family))
            if (!missing(organism)) margs <- c(margs, list(organism = organism))
            if (!missing(bkgsites)) margs <- c(margs, list(bkgsites = bkgsites))
            if (!missing(strand)) margs <- c(margs, list(strand = strand))
            if (!missing(pval)) margs <- c(margs, list(pval = pval))
            if (!missing(qval)) margs <- c(margs, list(qval = qval))
            if (!missing(eval)) margs <- c(margs, list(eval = eval))
            if (!missing(extrainfo)) margs <- c(margs, list(extrainfo = extrainfo))
            if (length(consensus.all) > 1) {
              if (alphabet == "DNA") {
                consensus <- lapply(consensus.all, DNAString)
                consensus <- DNAStringSet(consensus)
              } else if (alphabet == "RNA") {
                consensus <- lapply(consensus.all, RNAString)
                consensus <- RNAStringSet(consensus)
              } else if (alphabet == "AA") {
                consensus <- lapply(consensus.all, AAString)
                consensus <- AAStringSet(consensus)
              } else {
                consensus <- lapply(consensus.all, BString)
                consensus <- BStringSet(consensus)
                if (alphabet != "custom") {
                  alph.deparsed <- strsplit(alphabet, "")[[1]]
                  if (any(!rownames(consensusMatrix(consensus)) %in%
                          alph.deparsed)) {
                    stop("consensus string does not match provided alphabet")
                  }
                } else {
                  alphabet <- paste(rownames(consensusMatrix(consensus)),
                                    collapse = "")
                }
              }
              if (!missing(type)) margs <- c(margs, list(type = type))
              motif <- do.call(create_motif,
                               c(list(input = consensus), margs,
                                 list(alphabet = alphabet)))
              return(motif)
            }
            motif <- apply(motif, 2, pcm_to_ppm, pseudoweight = 0)
            motif <- do.call(universalmotif, c(list(motif = motif),
                                               list(alphabet = alphabet),
                                               list(type = "PPM"),
                                               margs))
            if (!missing(type)) {
              motif <- convert_type(motif, type = type)
            } else {
              motif <- convert_type(motif, type = "PCM")
            }
            motif
          })

#' @describeIn create_motif Create motif from a matrix.
#' @export
setMethod("create_motif", signature(input = "matrix"),
          definition = function(input, alphabet, type, name, pseudoweight,
                                bkg, nsites, 
                                altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo) {
            matrix <- input
            if (!missing(alphabet) &&
                !alphabet %in% c("DNA", "RNA", "AA", "custom")) {
              alph.deparsed <- strsplit(alphabet, "")[[1]]
              if (any(!rownames(matrix) %in% alph.deparsed)) {
                stop("rownames do not match provided alphabet")
              }
              if (length(alph.deparsed) != nrow(matrix)) {
                stop("alphabet length does not match number of rows")
              }
            } else if (is.null(rownames(matrix)) && missing(alphabet)) {
              alphabet <- "custom"
            } else if (all(rownames(matrix) %in% DNA_BASES) &&
                       missing(alphabet) && nrow(matrix) == 4) {
              alphabet  <- "DNA"
            } else if (all(rownames(matrix) %in% RNA_BASES) &&
                       missing(alphabet) && nrow(matrix) == 4) {
              alphabet <- "RNA"
            } else if (nrow(matrix) == 20 && missing(alphabet)) {
              alphabet <- "AA" 
            } else if (all(rownames(matrix) %in% AA_STANDARD) &&
                       missing(alphabet) && nrow(matrix) == 20) {
              alphabet <- "AA"
            } else if (!is.null(rownames(matrix))) {
              alphabet <- paste(rownames(matrix), collapse = "") 
            } else if (nrow(matrix == 4) && missing(alphabet)) {
              alphabet <- "DNA" 
            } else if (missing(alphabet)) {
              alphabet <- paste(rownames(matrix), collapse = "")
            }

            if (alphabet %in% c("DNA", "RNA")) {
              if (nrow(matrix) != 4) {
                stop("incorrect number of rows")
              }
            } else if (alphabet == "AA") {
              if (nrow(matrix) != 20) {
                stop("incorrect number of rows")
              }
            }

            margs <- list(name = name, pseudoweight = pseudoweight)
            if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))
            if (!missing(nsites)) margs <- c(margs, list(nsites = nsites))
            if (!missing(altname)) margs <- c(margs, list(altname = altname))
            if (!missing(family)) margs <- c(margs, list(family = family))
            if (!missing(organism)) margs <- c(margs, list(organism = organism))
            if (!missing(bkgsites)) margs <- c(margs, list(bkgsites = bkgsites))
            if (!missing(strand)) margs <- c(margs, list(strand = strand))
            if (!missing(pval)) margs <- c(margs, list(pval = pval))
            if (!missing(qval)) margs <- c(margs, list(qval = qval))
            if (!missing(eval)) margs <- c(margs, list(eval = eval))
            if (!missing(extrainfo)) margs <- c(margs, list(extrainfo = extrainfo))

            motif <- apply(matrix, 2, pcm_to_ppm, pseudoweight = 0)

            motif <- do.call(universalmotif, c(list(motif = motif), margs,
                                               list(type = "PPM"),
                                               list(alphabet = alphabet)))
            if (missing(nsites))  motif["nsites"] <- sum(input[, 1])
            if (!missing(type)) {
              motif <- convert_type(motif, type = type)
            } else {
              motif <- convert_type(motif, type = "PCM")
            }
            motif
          })

#' @describeIn create_motif Create motif from an XStringSet object.
#' @export
setMethod("create_motif", signature(input = "XStringSet"),
          definition = function(input, alphabet, type, name, pseudoweight, 
                                bkg, nsites,
                                altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo) {

            sequences <- input

            if (length(unique(width(sequences))) != 1) {
              stop("all sequences must be the same width")
            }

            margs <- list(name = name, pseudoweight = pseudoweight)
            if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))
            if (!missing(nsites)) margs <- c(margs, list(nsites = nsites))
            if (!missing(altname)) margs <- c(margs, list(altname = altname))
            if (!missing(family)) margs <- c(margs, list(family = family))
            if (!missing(organism)) margs <- c(margs, list(organism = organism))
            if (!missing(bkgsites)) margs <- c(margs, list(bkgsites = bkgsites))
            if (!missing(strand)) margs <- c(margs, list(strand = strand))
            if (!missing(pval)) margs <- c(margs, list(pval = pval))
            if (!missing(qval)) margs <- c(margs, list(qval = qval))
            if (!missing(eval)) margs <- c(margs, list(eval = eval))
            if (!missing(extrainfo)) margs <- c(margs, list(extrainfo = extrainfo))

            alph <- sequences@elementType
            if (alph == "DNAString") {
              sequences <- consensusMatrix(sequences, baseOnly = TRUE)
              if (sum(sequences[5, ]) > 0) stop("only ACGT are accepted for DNA")
              motif <- apply(sequences[1:4, ], 2, pcm_to_ppm, pseudoweight = 0)
              motif <- do.call(universalmotif, c(list(motif = motif),
                                                 list(type = "PPM"),
                                                 margs,
                                                 list(alphabet = "DNA")))
            } else if (alph == "RNAString") {
              sequences <- consensusMatrix(sequences, baseOnly = TRUE)
              if (sum(sequences[5, ]) > 0) stop("only ACGU are accepted for RNA")
              motif <- apply(sequences[1:4, ], 2, pcm_to_ppm, pseudoweight = 0)
              motif <- do.call(universalmotif, c(list(motif = motif),
                                                 list(type = "PPM"),
                                                 margs,
                                                 list(alphabet = "RNA")))
            } else if (alph == "AAString") {
              sequences <- consensusMatrix(sequences)
              motif <- vector("list", 20)
              mot_len <- ncol(sequences)
              for (i in AA_STANDARD) {
                motif[[i]] <- sequences[rownames(sequences) == i, ]
                if (length(motif[[i]]) == 0) motif[[i]] <- rep(0, mot_len)
              }
              motif <- matrix(unlist(motif), ncol = mot_len, byrow = TRUE)
              motif <- apply(motif, 2, pcm_to_ppm, pseudoweight = 0)
              motif <- do.call(universalmotif, c(list(motif = motif),
                                                 list(type = "PPM"),
                                                 margs,
                                                 list(alphabet = "AA")))
            } else if (alph == "custom" || missing(alphabet)) {
              sequences <- consensusMatrix(sequences)
              motif <- apply(sequences, 2, pcm_to_ppm, pseudoweight = 0)
              motif <- do.call(universalmotif, c(list(motif = motif),
                                                 list(type = "PPM"),
                                                 margs,
                                                 list(alphabet = "custom")))
            } else {
              sequences <- consensusMatrix(sequences)
              alph.split <- strsplit(alphabet, "")[[1]]
              motif <- vector("list", length(alph.split))
              mot_len <- ncol(sequences)
              for (i in alph.split) {
                motif[[i]] <- sequences[rownames(sequences) == i, ]
                if (length(motif[[i]]) == 0) motif[[i]] <- rep(0, mot_len)
              }
              motif <- matrix(unlist(motif), ncol = mot_len, byrow = TRUE)
              motif <- apply(motif, 2, pcm_to_ppm, pseudoweight = 0)
              motif <- do.call(universalmotif, c(list(motif = motif),
                                                 list(type = "PPM"),
                                                 margs,
                                                 list(alphabet = alphabet)))
            }

            if (missing(nsites))  motif["nsites"] <- length(input)
            if (!missing(type)) {
              motif <- convert_type(motif, type = type)
            } else {
              motif <- convert_type(motif, type = "PCM")
            }
            motif

          })

#' @describeIn convert_motifs Convert a list of motifs.
#' @export
setMethod("convert_motifs", signature(motifs = "list"),
          definition = function(motifs, class, BPPARAM) {
            bplapply(motifs, function(x) convert_motifs(x, class = class,
                                                        BPPARAM = BPPARAM))
          })

#' @describeIn convert_motifs Convert a universalmotif object.
#' @export
setMethod("convert_motifs", signature(motifs = "universalmotif"),
          definition = function(motifs, class, BPPARAM) {
            
            out_class <- strsplit(class, "-")[[1]][2]
            out_class_pkg <- strsplit(class, "-")[[1]][1]

            if (out_class_pkg == "universalmotif") {
              return(motifs)
            }

            # MotIV-pwm2
            if (out_class_pkg == "MotIV" && out_class == "pwm2") {
              motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)
              motifs <- makePWM(motifs["motif"],
                                alphabet = motifs["alphabet"])
              return(motifs)
            }

            # TFBSTools- PFMatrix, PWMatrix, and ICMatrix
            if (out_class_pkg == "TFBSTools") {
              motifs <- convert_type(motifs, "PCM")
              bkg <- motifs["bkg"]
              names(bkg) <- DNA_BASES
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
              rownames(bio_mat) <- DNA_BASES
              bio_priors <- motifs["bkg"]
              names(bio_priors) <- DNA_BASES
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
              motifs <- convert_type(motifs, "PCM")
              bio_mat <- matrix(as.integer(motifs["motif"]), byrow = FALSE,
                                nrow = 4)
              rownames(bio_mat) <- DNA_BASES
              bio_priors <- motifs["bkg"]
              names(bio_priors) <- DNA_BASES
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
#' @export
setMethod("convert_motifs", signature(motifs = "ANY"),
          definition = function(motifs, class, BPPARAM) {
          
            success <- FALSE

            ## convert to universalmotif
            in_class <- class(motifs)[1]
            in_class_pkg <- attributes(class(motifs))$package

            if (paste(in_class_pkg, in_class, sep = "-") == class) {
              return(motifs)
            }

            # MotifDb-MotifList
            if (in_class_pkg == "MotifDb" && in_class == "MotifList") {
              motifdb_fun <- function(i) {
                x <- motifs[i]
                universalmotif(name = x@elementMetadata@listData$providerName,
                               altname = x@elementMetadata@listData$geneSymbol,
                               family = x@elementMetadata@listData$tfFamily,
                               organism = x@elementMetadata@listData$organism,
                               motif = x@listData[[1]], alphabet = "DNA",
                               type = "PPM")
              }
              motifs_out <- vector("list", length(motifs))
              motifs_out <- bplapply(seq_len(length(motifs)), motifdb_fun,
                                     BPPARAM = BPPARAM)
              motifs <- motifs_out
              success <- TRUE
            }

            # TFBSTools-PFMatrix
            if (in_class_pkg == "TFBSTools" && in_class == "PFMatrix") {
              if (all(names(motifs@bg) %in% DNA_BASES)) {
                alphabet <- "DNA"
              } else alphabet  <- "RNA"
              motifs <- universalmotif(name = motifs@name, altname = motifs@ID,
                                       family = motifs@tags$family,
                                       organism = motifs@tags$species,
                                       motif = motifs@profileMatrix,
                                       alphabet = alphabet, type = "PCM",
                                       bkg = motifs@bg,
                                       strand = paste0(motifs@strand,
                                                       collapse = ""))
              success <- TRUE
            }

            # TFBSTools-PWMatrix
            if (in_class_pkg == "TFBSTools" && in_class == "PWMatrix") {
              if (all(names(motifs@bg) %in% DNA_BASES)) {
                alphabet <- "DNA"
              } else alphabet  <- "RNA"
              motifs <- universalmotif(name = motifs@name, altname = motifs@ID,
                                       family = motifs@tags$family,
                                       organism = motifs@tags$species,
                                       motif = motifs@profileMatrix,
                                       alphabet = alphabet, type = "PWM",
                                       bkg = motifs@bg,
                                       strand = paste0(motifs@strand,
                                                       collapse = ""))
              success <- TRUE
            }

            # TFBSTools-ICMatrix
            if (in_class_pkg == "TFBSTools" && in_class == "ICMatrix") {
              if (all(names(motifs@bg) %in% DNA_BASES)) {
                alphabet <- "DNA"
              } else alphabet  <- "RNA"
              motifs <- universalmotif(name = motifs@name, altname = motifs@ID,
                                       family = motifs@tags$family,
                                       organism = motifs@tags$species,
                                       motif = motifs@profileMatrix,
                                       alphabet = alphabet, type = "ICM",
                                       bkg = motifs@bg,
                                       strand = paste0(motifs@strand,
                                                       collapse = ""))
              success <- TRUE
            }

            # TFBSTools-PFMatrixList and -PWMatrixList and -ICMatrixList
            if (in_class_pkg == "TFBSTools" && in_class %in% c("PFMatrixList",
                                                              "PWMatrixList",
                                                              "ICMatrixList")) {
              motif_num <- length(motifs@listData)
              motifs_out <- bplapply(seq_len(motif_num),
                                     function(i) {
                                       motifs@listData[[i]]
                                     }, BPPARAM = BPPARAM)
              motif_names <- unlist(bplapply(seq_len(motif_num),
                                             function(i) {
                                              motifs@listData[[i]]@name
                                             }, BPPARAM = BPPARAM))
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
              if (all(names(motifs@pwm) %in% DNA_BASES)) {
                alphabet <- "DNA"
              } else alphabet <- "RNA"
              motifs <- universalmotif(name = motifs@name, motif = motifs@pwm,
                                       type = "PWM", alphabet = alphabet,
                                       bkg = motifs@prior.params,
                                       altname = motifs@id)
              success <- TRUE
            }

            # motifRG-Motif
            if (in_class_pkg == "motifRG" && in_class == "Motif") {
              motifs <- universalmotif(name = motifs@pattern,
                                       nsites = sum(motifs@count),
                                       alphabet = "DNA",
                                       type = "PCM",
                                       extrainfo = c(score = motifs@score),
                               strand = paste(unique(motifs@match$match.strand),
                                              collapse = ""),
          motif = create_motif(input = DNAStringSet(motifs@match$pattern)))
              success <- TRUE
            }

            paste(in_class_pkg)
            paste(in_class)
            if (!success) stop("unrecognized class")

            ## convert to desired class
            motifs <- convert_motifs(motifs, class = class)

            motifs

          })
