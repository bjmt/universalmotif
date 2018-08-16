#' Create a motif.
#'
#' Create a motif from a set of sequences, a matrix, or generate a random
#' motif.
#'
#' @param input \code{character}, \code{numberic}, \code{matrix},
#'              \linkS4class{XStringSet}, \code{missing}
#' @param alphabet \code{character(1)} One of \code{c('DNA', 'RNA', 'AA', 'custom')},
#'    or a combined string representing the letters.
#' @param type \code{character(1)} One of \code{c('PCM', 'PPM', 'PWM', 'ICM')}.
#' @param name \code{character(1)} Motif name.
#' @param pseudocount \code{numeric(1)} Correction to be applied to prevent \code{-Inf}
#'                     from apearing in PWM matrices.
#' @param bkg \code{numeric} Must sum to 1 and be equal in length to the alphabet
#'            length.
#' @param nsites \code{numeric(1)} Number of sites the motif was constructed from.
#' @param altname \code{character(1)} Alternate motif name.
#' @param family \code{character(1)} Transcription factor family.
#' @param organism \code{character(1)} Species of origin.
#' @param bkgsites \code{numeric(1)} Total number of sites used to find the motif.
#' @param strand |code{character(1)} Whether the motif is specific to a certain strand.
#' @param pval \code{numeric(1)} P-value associated with motif.
#' @param qval \code{numeric(1)} Adjusted P-value associated with motif.
#' @param eval \code{numeric(1)} E-value associated with motif.
#' @param extrainfo \code{character} Any other extra information, represented as
#'                  a named character vector.
#' @param add.multifreq \code{numeric(1)} If the motif is created from a set of
#'    sequences, then the \code{\link{add_multifreq}} function can be
#'    run at the same type.
#'
#' @return \linkS4class{universalmotif} object.
#'
#' @details
#'    The aim of this function is provide an easy interface to creating
#'    \linkS4class{universalmotif} motifs, as an alternative to the
#'    default class constructor (i.e. \code{new('universalmotif', name=...)}).
#'    See examples for potential use cases.
#'
#' @seealso \code{\link{convert_type}}, \code{\link{add_multifreq}},
#' \code{\link{create_sequences}}, \code{\link{shuffle_motifs}}.
#'
#' @examples
#' ##### create motifs from a single string
#'
#' # motif is by default generated as a PPM; change final type as desired
#' DNA.motif <- create_motif("TATAWAW")
#' DNA.motif <- create_motif("TATAWAW", type = "PCM")
#'
#' # nsites will be set to the number of input sequences unless specified 
#' DNA.motif <- create_motif("TTTTTTT", nsites = 10)
#'
#' # if ambiguity letters are found and nsites is not specified, nsites will
#' # be set to the minimum required to respect amibiguity letters
#' DNA.motif <- create_motif("TATAWAW")
#' DNA.motif <- create_motif("NNVVWWAAWWDDN")
#'
#' # be careful about setting nsites when using ambiguity letters!
#' DNA.motif <- create_motif("NNVVWWAAWWDDN", nsites = 1)
#'
#' RNA.motif <- create_motif("UUUCCG")
#' 
#' # 'create_motif' will try to detect the alphabet type; this can be 
#' # unreliable for AA and custom alphabets as DNA and RNA alphabets are
#' # detected first
#' AA.motif <- create_motif("AVLK", alphabet = "AA")
#' 
#' custom.motif <- create_motif("QWER", alphabet = "custom")
#' # specify custom alphabet
#' custom.motif <- create_motif("QWER", alphabet = "QWERASDF")
#'
#' ###### create motifs from multiple strings of equal length
#'
#' DNA.motif <- create_motif(c("TTTT", "AAAA", "AACC", "TTGG"), type = "PPM")
#' DNA.motif <- create_motif(c("TTTT", "AAAA", "AACC", "TTGG"), nsites = 20)
#' RNA.motif <- create_motif(c("UUUU", "AAAA", "AACC", "UUGG"), type = "PWM")
#' AA.motif <- create_motif(c("ARNDCQ", "EGHILK", "ARNDCQ"), alphabet = "AA")
#' custom.motif <- create_motif(c("POIU", "LKJH", "POIU", "CVBN"),
#'                              alphabet = "custom")
#'
#' # ambiguity letters are only allowed for single consensus strings; the
#' # following fails
#' \dontrun{
#' create_motif(c("WWTT", "CCGG"))
#' create_motif(c("XXXX", "XXXX"), alphabet = "AA")
#' }
#'
#' ##### create motifs from XStringSet objects
#'
#' library(Biostrings)
#' 
#' DNA.set <- DNAStringSet(c("TTTT", "AAAA", "AACC", "TTGG"))
#' DNA.motif <- create_motif(DNA.set)
#' RNA.set <- RNAStringSet(c("UUUU", "AACC", "UUCC"))
#' RNA.motif <- create_motif(RNA.set)
#' AA.set <- AAStringSet(c("VVVLLL", "AAAIII"))
#' AA.motif <- create_motif(AA.set)
#' 
#' # custom motifs can be created from BStringSet objects
#' B.set <- BStringSet(c("QWER", "ASDF", "ZXCV", "TYUI"))
#' custom.motif <- create_motif(B.set)
#'
#' ##### create motifs with filled 'multifreq' slot
#'
#' DNA.motif.k2 <- create_motif(DNA.set, add.multifreq = 2)
#'
#' ##### create motifs from matrices
#'
#' mat <- matrix(c(1, 1, 1, 1,
#'                 2, 0, 2, 0,
#'                 0, 2, 0, 2,
#'                 0, 0, 0, 0),
#'                 nrow = 4, byrow = TRUE)
#' DNA.motif <- create_motif(mat, alphabet = "DNA")
#' RNA.motif <- create_motif(mat, alphabet = "RNA", nsites = 20)
#' custom.motif <- create_motif(mat)
#'
#' # specify custom alphabet
#' custom.motif <- create_motif(mat, alphabet = "QWER")
#'
#' # alphabet can be detected from rownames
#' rownames(mat) <- DNA_BASES
#' DNA.motif <- create_motif(mat)
#' rownames(mat) <- c("Q", "W", "E", "R")
#' custom.motif <- create_motif(mat)
#'
#' # matrices can also be used as input
#' mat.ppm <- matrix(c(0.1, 0.1, 0.1, 0.1,
#'                     0.5, 0.5, 0.5, 0.5,
#'                     0.1, 0.1, 0.1, 0.1,
#'                     0.3, 0.3, 0.3, 0.3),
#'                     nrow = 4, byrow = TRUE)
#'
#' DNA.motif <- create_motif(mat.ppm, alphabet = "DNA", type = "PPM")
#'
#' ##### create random motifs
#'
#' # these are generated as PPMs with 10 positions
#'
#' DNA.motif <- create_motif()
#' RNA.motif <- create_motif(alphabet = "RNA")
#' AA.motif <- create_motif(alphabet = "AA")
#' custom.motif <- create_motif(alphabet = "QWER")
#'
#' # the number of positions can be specified
#'
#' DNA.motif <- create_motif(5)
#'
#' # if the background frequencies are not provided, they are assumed to be
#' # uniform; if different background frequencies are used, then at each
#' # position \code{rdirichlet(1, bkg)} is used 
#'
#' DNA.motif <- create_motif(bkg = c(0.3, 0.2, 0.2, 0.3))
#' DNA.motif <- create_motif(10, bkg = c(0.1, 0.4, 0.4, 0.1))
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("create_motif", function(input, alphabet, type = "PPM",
                                    name = "motif", pseudocount = 0,
                                    bkg, nsites, altname, family,
                                    organism, bkgsites, strand, pval, qval,
                                    eval, extrainfo, add.multifreq)
           standardGeneric("create_motif"))

#' @describeIn create_motif Create a random motif of length 10.
#' @include universalmotif-class.R
#' @export
setMethod("create_motif", signature(input = "missing"),
          definition = function(input, alphabet, type, name, pseudocount, 
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

            margs <- list()
            if (!missing(nsites)) margs <- c(margs, list(nsites = nsites))
            if (!missing(name)) margs <- c(margs, list(name = name))
            if (!missing(pseudocount)) margs <- c(margs, list(pseudocount = pseudocount))
            if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))
            if (!missing(alphabet)) margs <- c(margs, list(alphabet = alphabet))
            if (!missing(type)) margs <- c(margs, list(type = type))
            if (!missing(altname)) margs <- c(margs, list(altname = altname))
            if (!missing(family)) margs <- c(margs, list(family = family))
            if (!missing(organism)) margs <- c(margs, list(organism = organism))
            if (!missing(bkgsites)) margs <- c(margs, list(bkgsites = bkgsites))
            if (!missing(strand)) margs <- c(margs, list(strand = strand))
            if (!missing(pval)) margs <- c(margs, list(pval = pval))
            if (!missing(qval)) margs <- c(margs, list(qval = qval))
            if (!missing(eval)) margs <- c(margs, list(eval = eval))
            if (!missing(extrainfo)) margs <- c(margs, list(extrainfo = extrainfo))
            
            motif <- do.call(create_motif, c(list(input = 10), margs))
            if (!is.null(motif@motif)) {
              motif@motif <- motif@motif[order(rownames(motif@motif)), ]
            }
            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif

          })

#' @describeIn create_motif Create a random motif with a specified length.
#' @export
setMethod("create_motif", signature(input = "numeric"),
          definition = function(input, alphabet, type, name, pseudocount, 
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

            if (length(input) != 1) stop("input must be a single number")
            if (round(input) != input) stop("input must be a whole number")
            if (input == 0 ) stop("input must be greater than 0")

            if (missing(alphabet)) alphabet <- "DNA"
            if (alphabet %in% c("DNA", "RNA")) {
              alph_len <- 4
            } else if (alphabet == "AA") {
              alph_len <- 20
            } else {
              alph_len <- length(strsplit(alphabet, "")[[1]])
            }

            mot <- matrix(rep(NA, alph_len * input), nrow = alph_len)

            if (missing(bkg)) {
              for (i in seq_len(input)) {
                mot[, i] <- runif(alph_len)
                mot[, i] <- mot[, i] / sum(mot[, i])
              }
            } else {
              if (sum(bkg) < 0.99 || sum(bkg) > 1.01) stop("bkg must sum to 1")
              for (i in seq_len(input)) {
                mot[, i] <- rdirichlet(1, bkg)
              }
            }

            if (missing(type) && missing(nsites)) {
              type <- "PPM"
              nsites <- numeric(0)
            } else if (missing(type)) {
              type <- "PPM"
            } else if (missing(nsites) && type == "PCM") {
              nsites <- 100
            } else nsites <- numeric(0)

            margs <- list(type = type, nsites = nsites)
            if (!missing(name)) margs <- c(margs, list(name = name))
            if (!missing(altname)) margs <- c(margs, list(altname = altname))
            if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))
            if (!missing(pseudocount)) margs <- c(margs, list(pseudocount = pseudocount))
            if (!missing(family)) margs <- c(margs, list(family = family))
            if (!missing(organism)) margs <- c(margs, list(organism = organism))
            if (!missing(bkgsites)) margs <- c(margs, list(bkgsites = bkgsites))
            if (!missing(strand)) margs <- c(margs, list(strand = strand))
            if (!missing(pval)) margs <- c(margs, list(pval = pval))
            if (!missing(qval)) margs <- c(margs, list(qval = qval))
            if (!missing(eval)) margs <- c(margs, list(eval = eval))
            if (!missing(extrainfo)) margs <- c(margs, list(extrainfo = extrainfo))
            
            motif <- do.call(create_motif, c(list(input = mot), margs,
                                             list(alphabet = alphabet)))
            if (!is.null(motif@motif)) {
              motif@motif <- motif@motif[order(rownames(motif@motif)), ]
            }
            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif

          })

#' @describeIn create_motif Create motif from a consensus string.
#' @export
setMethod("create_motif", signature(input = "character"),
          definition = function(input, alphabet, type, name, pseudocount, 
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {
            if (missing(alphabet)) alphabet <- "missing"
            consensus <- input
            consensus.all <- consensus
            if (length(consensus) > 1) {
              consensus <- consensus[1]
            }
            consensus <- strsplit(consensus, split = "")[[1]]
            if (alphabet %in% c("DNA", "RNA") && length(consensus.all) == 1) {
              motif <- vapply(consensus, consensus_to_ppmC, numeric(4))
            } else if (alphabet == "AA" && length(consensus.all) == 1) {
              motif <- vapply(consensus, consensus_to_ppmAAC, numeric(20))
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
                motif <- vapply(consensus, consensus_to_ppmAAC, numeric(20))
                alphabet <- "AA"
              } else if (any(consensus == "U") &&
                         !any(consensus %in% c("E", "F", "I", "J", "L", "O",
                                               "P", "Q", "T", "X", "Z"))) {
                alphabet <- "RNA" 
                motif <- vapply(consensus, consensus_to_ppmC, numeric(4))
              } else if (any(consensus %in% DNA_ALPHABET[-c(16:18)]) &&
                         !any(consensus %in% c("E", "F", "I", "J", "L", "O",
                                               "P", "Q", "X", "Z", "U"))) {
                alphabet <- "DNA"
                motif <- vapply(consensus, consensus_to_ppmC, numeric(4))
              } else if (length(consensus.all) == 1) {
                alphabet <- paste(sort(unique(consensus)), collapse = "")
                motif <- consensusMatrix(paste(consensus, collapse = ""))
              }
            }
            if (alphabet == "AA" && length(consensus.all) > 1) {
              motif2 <- vector("list", 20)
              mot_len <- ncol(motif)
              for (i in AA_STANDARD) {
                motif2[[i]] <- motif[rownames(motif) == i, ]
                if (length(motif2[[i]]) == 0) motif2[[i]] <- rep(0, mot_len)
              }
              motif <- matrix(unlist(motif2), ncol = mot_len, byrow = TRUE)
            }
            margs <- list(name = name, pseudocount = pseudocount)
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
              if (!is.null(motif@motif)) {
                motif@motif <- motif@motif[order(rownames(motif@motif)), ]
              }
              msg <- validObject_universalmotif(motif)
              if (length(msg) > 0) stop(msg)
              return(motif)
            }
            motif <- apply(motif, 2, pcm_to_ppmC, pseudocount = 0)
            motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                               list(alphabet = alphabet),
                                               list(type = "PPM"),
                                               margs))
            if (length(consensus.all) == 1 && missing(nsites) &&
                alphabet %in% c("DNA", "RNA")) {
              input.split <- strsplit(input, "")[[1]]
              if ("N" %in% input.split) {
                motif["nsites"] <- 4
                if (any(c("H", "B", "V", "D") %in% input.split)) {
                  motif["nsites"] <- 12
                }
              } else if (any(c("H", "B", "V", "D") %in% input.split)) {
                motif["nsites"] <- 3
                if (any(c("M", "R", "W", "S", "Y", "K") %in%
                        input.split)) motif["nsites"] <- 6
              } else if (any(c("M", "R", "W", "S", "Y", "K") %in%
                             input.split)) {
                motif["nsites"] <- 2
              }
            } else if (length(consensus.all) == 1 && missing(nsites) &&
                       alphabet == "AA") {
              input.split <- strsplit(input, "")[[1]]
              if ("X" %in% input.split) {
                motif["nsites"] <- 20
              } else if (any(c("B", "Z", "J") %in% input.split)) {
                motif["nsites"] <- 2
              }
            }

            Ncheck <- apply(motif@motif, 2, function(x) all(x == x[1]))
            if (type == "PPM") motif <- convert_type(motif, "PCM")
            if (any(Ncheck)) {
              print(motif)
              motif@motif[, Ncheck] <- rep(1 / nrow(motif@motif),
                                           nrow(motif@motif))
            }
            motif <- convert_type(motif, type = type)

            if (alphabet %in% c("DNA", "RNA") && length(input) > 1 &&
                !missing(add.multifreq)) {
              for (i in add.multifreq) {
                motif@multifreq[[as.character(i)]] <- add_multi(motif@bkg,
                                                                DNAStringSet(input),
                                                                i)
              }
              if (alphabet == "RNA") {
                for (i in names(motif@multifreq)) {
                  rownames(motif@multifreq[[i]]) <- gsub("T", "U",
                                                         rownames(motif@multifreq[[i]]))
                }
              }
            } else if (length(input) > 1 && !missing(add.multifreq)) {
              for (i in add.multifreq) {
                motif@multifreq[[as.character(i)]] <- add_multi_ANY(BStringSet(input),
                                                                    i,
                                                                    rownames(motif["motif"]))
              }
            }

            if (!is.null(rownames(motif@motif))) {
              motif@motif <- motif@motif[order(rownames(motif@motif)), ]
            }

            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif

          })

#' @describeIn create_motif Create motif from a matrix.
#' @export
setMethod("create_motif", signature(input = "matrix"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, 
                                altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {
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

            margs <- list(name = name, pseudocount = pseudocount)
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

            motif <- matrix

            motif <- do.call(universalmotif_cpp, c(list(motif = motif), margs,
                                               list(alphabet = alphabet)))
            if (missing(nsites)) {
              nsites <- sum(input[, 1])
              if (nsites == round(nsites) && nsites != 1 && abs(nsites) != Inf) {
                motif["nsites"] <- nsites
              }
            }
            motif <- convert_type(motif, type = type)
            if (!is.null(rownames(motif@motif))) {
              motif@motif <- motif@motif[order(rownames(motif@motif)), ]
            }
            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif
          })

#' @describeIn create_motif Create motif from a \linkS4class{DNAStringSet}.
#' @export
setMethod("create_motif", signature(input = "DNAStringSet"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {
            sequences <- input

            if (length(unique(width(sequences))) != 1) {
              stop("all sequences must be the same width")
            }

            margs <- list(name = name, pseudocount = pseudocount)
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

            sequences <- consensusMatrix(sequences, baseOnly = TRUE)
            if (sum(sequences[5, ]) > 0) stop("only ACGT are accepted for DNA")
            motif <- apply(sequences[1:4, ], 2, pcm_to_ppmC, pseudocount = 0)
            motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                               list(type = "PPM"),
                                               margs,
                                               list(alphabet = "DNA")))

            if (missing(nsites))  motif["nsites"] <- length(input)
            motif <- convert_type(motif, type = type)
            
            if (length(input) > 1 && !missing(add.multifreq)) {
              for (i in add.multifreq) {
                motif@multifreq[[as.character(i)]] <- add_multi(motif@bkg,
                                                         DNAStringSet(input),
                                                                i)
              }
            }

            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif

          })

#' @describeIn create_motif Create motif from a \linkS4class{RNAStringSet}.
#' @export
setMethod("create_motif", signature(input = "RNAStringSet"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {
            sequences <- input

            if (length(unique(width(sequences))) != 1) {
              stop("all sequences must be the same width")
            }

            margs <- list(name = name, pseudocount = pseudocount)
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

            sequences <- consensusMatrix(sequences, baseOnly = TRUE)
            if (sum(sequences[5, ]) > 0) stop("only ACGU are accepted for RNA")
            motif <- apply(sequences[1:4, ], 2, pcm_to_ppmC, pseudocount = 0)
            motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                               list(type = "PPM"),
                                               margs,
                                               list(alphabet = "RNA")))

            if (missing(nsites))  motif["nsites"] <- length(input)
            motif <- convert_type(motif, type = type)

            if (length(input) > 1 && !missing(add.multifreq)) {
              for (i in add.multifreq) {
                motif@multifreq[[as.character(i)]] <- add_multi(motif@bkg,
                                                        DNAStringSet(input),
                                                             i)
              }
              for (i in names(motif@multifreq)) {
                rownames(motif@multifreq[[i]]) <- gsub("T", "U",
                                                   rownames(motif@multifreq[[i]]))
              }
            }
            
            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif

          })

#' @describeIn create_motif Create motif from a \linkS4class{AAStringSet}.
#' @export
setMethod("create_motif", signature(input = "AAStringSet"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {
            sequences <- input

            if (length(unique(width(sequences))) != 1) {
              stop("all sequences must be the same width")
            }

            margs <- list(name = name, pseudocount = pseudocount)
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

            sequences <- consensusMatrix(sequences)
            if (any(!rownames(sequences) %in% AA_STANDARD)) {
              stop("only ACDEFGHIKLMNPQRSTVWY are accepted for AA")
            }
            motif <- vector("list", 20)
            mot_len <- ncol(sequences)
            for (i in AA_STANDARD) {
              motif[[i]] <- sequences[rownames(sequences) == i, ]
              if (length(motif[[i]]) == 0) motif[[i]] <- rep(0, mot_len)
            }
            motif <- matrix(unlist(motif), ncol = mot_len, byrow = TRUE)
            motif <- apply(motif, 2, pcm_to_ppmC, pseudocount = 0)
            motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                               list(type = "PPM"),
                                               margs,
                                               list(alphabet = "AA")))

            if (length(input) > 1 && !missing(add.multifreq)) {
              for (i in add.multifreq) {
                motif@multifreq[[as.character(i)]] <- add_multi_ANY(input,
                                                                    i,
                                                              AA_STANDARD)
              }
            }

            if (missing(nsites))  motif["nsites"] <- length(input)
            motif <- convert_type(motif, type = type)
            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif

          })

#' @describeIn create_motif Create motif from a \linkS4class{BStringSet}.
#' @export
setMethod("create_motif", signature(input = "BStringSet"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

            sequences <- input

            if (length(unique(width(sequences))) != 1) {
              stop("all sequences must be the same width")
            }

            margs <- list(name = name, pseudocount = pseudocount)
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

            if (missing(alphabet)) {
              sequences <- consensusMatrix(sequences)
              motif <- apply(sequences, 2, pcm_to_ppmC, pseudocount = 0)
              # rownames(motif) <- rownames(sequences)
              alphabet <- paste(rownames(sequences), collapse = "")
              motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                                 list(type = "PPM"),
                                                 margs,
                                                 list(alphabet = alphabet)))
            } else if (alphabet == "custom") {
              sequences <- consensusMatrix(sequences)
              motif <- apply(sequences, 2, pcm_to_ppmC, pseudocount = 0)
              rownames(motif) <- rownames(sequences)
              motif <- do.call(universalmotif_cpp, c(list(motif = motif),
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
              motif <- apply(motif, 2, pcm_to_ppmC, pseudocount = 0)
              motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                                 list(type = "PPM"),
                                                 margs,
                                                 list(alphabet = alphabet)))
            }
          
            if (length(input) > 1 && !missing(add.multifreq)) {
              for (i in add.multifreq) {
                motif@multifreq[[as.character(i)]] <- add_multi_ANY(input,
                                                                    i,
                                                 rownames(motif["motif"]))
              }
            }

            if (missing(nsites))  motif["nsites"] <- length(input)
            motif <- convert_type(motif, type = type)
            if (!is.null(motif@motif) && !is.null(rownames(motif@motif))) {
              motif@motif <- motif@motif[order(rownames(motif@motif)), ]
            }
            msg <- validObject_universalmotif(motif)
            if (length(msg) > 0) stop(msg)
            motif

          })
