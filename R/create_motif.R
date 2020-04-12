#' Create a motif.
#'
#' Create a motif from a set of sequences, a matrix, or generate a random
#' motif. See the "Motif import, export and manipulation" vignette for details.
#'
#' @param input `character`, `numeric`, `matrix`,
#'    \code{\link{XStringSet}}, or `missing`.
#' @param alphabet `character(1)` One of `c('DNA', 'RNA', 'AA')`,
#'    or a combined string representing the letters. If no alphabet is
#'    provided then it will try and guess the alphabet from the input.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#' @param name `character(1)` Motif name.
#' @param pseudocount `numeric(1)` Correction to be applied to prevent `-Inf`
#'   from appearing in PWM matrices. Defaults to 0.
#' @param bkg `numeric` A vector of probabilities, each between 0 and 1. If
#'    higher order backgrounds are provided, then the elements of the vector
#'    must be named. If unnamed, then the order of probabilities must be in the
#'    same order as the alphabetically sorted sequence alphabet.
#' @param nsites `numeric(1)` Number of sites the motif was constructed from. If
#'    blank, then `create_motif()` will guess the appropriate number if possible.
#'    To prevent this, provide `nsites = numeric()`.
#' @param altname `character(1)` Alternate motif name.
#' @param family `character(1)` Transcription factor family.
#' @param organism `character(1)` Species of origin.
#' @param bkgsites `numeric(1)` Total number of sites used to find the motif.
#' @param strand `character(1)` Whether the motif is specific to a certain strand.
#' @param pval `numeric(1)` P-value associated with motif.
#' @param qval `numeric(1)` Adjusted P-value associated with motif.
#' @param eval `numeric(1)` E-value associated with motif.
#' @param extrainfo `character` Any other extra information, represented as
#'    a named character vector.
#' @param add.multifreq `numeric` If the motif is created from a set of
#'    sequences, then the [add_multifreq()] function can be
#'    run at the same time (with `RC = FALSE`).
#'
#' @return [universalmotif-class] object.
#'
#' @details
#'    The aim of this function is provide an easy interface to creating
#'    [universalmotif-class] motifs, as an alternative to the
#'    default class constructor (i.e. `new('universalmotif', name=...)`).
#'    See examples for potential use cases.
#'
#'    Note: when generating random motifs, the `nsites` slot is also given a
#'    random value. Furthermore, be careful about the `nsites` slot when creating
#'    motifs from consensus strings: for example, the following call
#'    `create_motif("TAAAT")` generates a motif with `nsites = 1`.
#'
#'    See the `examples` section for more info on motif creation.
#'
#' @seealso [convert_type()], [add_multifreq()], [create_sequences()],
#'    [shuffle_motifs()].
#'
#' @examples
#' ##### create motifs from a single string
#'
#' # Motif is by default generated as a PPM: change final type as desired
#' DNA.motif <- create_motif("TATAWAW")
#' DNA.motif <- create_motif("TATAWAW", type = "PCM")
#'
#' # Nsites will be set to the number of input sequences unless specified
#' DNA.motif <- create_motif("TTTTTTT", nsites = 10)
#'
#' # If ambiguity letters are found and nsites is not specified, nsites will
#' # be set to the minimum required to respect amibiguity letters
#' DNA.motif <- create_motif("TATAWAW")
#' DNA.motif <- create_motif("NNVVWWAAWWDDN")
#'
#' # Be careful about setting nsites when using ambiguity letters!
#' DNA.motif <- create_motif("NNVVWWAAWWDDN", nsites = 1)
#'
#' RNA.motif <- create_motif("UUUCCG")
#'
#' # 'create_motif' will try to detect the alphabet type; this can be
#' # unreliable for AA and custom alphabets as DNA and RNA alphabets are
#' # detected first
#' AA.motif <- create_motif("AVLK", alphabet = "AA")
#'
#' custom.motif <- create_motif("QWER", alphabet = "QWER")
#' # Specify custom alphabet
#' custom.motif <- create_motif("QWER", alphabet = "QWERASDF")
#'
#' ###### Create motifs from multiple strings of equal length
#'
#' DNA.motif <- create_motif(c("TTTT", "AAAA", "AACC", "TTGG"), type = "PPM")
#' DNA.motif <- create_motif(c("TTTT", "AAAA", "AACC", "TTGG"), nsites = 20)
#' RNA.motif <- create_motif(c("UUUU", "AAAA", "AACC", "UUGG"), type = "PWM")
#' AA.motif <- create_motif(c("ARNDCQ", "EGHILK", "ARNDCQ"), alphabet = "AA")
#' custom.motif <- create_motif(c("POIU", "LKJH", "POIU", "CVBN"),
#'                              alphabet = "POIULKJHCVBN")
#'
#' # Ambiguity letters are only allowed for single consensus strings: the
#' # following fails
#' \dontrun{
#' create_motif(c("WWTT", "CCGG"))
#' create_motif(c("XXXX", "XXXX"), alphabet = "AA")
#' }
#'
#' ##### Create motifs from XStringSet objects
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
#' # Custom motifs can be created from BStringSet objects
#' B.set <- BStringSet(c("QWER", "ASDF", "ZXCV", "TYUI"))
#' custom.motif <- create_motif(B.set)
#'
#' ##### Create motifs with filled 'multifreq' slot
#'
#' DNA.motif.k2 <- create_motif(DNA.set, add.multifreq = 2)
#'
#' ##### Create motifs from matrices
#'
#' mat <- matrix(c(1, 1, 1, 1,
#'                 2, 0, 2, 0,
#'                 0, 2, 0, 2,
#'                 0, 0, 0, 0),
#'                 nrow = 4, byrow = TRUE)
#' DNA.motif <- create_motif(mat, alphabet = "DNA")
#' RNA.motif <- create_motif(mat, alphabet = "RNA", nsites = 20)
#' custom.motif <- create_motif(mat, alphabet = "QWER")
#'
#' # Specify custom alphabet
#' custom.motif <- create_motif(mat, alphabet = "QWER")
#'
#' # Alphabet can be detected from rownames
#' rownames(mat) <- DNA_BASES
#' DNA.motif <- create_motif(mat)
#' rownames(mat) <- c("Q", "W", "E", "R")
#' custom.motif <- create_motif(mat)
#'
#' # Matrices can also be used as input
#' mat.ppm <- matrix(c(0.1, 0.1, 0.1, 0.1,
#'                     0.5, 0.5, 0.5, 0.5,
#'                     0.1, 0.1, 0.1, 0.1,
#'                     0.3, 0.3, 0.3, 0.3),
#'                     nrow = 4, byrow = TRUE)
#'
#' DNA.motif <- create_motif(mat.ppm, alphabet = "DNA", type = "PPM")
#'
#' ##### Create random motifs
#'
#' # These are generated as PPMs with 10 positions
#'
#' DNA.motif <- create_motif()
#' RNA.motif <- create_motif(alphabet = "RNA")
#' AA.motif <- create_motif(alphabet = "AA")
#' custom.motif <- create_motif(alphabet = "QWER")
#'
#' # The number of positions can be specified
#'
#' DNA.motif <- create_motif(5)
#'
#' # If the background frequencies are not provided, they are generated
#' # using `rpois`; positions are created using `rdirichlet(1, bkg)`.
#' # (calling `create_motif()` creates motifs with an average
#' # positional IC of 1)
#'
#' DNA.motif <- create_motif(bkg = c(0.3, 0.2, 0.2, 0.3))
#' DNA.motif <- create_motif(10, bkg = c(0.1, 0.4, 0.4, 0.1))
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [create_sequences()]
#' @export
setGeneric("create_motif", function(input, alphabet, type = "PPM",
                                    name = "motif", pseudocount = 0,
                                    bkg, nsites, altname, family,
                                    organism, bkgsites, strand, pval, qval,
                                    eval, extrainfo, add.multifreq)
           standardGeneric("create_motif"))

# TODO: Organise the methods better, lots of repeat code

#' @describeIn create_motif Create a random motif of length 10.
#' @include universalmotif-class.R
#' @export
setMethod("create_motif", signature(input = "missing"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

  margs <- parse_args(as.list(environment()),
                      c("nsites", "name", "pseudocount", "bkg", "alphabet",
                        "type", "altname", "family", "organism", "bkgsites",
                        "strand", "pval", "qval", "eval", "extrainfo",
                        "add.multifreq"))

  motif <- do.call(create_motif, c(list(input = 10), margs))

  validObject_universalmotif(motif)
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
  if (as.integer(input) != input) stop("'input' must be a whole number")
  if (input <= 0 ) stop("input must be greater than 0")

  if (missing(alphabet)) alphabet <- "DNA"

  alph.split <- switch(alphabet, "DNA" = DNA_BASES,
                       "RNA" = RNA_BASES, "AA" = AA_STANDARD2,
                       "custom" = stop(wmsg("`alphabet = 'custom'` is no longer",
                                            " acceptable; please provide the ",
                                            "actual letters")),
                        safeExplode(alphabet))
  alph_len <- length(alph.split)

  if (missing(bkg)) {

    bkg <- rpois(alph_len, 1000 / alph_len) / 1000
    bkg <- bkg / sum(bkg)
    mot <- generate_motif(input, bkg)

  } else {

    if (is.null(names(bkg)) && length(bkg) > alph_len)
      stop("please provide a named vector for 'bkg'")
    else if (is.null(names(bkg)) && length(bkg) == alph_len)
      names(bkg) <- alph.split

    if (length(bkg) < alph_len)
      stop("'bkg' must be at least ", alph_len, " elements long")

    mot <- generate_motif(input, bkg[alph.split])

  }

  if (missing(type) && missing(nsites)) {
    type <- "PPM"
    nsites <- sample.int(201, 1) + 49
  } else if (missing(type)) type <- "PPM"
  else if (missing(nsites) && type == "PCM") nsites <- 100
  else nsites <- numeric(0)

  margs <- parse_args(as.list(environment()),
                      c("name", "altname", "bkg", "pseudocount", "family",
                        "organism", "bkgsites", "strand", "pval", "qval",
                        "eval", "extrainfo"))
  margs <- c(margs, list(type = type), list(nsites = nsites))

  motif <- do.call(create_motif, c(list(input = mot), margs,
                                   list(alphabet = alphabet)))

  validObject_universalmotif(motif)

  if (!missing(add.multifreq)) {
    motif <- add_multifreq(motif, sample_sites(motif, motif@nsites),
                           add.multifreq)
  }

  motif

})

#' @describeIn create_motif Create motif from a consensus string.
#' @export
setMethod("create_motif", signature(input = "character"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

  lets.uniq <- sort_unique_cpp(safeExplode(collapse_cpp(input)))
  if (!missing(alphabet) && length(alphabet) > 1)
    stop(wmsg("`alphabet` parameter should be a single string"), call. = FALSE)

  if (length(unique(nchar(input))) > 1) {
    stop(wmsg("All strings must be of equal length"))
  }
  if (any(nchar(input) == 0)) {
    stop(wmsg("Zero-length input strings as not allowed"))
  }

  if (missing(alphabet)) {
    if (any(tolower(lets.uniq) != lets.uniq) &&
        any(toupper(lets.uniq) != lets.uniq)) {
      message(wmsg(
        "Note: detected both lower and upper-case letters. These ",
        "will be treated as separate letters."
      ))
    }
    if (all(lets.uniq %in% DNA_ALPHABET[-(16:18)])) {
      alphabet <- "DNA"
    } else if (all(lets.uniq %in% RNA_BASES)) {
      alphabet <- "RNA"
    } else if (all(lets.uniq %in% AA_ALPHABET[-(27:30)])) {
      alphabet <- "AA"
    } else {
      message(wmsg(
        "Note: failed to auto-detect alphabet type, creating custom alphabet ",
        "based on input letters."
      ))
      alphabet <- collapse_cpp(lets.uniq)
    }
  } else if (alphabet == "custom") {
    stop(wmsg("`alphabet = 'custom'` is no longer acceptable; please provide ",
              "the actual letters"))
  } else if (alphabet %in% c("DNA", "RNA", "AA")) {
    switch(alphabet,
      "DNA" = if (any(!lets.uniq %in% DNA_BASES)) {
        stop(wmsg(
          "Detected non-DNA letters in input, despite `alphabet='DNA'` [",
          collapse_cpp(lets.uniq[!lets.uniq %in% DNA_BASES]), "]"
        ))
      },
      "RNA" = if (any(!lets.uniq %in% RNA_BASES)) {
        stop(wmsg(
          "Detected non-RNA letters in input, despite `alphabet='RNA'` [",
          collapse_cpp(lets.uniq[!lets.uniq %in% RNA_BASES]), "]"
        ))
      },
      "AA" = if (any(!lets.uniq %in% AA_STANDARD2)) {
        stop(wmsg(
          "Detected non-AA letters in input, despite `alphabet='AA'` [",
          collapse_cpp(lets.uniq[!lets.uniq %in% AA_STANDARD2]), "]"
        ))
      }
    )
  } else {
    if (nchar(alphabet) == 1) stop("alphabet string must be longer than 1 character")
    alphabet <- collapse_cpp(sort_unique_cpp(safeExplode(alphabet)))
  }

  if (length(input) == 1) {
    motif <- switch(alphabet,
      "DNA" = vapply(safeExplode(input), consensus_to_ppmC, numeric(4)),
      "RNA" = vapply(safeExplode(input), consensus_to_ppmC, numeric(4)),
      "AA" = vapply(safeExplode(input), consensus_to_ppmAAC, numeric(20)),
      {
        out <- consensusMatrix(input)
        if (nrow(out) != nchar(alphabet)) {
          missing.lets <- safeExplode(alphabet)
          missing.lets <- missing.lets[!missing.lets %in% rownames(out)]
          out2 <- rbind(out, matrix(0, ncol = ncol(out), nrow = length(missing.lets)))
          rownames(out2) <- c(rownames(out), missing.lets)
          out <- out2[sort_unique_cpp(rownames(out2)), ]
        }
        out
      }
    )
  } else {
    motif <- switch(alphabet,
      "DNA" = DNAStringSet(input),
      "RNA" = RNAStringSet(input),
      "AA" = AAStringSet(input),
      BStringSet(input)
    )
    if (alphabet %in% c("DNA", "RNA", "AA")) {
      switch(alphabet,
        "DNA" = if (sum(consensusMatrix(motif, baseOnly = TRUE)[5, ]))
          stop(wmsg("DNA ambiguity letters are not accepted if input is more than one string")),
        "RNA" = if (sum(consensusMatrix(motif, baseOnly = TRUE)[5, ]))
          stop(wmsg("RNA ambiguity letters are not accepted if input is more than one string")),
        "AA" = if (any(lets.uniq %in% c("U", "O", "B", "J", "Z", "X")))
          stop(wmsg("AA ambiguity letters are not accepted if input is more than one string"))
      )
    }
  }

  margs <- parse_args(
    as.list(environment()),
    c("bkg", "altname", "family", "organism", "bkgsites", "strand", "pval",
      "qval", "eval", "extrainfo")
  )
  margs <- c(margs, list(name = name), list(pseudocount = pseudocount))
  if (!missing(nsites)) {
    margs <- c(margs, list(nsites = nsites))
  } else {
    if (length(input) > 1) margs$nsites <- length(input)
  }

  if (is.matrix(motif)) {
    motif <- do.call(universalmotif_cpp, c(list(motif = motif, alphabet = alphabet, type = "PPM"), margs))
  } else {
    motif <- do.call(create_motif, c(list(input = motif, alphabet = alphabet, type = "PPM"), margs))
  }

  if (type == "PPM") motif <- convert_type_internal(motif, "PCM")
  motif <- convert_type_internal(motif, type = type)

  if (!missing(add.multifreq) && length(input) > 1) {

    for (i in add.multifreq) {
      motif@multifreq[[as.character(i)]] <- add_multi_cpp(input, i, rownames(motif@motif))
    }

    new.bkg <- get_bkg(BStringSet(input), k = add.multifreq,
                       pseudocount = pseudocount, list.out = FALSE,
                       alphabet = rownames(motif@motif))
    motif@bkg <- c(motif@bkg, new.bkg)

  }

  validObject_universalmotif(motif)

  motif

})

# setMethod("create_motif", signature(input = "character"),
#           definition = function(input, alphabet, type, name, pseudocount,
#                                 bkg, nsites, altname, family, organism,
#                                 bkgsites, strand, pval, qval, eval,
#                                 extrainfo, add.multifreq) {
#
#   if (missing(alphabet)) alphabet <- "missing"
#   # if (missing(alphabet))  {
#   #   alphabet <- collapse_cpp(sort_unique_cpp(safeExplode(collapse_cpp(input))))
#   # }
#   else if (alphabet == "custom")
#     stop(wmsg("`alphabet = 'custom'` is no longer acceptable; please provide ",
#               "the actual letters"))
#
#   consensus <- input
#   consensus.all <- consensus
#
#   if (length(consensus) > 1) {
#     consensus <- collapse_cpp(consensus)
#     if (length(unique(nchar(input))) != 1)
#       stop("all sequences must have the same number of characters")
#   }
#
#   if (nchar(consensus) < 1) stop("sequence must be at least one character")
#
#   consensus <- sort_unique_cpp(safeExplode(consensus))
#
#   if (alphabet %in% c("DNA", "RNA") && length(consensus.all) == 1) {
#     motif <- vapply(consensus, consensus_to_ppmC, numeric(4))
#   } else if (alphabet == "AA" && length(consensus.all) == 1) {
#     motif <- vapply(consensus, consensus_to_ppmAAC, numeric(20))
#   } else if (!missing(alphabet)) {
#     motif <- consensusMatrix(collapse_cpp(consensus))
#   }
#
#   if (!alphabet %in% c("DNA", "RNA", "AA", "missing") &&
#       length(consensus.all) == 1) {
#     alph.deparsed <- sort_unique_cpp(safeExplode(alphabet))
#     if (any(!consensus %in% alph.deparsed)) {
#       stop("consensus string does not match provided alphabet")
#     }
#     motif2 <- vector("list", length(alph.deparsed))
#     mot_len <- length(consensus)
#     for (i in alph.deparsed) {
#       motif2[[i]] <- motif[rownames(motif) == i, ]
#       if (length(motif2[[i]]) == 0) motif2[[i]] <- rep(0, mot_len)
#     }
#     motif <- matrix(unlist(motif2), ncol = mot_len, byrow = TRUE)
#   }
#
#   if (alphabet == "missing") {
#     if (any(toupper(consensus) != consensus) &&
#         any(tolower(consensus) != consensus)) {
#       message(wmsg(
#         "Note: detected both upper and lower case letters. These will",
#         " be treated as individual letters."
#       ))
#     }
#     if (any(consensus %in% c("E", "F", "I", "P", "Q", "X", "Z")) &&
#         !any(consensus %in% c("O", "U", letters, as.character(0:9)))) {
#       motif <- vapply(consensus, consensus_to_ppmAAC, numeric(20))
#       alphabet <- "AA"
#     } else if (any(consensus == "U") &&
#                !any(consensus %in% c("E", "F", "I", "J", "L", "O",
#                                      "P", "Q", "T", "X", "Z",
#                                      letters, as.character(0:9)))) {
#       alphabet <- "RNA"
#       motif <- vapply(consensus, consensus_to_ppmC, numeric(4))
#     } else if (any(consensus %in% DNA_ALPHABET[-c(16:18)]) &&
#                !any(consensus %in% c("E", "F", "I", "J", "L", "O",
#                                      "P", "Q", "X", "Z", "U",
#                                      letters, as.character(0:9)))) {
#       alphabet <- "DNA"
#       motif <- vapply(consensus, consensus_to_ppmC, numeric(4))
#     } else if (length(consensus.all) == 1) {
#       alphabet <- collapse_cpp(sort_unique_cpp(consensus))
#       motif <- consensusMatrix(collapse_cpp(consensus))
#     }
#   }
#
#   if (alphabet == "AA" && length(consensus.all) > 1) {
#     motif2 <- vector("list", 20)
#     mot_len <- ncol(motif)
#     for (i in AA_STANDARD2) {
#       motif2[[i]] <- motif[rownames(motif) == i, ]
#       if (length(motif2[[i]]) == 0) motif2[[i]] <- rep(0, mot_len)
#     }
#     motif <- matrix(unlist(motif2), ncol = mot_len, byrow = TRUE)
#   }
#
#   margs <- parse_args(as.list(environment()),
#                       c("bkg", "altname", "family", "organism", "bkgsites",
#                         "strand", "pval", "qval", "eval", "extrainfo"))
#
#   margs <- c(margs, list(name = name), list(pseudocount = pseudocount))
#   if (!missing(nsites)) margs <- c(margs, list(nsites = nsites))
#   else margs <- c(margs, list(nsites = length(consensus.all)))
#
#   if (length(consensus.all) > 1) {
#
#     switch(alphabet,
#       "DNA" = {
#         consensus <- DNAStringSet(lapply(consensus.all, DNAString))
#       },
#       "RNA" = {
#         consensus <- RNAStringSet(lapply(consensus.all, RNAString))
#       },
#       "AA" = {
#         consensus <- AAStringSet(lapply(consensus.all, AAString))
#       },
#       {
#         if (alphabet != "missing") {
#           consensus <- BStringSet(lapply(consensus.all, BString))
#           alph.deparsed <- sort_unique_cpp(safeExplode(alphabet))
#           if (any(!rownames(consensusMatrix(consensus)) %in%
#                   alph.deparsed))
#             stop("consensus string does not match provided alphabet")
#         } else {
#           alphabet <- collapse_cpp(consensus)
#         }
#       }
#     )
#
#     if (!missing(type)) margs <- c(margs, list(type = type))
#     motif <- do.call(create_motif,
#                      c(list(input = consensus), margs,
#                        list(alphabet = alphabet)))
#
#     validObject_universalmotif(motif)
#     return(motif)
#
#   }
#
#   motif <- apply(motif, 2, pcm_to_ppmC, pseudocount = 0)
#   if (nchar(alphabet) == 1) stop("alphabet must be longer than 1 character")
#   motif <- do.call(universalmotif_cpp, c(list(motif = motif),
#                                      list(alphabet = alphabet),
#                                      list(type = "PPM"), margs))
#
#   if (length(consensus.all) == 1 && missing(nsites) &&
#       alphabet %in% c("DNA", "RNA")) {
#
#     input.split <- safeExplode(input)
#     if ("N" %in% input.split) {
#       motif@nsites <- 4
#       if (any(c("H", "B", "V", "D") %in% input.split))
#         motif@nsites <- 12
#     } else if (any(c("H", "B", "V", "D") %in% input.split)) {
#       motif@nsites <- 3
#       if (any(c("M", "R", "W", "S", "Y", "K") %in% input.split))
#         motif@nsites <- 6
#     } else if (any(c("M", "R", "W", "S", "Y", "K") %in% input.split))
#       motif@nsites <- 2
#
#   } else if (length(consensus.all) == 1 && missing(nsites) &&
#              alphabet == "AA") {
#
#     input.split <- safeExplode(input)
#     if ("X" %in% input.split) motif@nsites <- 20
#     else if (any(c("B", "Z", "J") %in% input.split)) motif@nsites <- 2
#
#   }
#
#   if (type == "PPM") motif <- convert_type_internal(motif, "PCM")
#   motif <- convert_type_internal(motif, type = type)
#
#   if (!missing(add.multifreq) && length(input) > 1) {
#
#     for (i in add.multifreq) {
#       motif@multifreq[[as.character(i)]] <- add_multi_cpp(input, i, rownames(motif@motif))
#     }
#
#     new.bkg <- get_bkg(BStringSet(input), k = add.multifreq,
#                        pseudocount = pseudocount, list.out = FALSE,
#                        alphabet = rownames(motif@motif))
#     motif@bkg <- c(motif@bkg, new.bkg)
#
#   }
#
#   validObject_universalmotif(motif)
#   motif
#
# })

#' @describeIn create_motif Create motif from a matrix.
#' @export
setMethod("create_motif", signature(input = "matrix"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites,
                                altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

  if (anyNA(input))
    stop("matrix cannot have NA values")

  matrix <- input

  if (!is.null(rownames(matrix)))
    matrix <- matrix[order(rownames(matrix)), ]

  if (!missing(alphabet) && alphabet == "custom")
    stop(wmsg("`alphabet = 'custom'` is no longer acceptable; please provide ",
              "the actual letters"))

  if (!missing(alphabet) &&
      !alphabet %in% c("DNA", "RNA", "AA")) {

    alph.deparsed <- sort_unique_cpp(safeExplode(alphabet))

    if (any(!rownames(matrix) %in% alph.deparsed)) {
      stop("rownames do not match provided alphabet")
    }

    if (length(alph.deparsed) != nrow(matrix)) {
      stop("alphabet length does not match number of rows")
    }

  } else if (is.null(rownames(matrix)) && missing(alphabet))
    stop("Please provide the 'alphabet' arg or a matrix with rownames")
  else if (all(rownames(matrix) %in% DNA_BASES) &&
           missing(alphabet) && nrow(matrix) == 4)
    alphabet  <- "DNA"
  else if (all(rownames(matrix) %in% RNA_BASES) &&
           missing(alphabet) && nrow(matrix) == 4)
    alphabet <- "RNA"
  else if (nrow(matrix) == 20 && missing(alphabet))
    alphabet <- "AA"
  else if (all(rownames(matrix) %in% AA_STANDARD2) &&
           missing(alphabet) && nrow(matrix) == 20)
    alphabet <- "AA"
  else if (!is.null(rownames(matrix)))
    alphabet <- collapse_cpp(rownames(matrix))
  else if (nrow(matrix) == 4 && missing(alphabet))
    alphabet <- "DNA"
  else if (missing(alphabet))
    alphabet <- collapse_cpp(rownames(matrix))

  if (alphabet %in% c("DNA", "RNA")) {
    if (nrow(matrix) != 4) stop("incorrect number of rows")
  } else if (alphabet == "AA") {
    if (nrow(matrix) != 20) stop("incorrect number of rows")
  }

  margs <- parse_args(as.list(environment()),
                      c("bkg", "nsites", "altname", "family", "organism",
                        "bkgsites", "strand", "pval", "qval", "eval",
                        "extrainfo"))
  margs <- c(margs, list(name = name), list(pseudocount = pseudocount))

  motif <- matrix

  motif <- do.call(universalmotif_cpp, c(list(motif = motif), margs,
                                         list(alphabet = alphabet)))

  if (missing(nsites)) {
    nsites <- sum(input[, 1])
    if (nsites == round(nsites) && nsites != 1 && abs(nsites) != Inf)
      motif@nsites <- nsites
  }

  motif <- convert_type_internal(motif, type = type)

  validObject_universalmotif(motif)
  motif

})

#' @describeIn create_motif Create motif from a \code{\link{DNAStringSet}}.
#' @export
setMethod("create_motif", signature(input = "DNAStringSet"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

  sequences <- input

  if (length(unique(width(sequences))) != 1)
    stop("all sequences must be the same width")

  margs <- parse_args(as.list(environment()),
                      c("nsites", "altname", "family", "organism", "bkgsites",
                        "strand", "pval", "qval", "eval", "extrainfo"))
  margs <- c(margs, list(name = name), list(pseudocount = pseudocount))

  sequences <- consensusMatrix(sequences, baseOnly = TRUE)
  if (sum(sequences[5, ]) > 0) stop("only ACGT are accepted for DNA")
  motif <- apply(sequences[1:4, ], 2, pcm_to_ppmC, pseudocount = 0)

  if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))

  motif <- do.call(universalmotif_cpp,
                   c(list(motif = motif), list(type = "PPM"), margs,
                     list(alphabet = "DNA")))

  if (missing(nsites))  motif@nsites <- length(input)
  motif <- convert_type_internal(motif, type = type)

  if (length(input) > 1 && !missing(add.multifreq)) {
    for (i in add.multifreq) {
      motif@multifreq[[as.character(i)]] <- add_multi_cpp(as.character(input), i, DNA_BASES)
    }
    new.bkg <- get_bkg(DNAStringSet(input), k = add.multifreq, list.out = FALSE,
                       pseudocount = pseudocount, alphabet = DNA_BASES)
    motif@bkg <- c(motif@bkg, new.bkg)
  }

  validObject_universalmotif(motif)
  motif

})

#' @describeIn create_motif Create motif from a \code{\link{RNAStringSet}}.
#' @export
setMethod("create_motif", signature(input = "RNAStringSet"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

  sequences <- input

  if (length(unique(width(sequences))) != 1)
    stop("all sequences must be the same width")

  margs <- parse_args(as.list(environment()),
                      c("nsites", "altname", "family", "organism", "bkgsites",
                        "strand", "pval", "qval", "eval", "extrainfo"))
  margs <- c(margs, list(name = name), list(pseudocount = pseudocount))

  sequences <- consensusMatrix(sequences, baseOnly = TRUE)
  if (sum(sequences[5, ]) > 0) stop("only ACGU are accepted for RNA")
  motif <- apply(sequences[1:4, ], 2, pcm_to_ppmC, pseudocount = 0)

  if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))

  motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                         list(type = "PPM"),
                                         margs,
                                         list(alphabet = "RNA")))

  if (missing(nsites))  motif@nsites <- length(input)
  motif <- convert_type_internal(motif, type = type)

  if (length(input) > 1 && !missing(add.multifreq)) {
    for (i in add.multifreq) {
      motif@multifreq[[as.character(i)]] <- add_multi_cpp(as.character(input), 
                                                          i, RNA_BASES)
    }
    new.bkg <- get_bkg(RNAStringSet(input), k = add.multifreq, list.out = FALSE,
                       pseudocount = pseudocount, alphabet = RNA_BASES)
    motif@bkg <- c(motif@bkg, new.bkg)
  }

  validObject_universalmotif(motif)
  motif

})

#' @describeIn create_motif Create motif from a \code{\link{AAStringSet}}.
#' @export
setMethod("create_motif", signature(input = "AAStringSet"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

  sequences <- input

  if (length(unique(width(sequences))) != 1)
    stop("all sequences must be the same width")

  margs <- parse_args(as.list(environment()),
                      c("nsites", "altname", "family", "organism", "bkgsites",
                        "strand", "pval", "qval", "eval", "extrainfo", "bkg"))
  margs <- c(margs, list(name = name), list(pseudocount = pseudocount))

  sequences <- consensusMatrix(sequences)
  if (any(!rownames(sequences) %in% AA_STANDARD2))
    stop("only ACDEFGHIKLMNPQRSTVWY are accepted for AA")

  motif <- vector("list", 20)
  mot_len <- ncol(sequences)
  for (i in AA_STANDARD2) {
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
      motif@multifreq[[as.character(i)]] <- add_multi_cpp(as.character(input), i, AA_STANDARD2)
    }
    new.bkg <- get_bkg(AAStringSet(input), k = add.multifreq, list.out = FALSE,
                       pseudocount = pseudocount, alphabet = AA_STANDARD2)
    motif@bkg <- c(motif@bkg, new.bkg)
  }

  if (missing(nsites)) motif@nsites <- length(input)
  motif <- convert_type_internal(motif, type = type)

  validObject_universalmotif(motif)
  motif

})

#' @describeIn create_motif Create motif from a \code{\link{BStringSet}}.
#' @export
setMethod("create_motif", signature(input = "BStringSet"),
          definition = function(input, alphabet, type, name, pseudocount,
                                bkg, nsites, altname, family, organism,
                                bkgsites, strand, pval, qval, eval,
                                extrainfo, add.multifreq) {

  sequences <- input

  if (!missing(alphabet) && alphabet == "custom")
    stop(wmsg("`alphabet = 'custom'` is no longer acceptable; please provide ",
              "the actual letters"))

  if (length(unique(width(sequences))) != 1)
    stop("all sequences must be the same width")

  margs <- parse_args(as.list(environment()),
                      c("nsites", "altname", "family", "organism", "bkgsites",
                        "strand", "pval", "qval", "eval", "extrainfo", "bkg"))
  margs <- c(margs, list(name = name), list(pseudocount = pseudocount))

  if (missing(alphabet)) {

    sequences <- consensusMatrix(sequences)
    motif <- apply(sequences, 2, pcm_to_ppmC, pseudocount = 0)
    alphabet <- collapse_cpp(rownames(sequences))

    motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                           list(type = "PPM"),
                                           margs,
                                           list(alphabet = alphabet)))

  } else {

    sequences <- consensusMatrix(sequences)
    alph.split <- sort_unique_cpp(safeExplode(alphabet))
    motif <- vector("list", length(alph.split))
    mot_len <- ncol(sequences)
    for (i in alph.split) {
      motif[[i]] <- sequences[rownames(sequences) == i, ]
      if (length(motif[[i]]) == 0) motif[[i]] <- rep(0, mot_len)
    }
    motif <- matrix(unlist(motif), ncol = mot_len, byrow = TRUE)
    motif <- apply(motif, 2, pcm_to_ppmC, pseudocount = 0)
    if (!is.null(rownames(motif)))
      motif <- motif[order(rownames(motif)), ]

    motif <- do.call(universalmotif_cpp, c(list(motif = motif),
                                           list(type = "PPM"), margs,
                                           list(alphabet = alphabet)))

  }

  if (length(input) > 1 && !missing(add.multifreq)) {
    for (i in add.multifreq) {
      motif@multifreq[[as.character(i)]] <- add_multi_cpp(as.character(input), i,
                                                          rownames(motif@motif))
    }
    new.bkg <- get_bkg(BStringSet(input), k = add.multifreq,
                       pseudocount = pseudocount, list.out = FALSE,
                       alphabet = rownames(motif@motif))
    motif@bkg <- c(motif@bkg, new.bkg)
  }

  if (missing(nsites)) motif@nsites <- length(input)
  motif <- convert_type_internal(motif, type = type)

  validObject_universalmotif(motif)
  motif

})

parse_args <- function(args, names) {

  # param check --------------------------------------------
  all_checks <- character(0)
  char_check <- check_fun_params(list(alphabet = args$alphabet,
                                      type = args$type, name = args$name,
                                      altname = args$altname,
                                      family = args$family,
                                      organism = args$organism,
                                      strand = args$strand,
                                      extrainfo = args$extrainfo),
                                 c(rep(1, 7), 0), rep(TRUE, 8), TYPE_CHAR)
  num_check <- check_fun_params(list(pseudocount = args$pseudocount,
                                     bkg = args$bkg, nsites = args$nsites,
                                     bkgsites = args$bkgsites,
                                     pval = args$pval, qval = args$qval,
                                     eval = args$eval,
                                     add.multifreq = args$add.multifreq),
                                c(1, 0, 1, 1, 1, 1, 1, 0), rep(TRUE, 8),
                                TYPE_NUM)
  all_checks <- c(all_checks, char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  to.keep <- !vapply(args, is.symbol, logical(1))
  args <- args[to.keep]
  args <- args[names(args) %in% names]

  args

}
