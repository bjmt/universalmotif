#' Create a TFFM. 
#'
#' If the sequences of all binding sites are available, a transcription
#' factor flexible model \insertCite{tffm}{universalmotif} can be generated.
#' TFFMs are based on hidden Markov models and are a more flexible
#' representation of a transcription factor binding site. The code to generate
#' TFFMs using this function is based on the python module 'TFFM'
#' \insertCite{tffm}{universalmotif}.
#'
#' @param sequences DNAStringSet (or list of). All sequences must be of the
#'                  same width.
#' @param memefile Character. File path to MEME output.
#' @param type Character. Currently only 'first' is supported.
#' @param ID Character. Motif ID.
#' @param name Character. Motif name.
#' @param strand Character. Motif strand.
#' @param family Character. Transcription factor family.
#' @param bkg Numeric. Background frequencies.
#' @param pseudocount Numeric. Pseudocount to be added when calculating TFFM.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return TFFMFirst object, TFFMDetail object, or a list of such objects.
#'
#' @details
#'    Code based on TFFM python module \insertCite{tffm}{universalmotif}.
#'    Currently the only support for this type of motif in Bioconductor is
#'    via the TFBSTools package \insertCite{tfbstools}{universalmotif}. As a
#'    result, the TFFM is stored using the TFFM class from TFBSTools.
#'
#'    If a MEME output file is used to create the TFFM, then the 'strand'
#'    and 'bkg' fields will be taken from the MEME output. If there are
#'    multiple motifs in the output file, a TFFM will be created from each
#'    and a list of TFFM objects returned.
#'
#' @seealso \code{\link{create_motif}}
#'
#' @examples
#'    library(Biostrings)
#'    sites <- readDNAStringSet(system.file("extdata", "sites.txt",
#'                                          package = "universalmotif"))
#'    TFFM <- create_tffm(sites)
#'    TFBSTools::seqLogo(TFFM)
#'
#' @references
#'    \insertRef{tffm}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
create_tffm <- function(sequences, memefile, type = "first", ID = "Unknown",
                        name = "Unknown", strand = "+", family = "Unknown",
                        bkg = c(0.25, 0.25, 0.25, 0.25),
                        pseudocount = 1, BPPARAM = bpparam()) {

  if (!missing(sequences) && !missing(memefile)) {
    stop("please only enter one of 'sequences' or 'memefile'")
  }

  if (is.list(sequences)) {
    tffms <- bplapply(sequences,
                      function(x) create_tffm(sequences = x,
                                              type = type, ID = ID,
                                              name = name, strand = strand,
                                              family = family, bkg = bkg,
                                              pseudocount = pseudocount),
                      BPPARAM = BPPARAM)
    return(tffms)
  }

  if (!missing(memefile)) {
    memefile <- read_meme(memefile, readsites = TRUE, BPPARAM = BPPARAM)
    alph <- memefile[[1]][[1]]["alphabet"]
    if (alph != "DNA") stop("'create_tffm' can only handle DNA motifs")
    strand <- memefile[[1]][[1]]["strand"]
    bkg <- memefile[[1]][[1]]["bkg"]
    memefile <- memefile[[2]]
    if (length(memefile) == 0) stop("no sites found in MEME file")
    tffms <-  bplapply(memefile,
                       function(x) create_tffm(sequences = x, type = type, ID = ID,
                                               name = name, strand = strand,
                                               family = family, bkg = bkg,
                                               pseudocount = pseudocount),
                        BPPARAM = BPPARAM)
    if (length(tffms) == 1) return(tffms[[1]]) else return(tffms)
  }

  if (class(sequences) != "DNAStringSet") {
    stop("'sequences' must a DNAStringSet")
  }
  if (!is.character(type) || length(type) != 1) {
    stop("incorrect 'type'")
  }
  if (!is.character(ID) || length(ID) != 1) {
    stop("incorrect 'ID'")
  }
  if (!is.character(name) || length(name) != 1) {
    stop("incorrect 'name'")
  }
  if (!is.character(strand) || !length(strand) %in% 1:2) {
    stop("incorrect 'strand'")
  }
  if (!is.character(family) || length(family) != 1) {
    stop("incorrect 'family'")
  }
  if (!is.numeric(bkg) || length(bkg) != 4) {
    stop("'bkg' must be a numeric of length 4")
  }
  if (!is.numeric(pseudocount) || length(pseudocount) != 1) {
    stop("'pseudocount' must be a numeric of length 1")
  }

  pc1 <- pseudocount
  pc4 <- pseudocount * 4

  names(bkg) <- DNA_BASES

  motif.mat <- as.matrix(sequences)
  mot.len <- ncol(motif.mat)
  motif.pfm <- consensusMatrix(sequences, baseOnly = TRUE)
 
  num.hits <- nrow(motif.mat)
  first.letters <- rowSums(consensusMatrix(sequences))[1:4]
  # num.residues <- sum(first.letters)  # the python module uses num.hits * 100;
  num.residues <- num.hits * 100        # though apparently MEME has this info
  first.letters <- rep(num.residues / 4, 4)  # want bkg frequencies, not site..

  if (type == "first") {

    emissions <- list(vapply(first.letters, function(x) (x + pc1) / (num.residues + pc4),
                             numeric(1)))

    for (position in seq_len(mot.len)) {
      frequencies <- vector("numeric", 1)
      for (letter in 1:4) {
        freq <- (motif.pfm[letter, position] + pc1) / (num.hits + pc4)
        frequencies[letter] <- freq
      }
      emissions <- c(emissions, list(rep(frequencies, 4)))
    }

    transition.bkg <- list(c(0, 1, rep(0, mot.len)))

    bkg2bkg <- 1 - num.hits / (num.hits * num.residues)
    bkg2fg <- 1 - bkg2bkg

    transition.core <- list(c(0, bkg2bkg, bkg2fg, rep(0, mot.len - 1)))
    for (position in seq_len(mot.len - 1)) {
      transition.pos <- c(rep(0, position + 2), 1, rep(0, mot.len - position - 1))
      transition.core <- c(transition.core, list(transition.pos))
    }
    transition.final <- list(c(0, 1, rep(0, mot.len)))

    transition <- c(transition.bkg, transition.core, transition.final)
    transition <- matrix(unlist(transition), nrow = length(transition), byrow = TRUE)

    dimnames(transition) <- list(0:(nrow(transition) - 1),
                                 0:(ncol(transition) - 1))

    # initial.probabilities <- c(2, rep(0, num.hits - 1))

    TFFMFirst(emission = emissions, transition = transition[, -1],
              profileMatrix = matrix(motif.pfm[1:4, ],
                                     nrow = 4, byrow = TRUE,
                                     dimnames = list(c("A", "C", "G", "T"))),
              type = "First", name = name, ID = ID, matrixClass = family,
              bg = bkg)

  } else if (type == "detailed") {

    warning("type 'detailed' currently does not work as intended")
  
    emissions <- rep(list(c(1,0,0,0), c(0,1,0,0), c(0,0,1,0), c(0,0,0,1)), mot.len + 1)

    bkg2bkg <- 1 - num.hits / (num.hits * num.residues)
    bkg2fg <- 1 - bkg2bkg
    bkg2bkg <- bkg2bkg / 4

    transition.i <- vector("numeric", 4)
    for (letter in 1:4) {
      freq <- (motif.pfm[letter, 1] + pc1) / (num.hits + pc4)
      transition.i[letter] <- freq * bkg2fg
    }

    transitions <- c(rep(bkg2bkg, 4), transition.i, rep(0, 4 * (mot.len - 1)))
    transitions <- list(transitions, transitions, transitions, transitions)

    pfm <- vector()
    for (letter in 1:4) {
      for (position in seq_len(mot.len)) {
        pfm <- c(pfm, (motif.pfm[letter, position] + pc1) /
                      (num.hits + pc4))
      }
    }

    for (position in seq_len(mot.len * 4)) {
      transitions <- c(transitions, list(rep(0, 4 * (mot.len + 1))))
    }

    for (position in seq_len(mot.len - 1)) {
      for (line in (4 * (position - 1) + 1):(4 * (position - 1) + 5 - 1)) {
        for (column in (4 * position + 1):(4 * position + 5 - 1)) {
          index <- (column - (4 * position + 1)) * mot.len + position
          transitions[[line + 3]][column + 3] <- pfm[index]
        }
      }
    }

    for (index in 1:4) {
      state <- mot.len * 4 + (index - 1)
      transitions[[state]][1] <- 0.25
      transitions[[state]][2] <- 0.25
      transitions[[state]][3] <- 0.25
      transitions[[state]][4] <- 0.25
    }

    transitions <- matrix(unlist(transitions), nrow = length(transitions),
                          byrow = TRUE) + 0.00000001  # somethings wrong,
                                                      # shouldn't need this(?)
    dimnames(transitions) <- list(0:(nrow(transitions) - 1),
                                  0:(ncol(transitions) - 1))

    TFFMDetail(emission = emissions, transition = transitions,
              profileMatrix = matrix(motif.pfm[1:4, ],
                                     nrow = 4, byrow = TRUE,
                                     dimnames = list(c("A", "C", "G", "T"))),
              type = "Detailed", name = name, ID = ID, matrixClass = family,
              bg = c(A = bkg[1], C = bkg[2], G = bkg[3], T = bkg[4]))
  
  } else stop("unknown type")

}

# possible workflow:

#   create_tffm2 = create pHMM with aphid
#     - option to create from universalmotif obj? would need to train..
#     - create from meme output file?
#   train_tffm2 = train pHMM with aphid using BaumWelch algo
#   scan_tffm2 = scan pHMM with HMMER3 (writePHMM)
