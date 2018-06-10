#' Shuffle input sequences or generate random sequences.
#'
#' @param sequences XStringSet objects.
#' @param inter Logical. Shuffle letters between sequences? If \code{FALSE},
#'              then the alphabet frequencies for each sequence are maintained.
#'              Ignored if \code{sequences} is missing.
#' @param alphabet Character. Ignored if \code{sequences} is given. If not one
#'                 of 'DNA', 'RNA', or 'AA', then a sequence of letters
#'                 representing the alphabet to use (e.g. "QWERTY").
#' @param bkg Numeric. Ignored if \code{sequences} is given; otherwise generate
#'            random sequences with these alphabet frequencies. If missing,
#'            will generate sequences with even background frequencies.
#' @param numseqs Numeric. Ignored if \code{sequences} is given; otherwise
#'                the number of random sequences to generate.
#' @param seqlen Numeric. Ignored if \code{sequences} is given;otherwise the
#'               length of the random sequences.
#'
#' @return XStringSet object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
shuffle_sequences <- function(sequences, inter = FALSE, alphabet, bkg,
                              numseqs = 100, seqlen = 100) {

  if (!missing(sequences)) {

    alph <- sequences@elementType
    if (alph == "DNAString") {
      alph <- "DNA"
    } else if (alph == "RNAString") {
      alph <- "RNA"
    } else if (alph == "AAString") {
      alph <- "AA"
    } else alph <- "B"

    seq.names <- names(sequences)
    if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
    seqs <- as.character(sequences)
    widths <- width(sequences)
    seqs <- lapply(seqs, function(x) strsplit(x, "")[[1]])

    if (!inter) {
      seqs <- mapply(function(x, y) sample(x, y), seqs, widths,
                     SIMPLIFY = FALSE)
      seqs <- lapply(seqs, function(x) paste(x, collapse = ""))
    } else if (inter) {
      seqs.all <- unlist(seqs)
      seqs.all <- sample(seqs.all, length(seqs.all))
      seqs.i <- mapply(function(x, y) as.factor(rep(x, y)),
                       seq_len(length(sequences)), widths, SIMPLIFY = FALSE)
      seqs.i <- unlist(seqs.i)
      seqs <- split(seqs.all, seqs.i)
      seqs <- lapply(seqs, function(x) paste(x, collapse = ""))
    } else stop("'inter' must be TRUE or FALSE")

    seqs <- unlist(seqs)

  } else {

    if (missing(alphabet)) stop("one of 'sequences' or 'alphabet' must be input")

    if (alphabet == "DNA") {
      alph.letters <- DNA_BASES
    } else if (alphabet == "RNA") {
      alph.letters <- RNA_BASES
    } else if (alphabet == "AA") {
      alph.letters <- AA_STANDARD
    } else {
      alph.letters <- strsplit(alphabet, "")[[1]]
    }

    if (!missing(bkg) && length(bkg) != length(alph.letters)) {
      stop("'bkg' and 'alphabet' must be of the same length")
    }
    if (missing(bkg)) bkg <- rep(1 / length(alph.letters), length(alph.letters))
    
    seqs <- vector("list", numseqs)
    for (i in seq_len(numseqs)) {
      seqs[[i]] <- sample(alph.letters, seqlen, replace = TRUE, prob = bkg)
      seqs[[i]] <- paste(seqs[[i]], collapse = "")
    }

    seqs <- unlist(seqs)

    seq.names <- seq_len(numseqs)
    
  }

  if (alphabet == "DNA") {
    seqs <- DNAStringSet(seqs)
  } else if (alphabet == "RNA") {
    seqs <- RNAStringSet(seqs)
  } else if (alphabet == "AA") {
    seqs <- AAStringSet(seqs)
  } else seqs <- BStringSet(seqs)

  names(seqs) <- seq.names

  seqs

}
