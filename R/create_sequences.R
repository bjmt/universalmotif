#' Create random background sequences.
#'
#' @param alphabet Character. One of 'DNA', 'RNA', 'AA', or a string of letters
#'                 to be used as the alphabet.
#' @param bkg Numeric. Alphabet frequencies to use. If missing assumes uniform
#'            frequencies.
#' @param numseqs Numeric. Number of sequences to generate.
#' @param seqlen Numeric. Length of random sequences.
#' @param difreq Numeric. Dinucleotide frequencies. DNA only. Must be a
#'               named numeric vector of length 16.
#'
#' @return XStringSet object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
create_sequences <- function(alphabet = "DNA", bkg, numseqs = 100, seqlen = 100,
                             difreq) {

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
  if (missing(difreq)) {
    for (i in seq_len(numseqs)) {
      seqs[[i]] <- sample(alph.letters, seqlen, replace = TRUE, prob = bkg)
      seqs[[i]] <- paste(seqs[[i]], collapse = "")
    }
  } else {
    if (!alphabet == "DNA") {
      stop("if 'difreq' is provided, alphabet must be 'DNA'")
    }
    if (length(difreq) != 16) stop("'difreq' must be length 16")
    dinucs <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC",
                "GG", "GT", "TA", "TC", "TG", "TT")
    if (!all(names(difreq) %in% dinucs)) {
      stop("dinucleotide frequncies must be provided for ", dinucs)
    }
    probsA <- difreq[c("AA", "AC", "AG", "AT")]
    probsC <- difreq[c("CA", "CC", "CG", "CT")]
    probsG <- difreq[c("GA", "GC", "GG", "GT")]
    probsT <- difreq[c("TA", "TC", "TG", "TT")]
    ditrans <- matrix(c(probsA, probsC, probsG, probsT), nrow = 4)
    rownames(ditrans) <- alph.letters
    seqs.out <- vector("list", numseqs)
    for (i in seq_len(numseqs)) {
      seqs.out[[i]] <- rep(NA, seqlen)
      seqs.out[[i]][1] <- sample(alph.letters, 1, prob = bkg)
      for (j in 2:seqlen) {
        previous.nuc <- seqs.out[[i]][j - 1]
        curr.prob <- ditrans[previous.nuc, ]
        seqs.out[[i]][j] <- sample(alph.letters, 1, prob = curr.prob)
      }
    }
    seqs <- lapply(seqs.out, function(x) paste(x, collapse = ""))
  }

  seqs <- unlist(seqs)

  seq.names <- seq_len(numseqs)

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
