#' Create random background sequences.
#'
#' @param alphabet Character. One of 'DNA', 'RNA', 'AA', or a string of letters
#'                 to be used as the alphabet.
#' @param bkg Numeric. Alphabet frequencies to use. If missing assumes uniform
#'            frequencies. Not used if \code{difreq} or \code{trifreq} are
#'            input.
#' @param numseqs Numeric. Number of sequences to generate.
#' @param seqlen Numeric. Length of random sequences.
#' @param difreq Numeric. Dinucleotide frequencies. DNA/RNA only. Must be a
#'               named numeric vector of length 16.
#' @param trifreq Numeric. Trinucleotide frequencies. DNA/RNA only. Must be a 
#'                named numeric vector of length 64.
#'
#' @return XStringSet object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
create_sequences <- function(alphabet = "DNA", bkg, numseqs = 100, seqlen = 100,
                             difreq, trifreq) {

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
  if (missing(difreq) && missing(trifreq)) {
    for (i in seq_len(numseqs)) {
      seqs[[i]] <- sample(alph.letters, seqlen, replace = TRUE, prob = bkg)
      seqs[[i]] <- paste(seqs[[i]], collapse = "")
    }
  } else if (!missing(difreq)) {
    if (!alphabet %in% c("DNA", "RNA")) {
      stop("if 'difreq' is provided, alphabet must be 'DNA' or 'RNA'")
    }
    difreq <- gsub("U", "T", difreq)
    if (length(difreq) != 16) stop("'difreq' must be length 16")
    dinucs <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC",
                "GG", "GT", "TA", "TC", "TG", "TT")
    if (!all(names(difreq) %in% dinucs)) {
      if (alphabet == "DNA") {
        stop("dinucleotide frequncies must be provided for ", dinucs)
      } else {
        dinucs <- gsub("T", "U", dinucs)
        stop("dinucleotide frequncies must be provided for ", dinucs)
      }
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
      first.di <- sample(names(difreq), 1, difreq)
      first.di <- strsplit(first.di, "")[[1]]
      seqs.out[[i]][1] <- first.di[1]
      seqs.out[[i]][2] <- first.di[2]
      for (j in 3:seqlen) {
        previous.nuc <- seqs.out[[i]][j - 1]
        curr.prob <- ditrans[previous.nuc, ]
        curr.prob[is.na(curr.prob)] <- 0.01
        seqs.out[[i]][j] <- sample(alph.letters, 1, prob = curr.prob)
      }
    }
    seqs <- lapply(seqs.out, function(x) paste(x, collapse = ""))
  } else if (!missing(trifreq)) {
    if (!alphabet %in% c("DNA", "RNA")) {
      stop("if 'trifreq' is provided, alphabet must be 'DNA' or 'RNA'")
    }
    if (length(trifreq) != 64) stop("'trifreq' must be length 64")
    trinucs <- names(trinucleotideFrequency(DNAString("A")))
    trifreq <- gsub("U", "T", trifreq)
    if (!all(names(trifreq) %in% trinucs)) {
      if (alphabet == "DNA") {
        stop("trinucleotide frequencies must be provided for ", trinucs)
      } else {
        trinucs <- gsub("T", "U", trinucs)
        stop("trinucleotide frequencies must be provided for ", trinucs)
      }
    }
    trifreq <- trifreq[trinucs]
    tritrans <- matrix(trifreq, nrow = 16, byrow = TRUE)
    rownames(tritrans) <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                            "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
    seqs.out <- vector("list", numseqs)
    for (i in seq_len(numseqs)) {
      seqs.out[[i]] <- rep(NA, seqlen)
      first.tri <- sample(names(trifreq), 1, prob = trifreq)
      first.tri <- strsplit(first.tri, "")[[1]]
      seqs.out[[i]][1] <- first.tri[1]
      seqs.out[[i]][2] <- first.tri[2]
      seqs.out[[i]][3] <- first.tri[3]
      for (j in 4:seqlen) {
        previous.nuc1 <- seqs.out[[i]][j - 1]
        previous.nuc2 <- seqs.out[[i]][j - 2]
        previous.nuc <- paste0(previous.nuc2, previous.nuc1)
        curr.prob <- tritrans[previous.nuc, ]
        curr.prob[is.na(curr.prob)] <- 0
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
