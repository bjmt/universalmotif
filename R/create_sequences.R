#' Create random background sequences.
#'
#' @param alphabet Character. One of 'DNA', 'RNA', 'AA', or a string of letters
#'                 to be used as the alphabet.
#' @param monofreqs Numeric. Alphabet frequencies to use. If missing assumes uniform
#'            frequencies. Not used if \code{difreq} or \code{trifreq} are
#'            input.
#' @param seqnum Numeric. Number of sequences to generate.
#' @param seqlen Numeric. Length of random sequences.
#' @param difreqs Numeric. Dinucleotide frequencies. DNA/RNA only. Must be a
#'               named numeric vector of length 16.
#' @param trifreqs Numeric. Trinucleotide frequencies. DNA/RNA only. Must be a 
#'                named numeric vector of length 64.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return XStringSet object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
create_sequences <- function(alphabet = "DNA", seqnum = 100, seqlen = 100,
                             monofreqs, difreqs, trifreqs,
                             BPPARAM = bpparam()) {

  if (alphabet == "DNA") {
    alph.letters <- DNA_BASES
  } else if (alphabet == "RNA") {
    alph.letters <- RNA_BASES
  } else if (alphabet == "AA") {
    alph.letters <- AA_STANDARD
    if (!missing(difreqs) || !missing(trifreqs)) {
      stop("'difreqs' and 'trifreqs' can only be used for 'DNA' and 'RNA'")
    }
  } else {
    alph.letters <- strsplit(alphabet, "")[[1]]
    if (!missing(difreqs) || !missing(trifreqs)) {
      stop("'difreqs' and 'trifreqs' can only be used for 'DNA' and 'RNA'")
    }
  }

  if (!missing(monofreqs) && length(monofreqs) != length(alph.letters)) {
    stop("'monofreqs' and 'alphabet' must have the same number of letters")
  }
  if (!missing(difreqs) && length(difreqs) != 16) {
    stop("'difreqs' must be of length 16")
  }
  if (!missing(trifreqs) && length(trifreqs) != 64) {
    stop("'trifreqs' must be of length 64")
  }

  if (missing(monofreqs) && missing(difreqs) && missing(trifreqs)) {
    monofreqs <- rep(1 / length(alph.letters), length(alph.letters))
  }

  seqs <- vector("list", seqnum)
  if (!missing(monofreqs)) {
    seqs <- bplapply(seq_len(seqnum),
                     function(x) create_k1(alph.letters = alph.letters,
                                           seqlen = seqlen,
                                           bkg = monofreqs),
                     BPPARAM = BPPARAM)
  } else if (!missing(difreqs)) {
    difreqs <- gsub("U", "T", difreqs)
    seqs <- bplapply(seq_len(seqnum),
                     function(x) create_k2(alph.letters = alph.letters,
                                           seqlen = seqlen,
                                           difreq = difreqs),
                     BPPARAM = BPPARAM)
  } else if (!missing(trifreqs)) {
    trifreqs <- gsub("U", "T", trifreqs)
    seqs <- bplapply(seq_len(seqnum),
                     function(x) create_k3(alph.letters = alph.letters,
                                           seqlen = seqlen,
                                           trifreq = trifreqs),
                     BPPARAM = BPPARAM)
  }

  seqs <- unlist(seqs)

  if (alphabet == "DNA") {
    seqs <- DNAStringSet(seqs)
  } else if (alphabet == "RNA") {
    seqs <- RNAStringSet(seqs)
  } else if (alphabet == "AA") {
    seqs <- AAStringSet(seqs)
  } else seqs <- BStringSet(seqs)

  seqs

}

create_k1 <- function(alph.letters, seqlen, bkg) {
  seqout <- sample(alph.letters, seqlen, replace = TRUE, prob = bkg)
  seqout <- paste(seqout, collapse = "")
}

create_k2 <- function(alph.letters, seqlen, difreq) {
  probsA <- difreq[c("AA", "AC", "AG", "AT")]
  probsC <- difreq[c("CA", "CC", "CG", "CT")]
  probsG <- difreq[c("GA", "GC", "GG", "GT")]
  probsT <- difreq[c("TA", "TC", "TG", "TT")]
  ditrans <- matrix(c(probsA, probsC, probsG, probsT), nrow = 4)
  rownames(ditrans) <- alph.letters
  seqout <- rep(NA, seqlen)
  first.di <- sample(names(difreq), 1, difreq)
  first.di <- strsplit(first.di, "")[[1]]
  seqout[1] <- first.di[1]
  seqout[2] <- first.di[2]
  for (i in 3:seqlen) {
    previous.nuc <- seqout[i - 1]
    curr.prob <- ditrans[previous.nuc, ]
    curr.prob[is.na(curr.prob)] <- 0.00001
    seqout[i] <- sample(alph.letters, 1, prob = curr.prob)
  }
  seqout <- paste(seqout, collapse = "")
}

create_k3 <- function(alph.letters, seqlen, trifreq) {
  trinucs <- names(trinucleotideFrequency(DNAString("A")))
  trifreq <- trifreq[trinucs]
  tritrans <- matrix(trifreq, nrow = 16, byrow = 16)
  rownames(tritrans) <- DNA_DI
  seqout <- rep(NA, seqlen)
  first.tri <- sample(names(trifreq), 1, prob = trifreq)
  first.tri <- strsplit(first.tri, "")[[1]]
  seqout[1] <- first.tri[1]
  seqout[2] <- first.tri[2]
  seqout[3] <- first.tri[3]
  for (i in 4:seqlen) {
    previous.nuc1 <- seqout[i - 1]
    previous.nuc2 <- seqout[i - 2]
    previous.nuc <- paste0(previous.nuc2, previous.nuc1)
    curr.prob <- tritrans[previous.nuc, ]
    curr.prob[is.na(curr.prob)] <- 0.00001
    seqout[i] <- sample(alph.letters, 1, prob = curr.prob)
  }
  seqout <- paste(seqout, collapse = "")
}
