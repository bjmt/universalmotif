#' Create random sequences.
#'
#' Generate random sequences from any set of characters, represented as
#' \code{\link{XStringSet}} objects.
#'
#' @param alphabet `character(1)` One of `c('DNA', 'RNA', 'AA')`, or a string of
#'    characters to be used as the alphabet.
#' @param monofreqs `numeric` Alphabet frequencies to use. If missing assumes uniform
#'   frequencies. Not used if `difreq` or `trifreq` are
#'   input.
#' @param seqnum `numeric(1)` Number of sequences to generate.
#' @param seqlen `numeric(1)` Length of random sequences.
#' @param difreqs `numeric` Dinucleotide frequencies. DNA/RNA only. Must be a
#'   named numeric vector of length 16.
#' @param trifreqs `numeric` Trinucleotide frequencies. DNA/RNA only. Must be a 
#'   named numeric vector of length 64.
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [create_sequences()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for large jobs (such as
#'    `create_sequences(seqlen=100000,seqnum=100000)`). Furthermore,
#'    the behaviour of `progress = TRUE` is
#'    changed if `BP = TRUE`; the default \pkg{BiocParallel} progress bar will
#'    be shown (which unfortunately is much less informative).
#'
#' @return \code{\link{XStringSet}} The returned sequences are _unnamed_.
#'
#' @examples
#' ## create DNA sequences with slightly increased AT content:
#' sequences <- create_sequences(monofreqs = c(0.3, 0.2, 0.2, 0.3))
#' ## create custom sequences:
#' sequences.QWER <- create_sequences("QWER")
#' ## you can include non-alphabet characters are well, even spaces:
#' sequences.custom <- create_sequences("!@#$ ")
#'
#' @references
#'    \insertRef{biostrings}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [create_motif()], [shuffle_sequences()]
#' @export
create_sequences <- function(alphabet = "DNA", seqnum = 100, seqlen = 100,
                             monofreqs, difreqs, trifreqs, progress = FALSE,
                             BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(alphabet = args$alphabet),
                                 1, FALSE, "character")
  num_check <- check_fun_params(list(seqnum = args$seqnum, seqlen = args$seqlen,
                                     monofreqs = args$monofreqs,
                                     difreqs = args$difreqs,
                                     trifreqs = args$trifreqs),
                                c(1, 1, rep(0, 3)), c(FALSE, FALSE, rep(TRUE, 3)),
                                "numeric")
  logi_check <- check_fun_params(list(progress = args$progress, BP = args$BP),
                                 numeric(), logical(), "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (alphabet == "DNA" || alphabet == "RNA") {
    alph.letters <- DNA_BASES
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
    seqs <- lapply_(seq_len(seqnum),
                     function(x) create_k1(alph.letters = alph.letters,
                                           seqlen = seqlen,
                                           bkg = monofreqs),
                    BP = BP, PB = progress)
  } else if (!missing(difreqs)) {
    names(difreqs) <- gsub("U", "T", names(difreqs))
    seqs <- lapply_(seq_len(seqnum),
                     function(x) create_k2(alph.letters = DNA_BASES,
                                           seqlen = seqlen,
                                           difreq = difreqs),
                    BP = BP, PB = progress)
  } else if (!missing(trifreqs)) {
    names(trifreqs) <- gsub("U", "T", names(trifreqs))
    seqs <- lapply_(seq_len(seqnum),
                     function(x) create_k3(alph.letters = DNA_BASES,
                                           seqlen = seqlen,
                                           trifreq = trifreqs),
                    BP = BP, PB = progress)
  }

  seqs <- unlist(seqs)

  if (alphabet == "DNA") {
    seqs <- DNAStringSet(seqs)
  } else if (alphabet == "RNA") {
    seqs <- DNAStringSet(seqs)
    seqs <- RNAStringSet(seqs)
  } else if (alphabet == "AA") {
    seqs <- AAStringSet(seqs)
  } else seqs <- BStringSet(seqs)

  seqs

}

create_k1 <- function(alph.letters, seqlen, bkg) {
  seqout <- sample(alph.letters, seqlen, replace = TRUE, prob = bkg)
  seqout <- collapse_cpp(seqout)
  seqout
}

create_k2 <- function(alph.letters, seqlen, difreq) {
  probsA <- difreq[c("AA", "AC", "AG", "AT")]
  probsC <- difreq[c("CA", "CC", "CG", "CT")]
  probsG <- difreq[c("GA", "GC", "GG", "GT")]
  probsT <- difreq[c("TA", "TC", "TG", "TT")]
  ditrans <- matrix(c(probsA, probsC, probsG, probsT), nrow = 4)
  rownames(ditrans) <- alph.letters
  seqout <- rep(NA, seqlen)
  first.di <- sample(names(difreq), 1, prob = difreq)
  first.di <- strsplit(first.di, "")[[1]]
  seqout[1] <- first.di[1]
  seqout[2] <- first.di[2]
  for (i in 3:seqlen) {
    previous.nuc <- seqout[i - 1]
    curr.prob <- ditrans[previous.nuc, ]
    curr.prob[is.na(curr.prob)] <- 0.00001
    seqout[i] <- sample(alph.letters, 1, prob = curr.prob)
  }
  seqout <- collapse_cpp(seqout)
}

create_k3 <- function(alph.letters, seqlen, trifreq) {
  trifreq <- trifreq[DNA_TRI]
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
  seqout <- collapse_cpp(seqout)
}
