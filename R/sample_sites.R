#' Generate binding sites from a motif.
#'
#' Given probabilities for a sequence as represented by a motif, generate
#' random sequences with the same length as the motif.
#'
#' @param motif See [convert_motifs()] for acceptable formats.
#' @param n `numeric(1)` Number of sites to generate.
#' @param use.freq `numeric(1)` If one, use regular motif matrix. Otherwise,
#'    use respective `multifreq` matrix.
#'
#' @return \code{\link{XStringSet}} object.
#'
#' @examples
#' motif <- create_motif()
#' sites <- sample_sites(motif)
#'
#' @seealso [create_sequences()], [create_motif()], [add_multifreq()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
sample_sites <- function(motif, n = 100, use.freq = 1) {

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list(n = args$n), 1, FALSE, TYPE_NUM)
  all_checks <- c(num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motif <- convert_motifs(motif)
  motif <- convert_type_internal(motif, "PPM")

  if (use.freq == 1) {
    mot.mat <- motif@motif
  } else {
    mot.mat <- motif@multifreq[[as.character(use.freq)]]
  }

  alph <- rownames(motif@motif)

  if (use.freq == 1) {
    sites <- vapply(seq_len(n), function(x) sample_k1(mot.mat, alph, n),
                    character(1))
  } else {
    sites <- sample_motif(use.freq, mot.mat, alph, n)
  }

  alph <- motif@alphabet

  sites <- switch(alph, "DNA" = DNAStringSet(sites), "RNA" = RNAStringSet(sites),
                  "AA" = AAStringSet(sites), BStringSet(sites))

  sites

}

sample_k1 <- function(mat, alph, n) {
  collapse_cpp(apply(mat, 2, function(x) sample(alph, 1, prob = x)))
}

sample_motif <- function(k, kmat, alph, n) {

  seq.len <- ncol(kmat) + k - 1
  alph.len <- length(alph)

  k.lets <- rownames(kmat)
  k.lets.m1 <- get_klets(alph, k - 1)
  k.lets.m1 <- rep(k.lets.m1, each = alph.len)

  seqout <- vapply(seq_len(n),
                   function(x) sample_single(k.lets, kmat, seq.len, k, alph,
                                             k.lets.m1),
                   character(1))

  seqout

}

sample_single <- function(k.lets, kmat, seq.len, k, alph, k.lets.m1) {

  first.k <- sample(k.lets, 1, prob = kmat[, 1])
  seqout <- character(seq.len)
  seqout[seq_len(k)] <- safeExplode(first.k)

  for (i in seq_len(seq.len - k)) {

    j <- i + k

    prev.k <- collapse_cpp(seqout[seq(j - k + 1, j - 1)])
    seqout[j] <- sample(alph, 1, prob = kmat[k.lets.m1 == prev.k, i + 1])

  }

  collapse_cpp(seqout)

}

# > 1-letter steps (k = 3)
#
#        A C A | | | | | | |
#  i=1 . - C A A | | | | | | . l=2 . r=4
#  i=2 . --- A A T | | | | | . l=3 . r=5
#  i=3 . ----- A T G | | | | . l=4 . r=6
#  i=4 . ------- T G C | | | . l=5 . r=7
#  i=5 . --------- G C C | | . l=6 . r=8
#  i=6 . ----------- C C C | . l=7 . r=9
#  i=7 . ------------- C C G . l=8 . r=10
#        ===================
#        A C A A T G C C C G
#        | | | | | | | | | |
#        1 2 3 4 5 6 7 8 9 10
