#' Add multi-letter information to a motif.
#'
#' If the original sequences are available for a particular motif, then they
#' can be used to generate higher-order PPM matrices.
#'
#' @param motif See [convert_motifs()] for acceptable formats. If the
#'    motif is not a [universalmotif-class] motif, then it will be
#'    converted.
#' @param sequences \code{\link{XStringSet}} The alphabet must match
#'    that of the motif. If
#'    these sequences are all the same length as the motif, then they are all
#'    used to generate the multi-freq matrices. Otherwise
#'    [scan_sequences()] is first run to find the right sequence.
#' @param add.k `numeric(1)` The k-let lengths to add.
#' @param threshold `numeric(1)` Between 0 and 1. See [scan_sequences()].
#' @param threshold.type `character(1)` One of `c('logodds', 'pvalue')`.
#'    See [scan_sequences()].
#' @param RC `logical(1)` Check the reverse complement of a DNA sequence.
#'    See [scan_sequences()].
#' @param motifs.perseq `numeric(1)` If [scan_sequences()] is run,
#'    then this indicates how many hits from each sequence is to be used.
#'
#' @details
#'    At each position in the motif, then the probability of each k-let 
#'    covering from the initial position to ncol - 1 is calculated. Only
#'    positions within the motif are considered; this means that the
#'    final k-let probability matrix will have ncol - 1 fewer columns.
#'    Calculating k-let probabilities for the missing columns would be
#'    trivial however, as you would only need the background frequencies.
#'    Since these would not be useful for [scan_sequences()]
#'    though, they are not calculated.
#'
#'    Currently [add_multifreq()] does not try to stay faithful to the default
#'    motif matrix when generating multifreq matrices. This means that if the
#'    sequences used for training are completely different from the actual
#'    motif, the multifreq matrices will be as well. However this is only really
#'    a problem if you supply [add_multifreq()] with a set of sequences of the
#'    same length as the motif; in this case [add_multifreq()] is forced to
#'    create the multifreq matrices from these sequences. Otherwise
#'    [add_multifreq()] will scan the input sequences for the motif and use the
#'    best matches to construct the multifreq matrices.
#'
#'    This 'multifreq' representation is only really useful within the
#'    \pkg{universalmotif} enrivonment. Despite this, if you wish it can be
#'    preserved in text using [write_motifs()].
#'
#'    Note: the number of rows for each k-let matrix is n^k, with n being the
#'    number of letters in the alphabet being used. This means that the size
#'    of the k-let matrix can become quite large as k increases. For example,
#'    if one were to wish to represent a DNA motif of length 10 as a 10-let,
#'    this would require a matrix with 1,048,576 rows (though at this point
#'    if what you want is to search for exact sequence matches,
#'    the motif format itself is not very useful).
#'
#' @return A [universalmotif-class] object with filled `multifreq` slot.
#'
#' @examples
#' sequences <- create_sequences(seqlen = 10)
#' motif <- create_motif()
#' motif.trained <- add_multifreq(motif, sequences, add.k = 2:4)
#' ## peek at the 2-let matrix:
#' motif.trained["multifreq"]$`2`
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [scan_sequences()], [convert_motifs()], [write_motifs()] 
#' @export
add_multifreq <- function(motif, sequences, add.k = 2:3, RC = FALSE,
                          threshold = 0.01, threshold.type = "logodds",
                          motifs.perseq = 1) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(threshold.type = args$threshold.type), 1,
                                 FALSE, "character")
  num_check <- check_fun_params(list(add.k = args$add.k, threshold = args$threshold,
                                     motifs.perseq = args$motifs.perseq),
                                c(0, 1, 1), c(FALSE, FALSE, FALSE),
                                "numeric")
  logi_check <- check_fun_params(list(RC = args$RC), 1, FALSE, "logical")
  s4_check <- check_fun_params(list(sequences = args$sequences),
                               c(1, 1), c(FALSE, FALSE), "S4")
  all_checks <- c(char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motif <- convert_motifs(motif)
  motif <- convert_type(motif, "PPM")
  
  if (all(ncol(motif["motif"]) != unique(width(sequences)))) {

    seq.names <- names(sequences)
    if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
    seq.res <- scan_sequences(motif, sequences, threshold = threshold, RC = RC,
                              threshold.type = threshold.type, progress = FALSE,
                              verbose = 0)

    seqs.out <- vector("list", length(sequences))

    for (i in seq_len(length(sequences))) {
      seq.out <- seq.res[seq.res$sequence == seq.names[i], ]
      seq.out <- seq.out[order(seq.out$score, decreasing = TRUE), ]
      seq.out <- seq.out[seq_len(motifs.perseq), ]
      seqs.out[[i]] <- seq.out
    }

    seqs.out <- do.call(rbind, seqs.out)
    seqs.out <- seqs.out$match

  } else {

    seqs.out <- sequences

  }

  seqs.out <- seqs.out[!is.na(seqs.out)]

  if (length(seqs.out) == 0)
    stop("No motif matches found in sequences; consider lowering the minimum threshold")

  # multifreq <- vector("list", length(add.k))
  if (sequences@elementType == "DNAString") {
    sequences <- DNAStringSet(seqs.out)
    # for (i in seq_along(add.k)) {
      # multifreq[[i]] <- add_multi(motif["bkg"], sequences, add.k[i])
      # multifreq[[i]] <- add_multi_fix_freqs(motif["motif"], multifreq[[i]])
    # }
    multifreq <- lapply(add.k, function(x) add_multi(motif["bkg"], sequences, x))
  } else {
    alph <- rownames(motif["motif"])
    # for (i in seq_along(add.k)) {
      # multifreq[[i]] <- add_multi_ANY(sequences, add.k[i], alph)
      # multifreq[[i]] <- add_multi_fix_freqs(motif["motif"], multifreq[[i]])
    # }
    multifreq <- lapply(add.k, function(x) add_multi_ANY(sequences, x, alph))
  }

  names(multifreq) <- add.k
  prev.multifreq <- motif["multifreq"]
  if (length(prev.multifreq) > 0) {
    if (any(names(prev.multifreq) %in% names(multifreq))) {
      warning("Overwriting previous `multifreq`: ",
              paste(names(prev.multifreq)[names(prev.multifreq) %in% names(multifreq)],
                    collapse = ", "))
      prev.multifreq <- prev.multifreq[!names(prev.multifreq) %in% names(multifreq)]
      if (length(prev.multifreq) > 0) {
        multifreq <- c(prev.multifreq, multifreq)
        multifreq <- multifreq[sort(names(multifreq))]
      }
    } else {
      multifreq <- c(prev.multifreq, multifreq)
      multifreq <- multifreq[sort(names(multifreq))]
    }
  }
  motif@multifreq <- multifreq

  motif

}

# add_multi_fix_freqs <- function(mot.mat, multi.mat) {
  #
  # for (i in seq_len(ncol(multi.mat))) {
    # multi.mat[, i] <- add_multi_fix_freqs_sub(mot.mat[, i], multi.mat[, i])
  # }
#
  # multi.mat[is.na(multi.mat)] <- 0
  # multi.mat
#
# }

# add_multi_fix_freqs_sub <- function(mot.pos, multi.pos) {
#
  # n1 <- length(mot.pos)
  # n2 <- length(multi.pos)
#
  # split.i <- split(seq_len(n2), cut(seq_len(n2), n1, labels = FALSE))
  # multi.split <- vector("list", length = n1)
#
  # for (i in seq_along(multi.split)) {
    # multi.split[[i]] <- multi.pos[split.i[[i]]]
    # norm.i <- mot.pos[i] / sum(multi.split[[i]])
    # multi.split[[i]] <- multi.split[[i]] * norm.i
  # }
#
  # do.call(c, multi.split)
#
# }

add_multi <- function(bkg, sequences, k) {

  seq.width <- unique(width(sequences))
  if (seq.width < k - 1) {
    warning("motif is not long enough for k = ", k)
    return(matrix())
  }

  emissions <- matrix(nrow = 4^k, ncol = seq.width - k + 1)
  
  multi_rows <- matrix(nrow = k, ncol = 4^k)
  for (i in seq_len(k)) {
    j <- rep(DNA_BASES, each = 4^(k - i + 1) / 4)
    if (length(j) != 4^k) j <- rep(j, 4^k / length(j))
    multi_rows[i, ] <- j
  }
  multi_rows <- apply(multi_rows, 2, paste, collapse = "")

  rownames(emissions) <- multi_rows
  colnames(emissions) <- seq_len(ncol(emissions))

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])
  seqs.split <- t(seqs.split)

  for (i in seq_len(seq.width - k + 1)) {
    current.seqs <- seqs.split[, i:(i + k - 1)]
    current.seqs <- apply(current.seqs, 1, paste, collapse = "")
    current.seqs <- DNAStringSet(current.seqs)
    emissions.i <- colSums(oligonucleotideFrequency(current.seqs, k, 1))
    emissions.i <- emissions.i / sum(emissions.i)
    emissions[, i] <- emissions.i
  }

  emissions

}

add_multi_ANY <- function(sequences, k, alph) {

  seq.width <- unique(width(sequences))
  if (seq.width < k - 1) {
    warning("motif is not long enough for k = ", k)
    return(matrix())
  }

  alph.len <- length(alph)
  emissions <- matrix(rep(0, alph.len^k * (seq.width - k + 1)),
                      nrow = alph.len^k, ncol = seq.width - k + 1)

  alph.comb <- as.matrix(expand.grid(rep(list(alph), k)))
  alph.comb <- apply(alph.comb, 1, paste, collapse = "")
  alph.comb <- sort(alph.comb)

  rownames(emissions) <- alph.comb
  colnames(emissions) <- seq_len(ncol(emissions))

  seqs.split <- matrix(as.character(sequences), ncol = 1)
  seqs.split <- apply(seqs.split, 1, function(x) strsplit(x, "")[[1]])

  seq.list <- lapply(seq_len(ncol(seqs.split)),
                     function(x) single_to_k(seqs.split[, x], k))

  seq.list.i <- character(length(seq.list[[1]]))
  for (i in seq_along(seq.list[[1]])) {
    for (j in seq_along(seq.list)) {
      seq.list.i[j] <- seq.list[[j]][i]
    }
    seq.list.t <- as.matrix(table(seq.list.i))
    for (j in seq_len(nrow(seq.list.t))) {
      emissions[rownames(emissions) == rownames(seq.list.t)[j], i] <- seq.list.t[j, ]
    }
    emissions[, i] <- emissions[, i] / sum(emissions[, i])
  }

  emissions

}
