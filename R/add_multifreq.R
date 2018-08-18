#' Add multi-nucleotide information to a motif.
#'
#' If the original sequences are available for a particular motif, then they
#' can be used to generate higher-order PPM matrices.
#'
#' @param motif See \code{\link{convert_motifs}} for acceptable formats. If the
#'    motif is not a \linkS4class{universalmotif} motif, then it will be
#'    converted.
#' @param sequences \code{XStringSet} The alphabet must match that of the motif. If
#'    these sequences are all the same length as the motif, then they are all
#'    used to generate the multi-freq matrices. Otherwise
#'    \code{\link{scan_sequences}} is first run to find the right sequence.
#' @param add.k \code{numeric(1)} The k-let lengths to add.
#' @param threshold \code{numeric(1)} See \code{\link{scan_sequences}}.
#' @param threshold.type \code{character(1)} See \code{link{scan_sequences}}.
#' @param RC \code{logical(1)} See \code{link{scan_sequences}}.
#' @param motifs.perseq \code{numeric(1)} If \code{\link{scan_sequences}} is run, 
#'    then this indicates how many hits from each sequence is to be used.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @details
#'    At each position in the motif, then the probability of each k-let 
#'    covering from the initial position to k - 1 is calculated. Only
#'    positions within the motif are considered; this means that the
#'    final k-let probability matrix will have k - 1 fewer columns.
#'    Calculating k-let probabilities for the missing columns would be
#'    trivial however, as you would only need the background frequencies.
#'    Since these would not be useful for \code{\link{scan_sequences}}
#'    though, they are not calculated.
#'
#'    Note: the number of rows for each k-let matrix is n^k, with n being the
#'    number of letters in the alphabet being used. This means that the size
#'    of the k-let matrix can become quite large as k increases. For example,
#'    if one were to wish to represent a DNA motif of length 10 as a 10-let,
#'    this would require a matrix with 1,048,576 rows.
#'
#' @return A \linkS4class{universalmotif} object with filled 'multifreq' slot.
#'
#' @examples
#' sequences <- create_sequences(seqlen = 10)
#' motif <- create_motif()
#' motif.trained <- add_multifreq(motif, sequences, add.k = 2:4)
#' ## peak at the 2-let matrix:
#' motif.trained["multifreq"]$`2`
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso \code{\link{scan_sequences}}, \code{link{convert_motifs}} 
#' @export
add_multifreq <- function(motif, sequences, add.k = 2:3, RC = FALSE,
                          threshold = 0.01, threshold.type = "logodds",
                          motifs.perseq = 1, BPPARAM = SerialParam()) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(threshold.type = args$threshold.type), 1,
                                 FALSE, "character")
  num_check <- check_fun_params(list(add.k = args$add.k, threshold = args$threshold,
                                     motifs.perseq = args$motifs.perseq),
                                c(0, 1, 1), c(FALSE, FALSE, FALSE),
                                "numeric")
  logi_check <- check_fun_params(list(RC = args$RC), 1, FALSE, "logical")
  s4_check <- check_fun_params(list(sequences = args$sequences,
                                    BPPARAM = args$BPPARAM),
                               c(1, 1), c(FALSE, FALSE), "S4")
  all_checks <- c(char_check, num_check, logi_check, s4_check)
  all_checks <- paste(all_checks, collapse = "\n")
  if (length(all_checks) > 0 && all_checks[1] != "") stop(all_checks)
  #---------------------------------------------------------

  motif <- convert_motifs(motif, BPPARAM = BPPARAM)
  
  if (all(ncol(motif["motif"]) != unique(width(sequences)))) {

    seq.names <- names(sequences)
    if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
    seq.res <- scan_sequences(motif, sequences, threshold = threshold, RC = RC,
                              threshold.type = threshold.type, BPPARAM = BPPARAM,
                              verbose = FALSE)

    seqs.out <- vector("list", length(sequences))

    for (i in seq_len(length(sequences))) {
      seq.out <- seq.res[seq.res$sequence == seq.names[i], ]
      seq.out <- seq.out[order(seq.out$score, decreasing = TRUE), ]
      seq.out <- seq.out[seq_len(motifs.perseq), ]
      seqs.out[[i]] <- seq.out
    }

    seqs.out <- do.call(cbind, seqs.out)
    seqs.out <- seqs.out$match
    seqs.out <- levels(seqs.out)

  } else {

    seqs.out <- sequences

  }

  multifreq <- vector("list", length(add.k))
  if (sequences@elementType == "DNAString") {
    sequences <- DNAStringSet(seqs.out)
    for (i in seq_along(add.k)) {
      multifreq[[i]] <- add_multi(motif["bkg"], sequences, add.k[i])
    }
  } else {
    alph <- rownames(motif["motif"])
    for (i in seq_along(add.k)) {
      multifreq[[i]] <- add_multi_ANY(sequences, add.k[i], alph)
    }
  }

  names(multifreq) <- add.k
  prev.multifreq <- motif["multifreq"]
  if (length(prev.multifreq) > 0) {
    if (any(names(prev.multifreq) %in% names(multifreq))) {
      warning("Overwriting previous `multifreq`: ",
              names(prev.multifreq)[names(prev.multifreq) %in% names(multifreq)])
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
