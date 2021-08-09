#' Add multi-letter information to a motif.
#'
#' If the original sequences are available for a particular motif, then they
#' can be used to generate higher-order PPM matrices. See the "Motif import,
#' export, and manipulation" vignette for more information.
#'
#' @param motif See [convert_motifs()] for acceptable formats. If the
#'    motif is not a [universalmotif-class] motif, then it will be
#'    converted.
#' @param sequences \code{\link{XStringSet}} The alphabet must match
#'    that of the motif. If
#'    these sequences are all the same length as the motif, then they are all
#'    used to generate the multi-freq matrices. Otherwise
#'    [scan_sequences()] is first run to find the best sequence stretches
#'    within these.
#' @param add.k `numeric(1)` The k-let lengths to add.
#' @param motifs.perseq `numeric(1)` If [scan_sequences()] is run,
#'    then this indicates how many hits from each sequence is to be used.
#' @param add.bkg `logical(1)` Indicate whether to add corresponding higher
#'    order background information to the motif. Can sometimes be detrimental
#'    when the input consists of few short sequences, which can increase
#'    the likelihood of adding zero or near-zero probabilities.
#'
#' @details
#'    See [scan_sequences()] for more info on scanning parameters.
#'
#'    At each position in the motif, then the probability of each k-let
#'    covering from the initial position to `ncol - 1` is calculated. Only
#'    positions within the motif are considered: this means that the
#'    final k-let probability matrix will have `ncol - 1` fewer columns.
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
#'    same length as the motif. In this case [add_multifreq()] is forced to
#'    create the multifreq matrices from these sequences. Otherwise
#'    [add_multifreq()] will scan the input sequences for the motif and use the
#'    best matches to construct the multifreq matrices.
#'
#'    This 'multifreq' representation is only really useful within the
#'    \pkg{universalmotif} environment. Despite this, if you wish it can be
#'    preserved in text using [write_motifs()].
#'
#'    ## A note on motif size
#'    The number of rows for each k-let matrix is `n^k`, with `n` being the
#'    number of letters in the alphabet being used. This means that the size
#'    of the k-let matrix can become quite large as k increases. For example,
#'    if one were to wish to represent a DNA motif of length 10 as a 10-let,
#'    this would require a matrix with 1,048,576 rows (though at this point
#'    if what you want is to search for exact sequence matches,
#'    the motif format itself is not very useful).
#'
#' @return A [universalmotif-class] object with filled `multifreq` slot. The
#'    `bkg` slot is also expanded with corresponding higher order probabilities
#'    if `add.bkg = TRUE`.
#'
#' @examples
#' sequences <- create_sequences(seqlen = 10)
#' motif <- create_motif()
#' motif.trained <- add_multifreq(motif, sequences, add.k = 2:4)
#' ## peek at the 2-let matrix:
#' motif.trained["multifreq"]$`2`
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @seealso [scan_sequences()], [convert_motifs()], [write_motifs()]
#' @inheritParams scan_sequences
#' @export
add_multifreq <- function(motif, sequences, add.k = 2:3, RC = FALSE,
                          threshold = 0.001, threshold.type = "pvalue",
                          motifs.perseq = 1, add.bkg = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!threshold.type %in% c("logodds", "pvalue", "logodds.abs")) {
    threshold.type_check <- paste0(" * Incorrect 'threshold.type': expected ",
                                   "`logodds`, `logodds.abs`, or `pvalue`; got `",
                                   threshold.type, "`")
    all_checks <- c(all_checks, threshold.type_check)
  }
  char_check <- check_fun_params(list(threshold.type = args$threshold.type), 1,
                                 FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(add.k = args$add.k, threshold = args$threshold,
                                     motifs.perseq = args$motifs.perseq),
                                c(0, 1, 1), c(FALSE, FALSE, FALSE),
                                TYPE_NUM)
  logi_check <- check_fun_params(list(RC = args$RC), 1, FALSE, TYPE_LOGI)
  s4_check <- check_fun_params(list(sequences = args$sequences),
                               1, FALSE, TYPE_S4)
  all_checks <- c(all_checks, char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motif <- convert_motifs(motif)
  motif <- convert_type_internal(motif, "PPM")

  if (all(ncol(motif@motif) != unique(width(sequences)))) {

    seq.names <- names(sequences)
    if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
    seq.res <- scan_sequences(motif, sequences, threshold = threshold, RC = RC,
                              threshold.type = threshold.type,
                              verbose = 0)

    seqs.out <- vector("list", length(sequences))

    for (i in seq_len(length(sequences))) {
      seq.out <- seq.res[seq.res$sequence == seq.names[i], ]
      seq.out <- seq.out[order(seq.out$score, decreasing = TRUE), ]
      if (nrow(seq.out) > motifs.perseq)
        seq.out <- seq.out[seq_len(motifs.perseq), ]
      seqs.out[[i]] <- seq.out
    }

    seqs.out <- do.call(rbind, seqs.out)
    seqs.out <- seqs.out$match

  } else {

    seqs.out <- as.character(sequences)

  }

  seqs.out <- seqs.out[!is.na(seqs.out)]

  seqlen <- unique(vapply(seqs.out, nchar, integer(1)))
  if (length(seqlen) > 1)
    stop(wmsg("something went wrong with extracting motif matches, ",
              "check input carefully"))

  counter  <- 1
  for (i in seq_along(add.k)) {
    if (seqlen < add.k[i]) {
      warning("motif is not long enough for use.freq>=", add.k[i],
              immediate. = TRUE)
      add.k <- add.k[seq_len(counter)]
    }
    counter <- counter + 1
  }


  if (length(seqs.out) == 0)
    stop("No motif matches found in sequences; consider lowering the minimum threshold")

  alph <- rownames(motif@motif)
  multifreq <- lapply(add.k, function(x) add_multi_cpp(seqs.out, x, alph))

  names(multifreq) <- add.k
  prev.multifreq <- motif@multifreq
  if (length(prev.multifreq) > 0) {
    if (any(names(prev.multifreq) %in% names(multifreq))) {
      warning("Overwriting previous `multifreq`: ",
              paste0(names(prev.multifreq)[names(prev.multifreq) %in% names(multifreq)],
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

  if (add.bkg) {
    new.bkg <- get_bkg(sequences, k = add.k, RC = RC)
    new.bkg <- structure(new.bkg$probability, names = new.bkg$klet)
    motif@bkg <- c(motif@bkg, new.bkg)
  }

  motif

}
