#' Scan sequences for matches to input motifs.
#'
#' For sequences of any alphabet, scan them using the PWM matrices of
#' a set of input motifs.
#'
#' @param motifs See `convert_motifs()` for acceptable motif formats.
#' @param sequences \code{\link{XStringSet}} Sequences to scan. Alphabet
#'    should match motif.
#' @param threshold `numeric(1)` Between 0 and 1. See details.
#' @param threshold.type `character(1)` One of `c('logodds', 'pvalue')`.
#'    See details.
#' @param RC `logical(1)` If `TRUE`, check reverse complement of input
#'    sequences.
#' @param use.freq `numeric(1)` The default, 1, uses the motif matrix (from
#'    the `motif['motif']` slot) to search for sequences. If a higher
#'    number is used, then the matching k-let matrix from the
#'    `motif['multifreq']` slot is used. See [add_multifreq()].
#' @param verbose `numeric(1)` Describe progress, from none (`0`) to very
#'    verbose (`3`).
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#'    Set to `FALSE` if `verbose = 0`.
#' @param BP `logical(1)` Allows for the use of \pkg{BiocParallel} within
#'    [scan_sequences()]. See [BiocParallel::register()] to change the
#'    default backend. Setting `BP = TRUE` is only recommended for
#'    exceptionally large jobs. Keep in mind however that this function
#'    will not attempt to limit its memory usage. Furthermore, the
#'    behaviour of `porgress = TRUE` is changed if `BP = TRUE`; the
#'    default \pkg{BiocParallel} progress bar will be shown (which
#'    unfortunately is much less informative).
#'
#' @return `data.frame` with each row representing one hit; if the input
#'    sequences are \code{\link{DNAStringSet}} or
#'    \code{\link{RNAStringSet}}, then an
#'    additional column with the strand is included.
#'
#' @details
#'    Similar to [Biostrings::matchPWM()], the scanning method uses
#'    logodds scoring. (To see the scoring matrix for any motif, simply
#'    run `convert_type(motif, "PWM")`; for a `multifreq` scoring
#'    matrix: `apply(motif["multifreq"]$`2`, 2, ppm_to_pwm)`). In order
#'    to score a sequence, at each position within a sequence of length equal
#'    to the length of the motif, the scores for each base are summed. If the
#'    score sum is above the desired threshold, it is kept.
#'
#'    If `threshold.type = 'logodds'`, then to calculate the minimum
#'    allowed score the total possible score for a motif is multiplied
#'    by the value set by `threshold`. To determine the maximum and minimum
#'    possible scores a motif (of type PWM), run
#'    `sum(apply(motif['motif'], 2, max))` and 
#'    `sum(apply(motif['motif'], 2, min))`. If \code{threshold.type = 'pvalue'},
#'    then threshold logodds scores are generated using [motif_pvalue()].
#'
#'    Non-standard letters (such as "N", "+", "-", ".", etc in `DNAString`
#'    objects) will be safely ignored, resulting only in a warning and a very
#'    minor performance cost. This can used to scan
#'    masked sequences. See \code{\link[Biostrings:maskMotif]{Biostrings::mask()}}
#'    for masking sequences
#'    (generating `MaskedXString` objects), and [Biostrings::injectHardMask()]
#'    to recover masked `XStringSet` objects for use with [scan_sequences()].
#'
#' @examples
#' ## any alphabet can be used
#' \dontrun{
#' set.seed(1)
#' alphabet <- paste(c(letters), collapse = "")
#' motif <- create_motif("hello", alphabet = alphabet)
#' sequences <- create_sequences(alphabet, seqnum = 1000, seqlen = 100000)
#' scan_sequences(motif, sequences)
#' }
#'
#' ## Sequence masking:
#' library(Biostrings)
#' data(ArabidopsisMotif)
#' data(ArabidopsisPromoters)
#' seq <- ArabidopsisPromoters[[1]]  # Only works for XString, not XStringSet
#' seq <- mask(seq, pattern = "AAAA")  # MaskedDNAString class
#' seq <- injectHardMask(seq, letter = "+")  # Recover XString
#' seq <- DNAStringSet(seq)  # scan_sequences() needs XStringSet
#' scan_sequences(ArabidopsisMotif, seq, verbose = 0, progress = FALSE)
#' # A warning regarding the presence of non-standard letters will be given,
#' # but can be safely ignored in this case.
#'
#' @references
#'    \insertRef{biostrings}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [add_multifreq()], [Biostrings::matchPWM()],
#'    [enrich_motifs()], [motif_pvalue()]
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.001,
                           threshold.type = "pvalue", RC = FALSE,
                           use.freq = 1, verbose = 1,
                           progress = TRUE, BP = FALSE) {

  # TODO: Work with Masked*String objects. Masked letters show up as "#" after
  #       as.character() calls, which should just cause scan_sequences() to
  #       ignore these and work as intended. For now, just using "-", "." or
  #       "+" within DNA/RNA/AAStringSet objects will work the same.

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!threshold.type %in% c("logodds", "pvalue")) {
    threshold.type_check <- paste0(" * Incorrect 'threshold.type': expected",
                                   "`logodds` or `pvalue`; got `",
                                   threshold.type, "`")
    all_checks <- c(all_checks, threshold.type_check)
  }
  char_check <- check_fun_params(list(threshold.type = args$threshold.type),
                                 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(threshold = args$threshold,
                                     use.freq = args$use.freq,
                                     verbose = args$verbose),
                                c(0, 1, 1), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(RC = args$RC),
                                 numeric(), logical(), TYPE_LOGI)
  s4_check <- check_fun_params(list(sequences = args$sequences), numeric(),
                               logical(), TYPE_S4)
  all_checks <- c(all_checks, char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (verbose <= 0) progress <- FALSE

  if (verbose > 2) {
    cat(" * Input parameters\n")
    cat("   * motifs:              ", deparse(substitute(motifs)), "\n")
    cat("   * sequences:           ", deparse(substitute(sequences)), "\n")
    cat("   * threshold:           ", ifelse(length(threshold) > 1, "...",
                                             threshold), "\n")
    cat("   * threshold.type:      ", threshold.type, "\n")
    cat("   * RC:                  ", RC, "\n")
    cat("   * use.freq:            ", use.freq, "\n")
    cat("   * verbose:             ", verbose, "\n")
  }

  if (missing(motifs) || missing(sequences)) {
    stop("need both motifs and sequences")
  }

  if (verbose > 0) cat(" * Processing motifs\n")

  if (verbose > 1) cat("   * Scanning", length(motifs),
                       ifelse(length(motifs) > 1, "motifs\n", "motif\n"))

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  motifs <- lapply(motifs, normalize)
  motifs <- convert_type_internal(motifs, "PWM")

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  mot.pwms <- lapply(motifs, function(x) x@motif)
  mot.lens <- vapply(mot.pwms, ncol, numeric(1))
  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) != 1) stop("can only scan using one alphabet")
  mot.alphs <- unique(mot.alphs)
  if (verbose > 1) cat("   * Motif alphabet:", mot.alphs, "\n")
  alph <- mot.alphs
  mot.bkgs <- lapply(motifs, function(x) x@bkg[rownames(x@motif)])
  seq.lens <- width(sequences)

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))

  seq.alph <- seqtype(sequences)
  if (seq.alph != "B" && seq.alph != mot.alphs)
    stop("Motif and Sequence alphabets do not match")
  else if (seq.alph == "B")
    seq.alph <- mot.alphs

  if (use.freq > 1) {
    if (any(vapply(motifs, function(x) length(x@multifreq) == 0, logical(1))))
      stop("missing multifreq slots")
    check_multi <- vapply(motifs,
                          function(x) any(names(x@multifreq) %in%
                                          as.character(use.freq)),
                          logical(1))
    if (!any(check_multi)) stop("not all motifs have correct multifreqs")
  }

  for (i in seq_along(motifs)) {
    if (motifs[[i]]@pseudocount == 0) {
      if (verbose > 1)
        cat("   * Setting 'pseudocount' to 1 for motif:", mot.names[i], "\n")
      motifs[[i]]["pseudocount"] <- 1
    }
    if (length(motifs[[i]]@nsites) == 0) {
      if (verbose > 1)
        cat("   * Setting 'nsites' to 100 for motif:", mot.names[i], "\n")
      motifs[[i]]["nsites"] <- 100
    }
  }

  if (use.freq == 1) {
    score.mats <- mot.pwms
  } else {
    score.mats <- lapply(motifs,
                         function(x) x@multifreq[[as.character(use.freq)]])
    for (i in seq_along(score.mats)) {
      score.mats[[i]] <- apply(score.mats[[i]], 2, ppm_to_pwmC,
                               nsites = motifs[[i]]@nsites,
                               pseudocount = motifs[[i]]@pseudocount)
    }
  }

  max.scores <- vapply(score.mats, function(x) sum(apply(x, 2, max)), numeric(1))
  min.scores <- vapply(score.mats, function(x) sum(apply(x, 2, min)), numeric(1))

#-------------------------------------------------------------------------------
# PVALUES --> THRESHOLDS

  switch(threshold.type,

    "logodds" = {

      thresholds <- ((abs(max.scores) + abs(min.scores)) * threshold) - abs(min.scores)

    },

    "pvalue" = {

      if (progress && !BP && verbose > 0)
        cat(" * Converting P-values to logodds thresholds ...")
      else if ((progress && BP && verbose > 0) || verbose > 0)
        cat(" * Converting P-values to logodds thresholds\n")
      thresholds <- vector("numeric", length(motifs))
      thresholds <- motif_pvalue(motifs, pvalue = threshold, use.freq = use.freq,
                                 k = 6, progress = progress, BP = BP)
      for (i in seq_along(thresholds)) {
        if (thresholds[i] > max.scores[i]) thresholds[i] <- max.scores[i]
      }
      if (verbose > 3) {
        for (i in seq_along(thresholds)) {
          cat("   * Motif ", mot.names[i], ": max.score = ", max.scores[i],
              ", threshold = ", thresholds[i], "\n", sep = "")
        }
      }
      thresholds <- unlist(thresholds)

    },

    stop("unknown 'threshold.type'")

  )

#-------------------------------------------------------------------------------
# SEQS --> SPLIT SEQS

  if (verbose > 0) cat(" * Processing sequences\n")

  if (verbose > 1) {
    cat("   * Number of sequences:", length(sequences), "\n")
    cat("   * Mean sequence width:", mean(width(sequences)), "\n")
  }

  seqs.aschar <- as.character(sequences)

  if (progress && !BP && verbose > 0)
    cat("   * Splitting up sequences ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Splitting up sequences\n")

  seqs.aschar <- lapply_(seqs.aschar, safeExplode, BP = BP, PB = progress)

#-------------------------------------------------------------------------------
# SEQS --> SEQ MATRICES

  seq.lens <- width(sequences)
  seq.matrices <- lapply(seq.lens, function(x) matrix(ncol = x - use.freq + 1,
                                                      nrow = use.freq))

  if (progress && !BP && verbose > 0)
    cat("   * Creating sequence matrices ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Creating sequence matrices\n")

  seq.matrices <- mapply_(scan_process_seqs, seq.matrices, seqs.aschar, BP = BP,
                           MoreArgs = list(k = use.freq), SIMPLIFY = FALSE,
                           PB = progress)

  switch(mot.alphs,
    "DNA" = {
      mot.alphs <- DNA_BASES
    },
    "RNA" = {
      mot.alphs <- RNA_BASES
    },
    "AA" = {
      mot.alphs <- AA_STANDARD
    },
    "custom" = {
      if (RC) stop("RC search is only available for DNA/RNA")
      mot.alphs <- lapply(seqs.aschar, unique)
      mot.alphs <- unique(unlist(mot.alphs))
    },
    {
      if (RC) stop("RC search is only available for DNA/RNA")
      mot.alphs <- safeExplode(mot.alphs)
    }
  )

  alph.int <- as.integer(seq_along(mot.alphs))

#-------------------------------------------------------------------------------
# STRINGS --> FACTORS

  if (progress && !BP && verbose > 0)
    cat("   * Converting sequences to factors ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Converting sequences to factors\n")

  seq.matrices <- lapply_(seq.matrices,
                          function(x) string_to_factor(x, mot.alphs),
                          BP = BP, PB = progress)

#-------------------------------------------------------------------------------
# CHECK FOR NAs

  ## BUG FIX: can't deal with NAs generated from non-standard DNA letters
  if (progress && !BP && verbose > 0)
    cat("   * Checking for non-standard letters ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Checking for non-standard letters\n")

  na.check <- lapply_(seq.matrices,
                      function(x) any(is.na(x)),
                      BP = BP, PB = progress)
  na.check <- any(as.logical(na.check))
  if (na.check) warning("Non-standard letters found, these will be ignored",
                        immediate. = TRUE, call. = FALSE)

#-------------------------------------------------------------------------------
# FACTORS --> INTEGERS

  if (progress && !BP && verbose > 0)
    cat("   * Converting sequences to integers ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Converting sequences to integers\n")

  seq.ints <- lapply_(seq.matrices,
                      function(x) LETTER_to_int(as.integer(x) - 1,
                                                use.freq, alph.int),
                      BP = BP, PB = progress)

#-------------------------------------------------------------------------------
# SCANNING

  if (progress && !RC && !BP && verbose > 0)
    cat(" * Scanning sequences for motifs ...")
  else if ((progress && RC && verbose > 0) || verbose > 0)
    cat(" * Scanning sequences for motifs\n")

  score.mats <- lapply(score.mats, numeric_to_integer_matrix)
  thresholds.int <- as.integer(thresholds * 1000)

  # Making the progress be per motif makes it kinda lame when scanning with a
  # single motif. Inversely though, scanning a single sequence with multiple
  # motifs would be awkard. Perhaps decide which to do based on whether
  # there are more motifs or sequences?

  if (RC && progress && !BP && verbose > 0)
    cat("   * Forward strand ...")
  else if (RC && verbose > 0)
    cat("   * Forward strand\n")

  if (!na.check) {
    to.keep <- lapply_(seq_along(score.mats),
                       function(x) scan_score_motif(seq.ints, score.mats[[x]],
                                                    thresholds.int[x]),
                       BP = BP, PB = progress)
  } else {
    to.keep <- lapply_(seq_along(score.mats),
                       function(x) scan_score_motif2(seq.ints, score.mats[[x]],
                                                     thresholds.int[x]),
                       BP = BP, PB = progress)
  }

  if (verbose > 2) {
    num.matches <- sum(sapply(to.keep, function(x) vapply(x, sum, integer(1))))
    cat("     * Found", num.matches,
        ifelse(num.matches == 1, "match\n", "matches\n"))
  }

  # TODO: Add a step to prune matches if they exceed a certain number. This will
  #       prevent keeping too many scores in memory, bogging down the
  #       get_res_cpp() step. 
  #       *Would this really work? Either way a big to.keep object is made..
  #
  #   matches <- lapply(to.keep, function(x) lapply(x, res_to_index))
  #   total.matches <- sum(sapply(matches, vapply(x, length, integer(1))))
  #
  #   if (total.matches > max.matches) {
  #
  #     scores <- vector("list", length(to.keep))
  #     for (i in seq_along(to.keep)) {
  #       scores[[i]] <- mapply(
  #          function(x, y) create_col_score(x, vapply(matches[[i]], length, integer(1)),
  #                                          length(to.keep[[i]]),
  #                                          sum(vapply(matches[[i]], length, integer(1)),
  #                                          seq.ints, score.mats[[i]], ncol(score.mats[[i]]),
  #                                          use.freq),
  #          to.keep[[i]], score.mats[[i]], SIMPLIFY = FALSE
  #       )
  #     }
  #
  #   }
  # 
  # Then sort and keep top X scores. If RC=TRUE, then after the RC scanning
  # combine the forward and reverse score then repeat keeping top X scores.
  # The calculated scores could be used in get_res_cpp() so that they don't have
  # to be re-calculated.

  if (RC) {

    score.mats.rc <- lapply(score.mats,
                            function(x) matrix(rev(as.numeric(x)),
                                               ncol = ncol(x)))

    if (progress && !BP && verbose > 0)
      cat("   * Reverse strand ...")
    else if ((progress && BP && verbose > 0) || verbose > 0)
      cat("   * Reverse strand\n")

    if (!na.check) {
      to.keep.rc <- lapply_(seq_along(score.mats.rc),
                            function(x) scan_score_motif(seq.ints, score.mats.rc[[x]],
                                                     thresholds.int[x]),
                            BP = BP, PB = progress)
    } else {
      to.keep.rc <- lapply_(seq_along(score.mats.rc),
                            function(x) scan_score_motif2(seq.ints, score.mats.rc[[x]],
                                                      thresholds.int[x]),
                            BP = BP, PB = progress)
    }

    if (verbose > 2) {
      num.matches.rc <- sum(sapply(to.keep.rc, function(x) vapply(x, sum, integer(1))))
      cat("     * Found", num.matches.rc,
          ifelse(num.matches.rc == 1, "match\n", "matches\n"))
    }

  }

#-------------------------------------------------------------------------------
# MAKE RESULTS

  if (progress && !RC && !BP && verbose > 0)
    cat(" * Processing results ...")
  else if ((progress && RC && verbose > 0) || verbose > 0)
    cat(" * Processing results\n")

  if (RC && progress && !BP && verbose > 0)
    cat("   * Forward strand ...")
  else if (RC && verbose > 0)
    cat("   * Forward strand\n")

  if (progress) print_pb(0)
  res <- vector("list", length(to.keep))
  res.len <- length(res)
  for (i in seq_along(res)) {
    to.keep[[i]] <- lapply(to.keep[[i]], res_to_index)
    res[[i]] <- get_res_cpp(to.keep[[i]], seqs.aschar, seq.ints, mot.lens[i],
                            min.scores[i], max.scores[i], mot.names[i],
                            seq.names, score.mats[[i]], "+",
                            seq.lens, use.freq)
    if (progress) update_pb(i, res.len)
  }

  if (RC) {

    if (progress && !BP && verbose > 0)
      cat("   * Reverse strand ...")
    else if (verbose > 0 || (progress && verbose > 0))
      cat("   * Reverse strand\n")

    if (progress) print_pb(0)
    res.rc <- vector("list", length(to.keep.rc))
    res.rc.len <- length(res.rc)
    for (i in seq_along(res.rc)) {
      to.keep.rc[[i]] <- lapply(to.keep.rc[[i]], res_to_index)
      res.rc[[i]] <- get_res_cpp(to.keep.rc[[i]], seqs.aschar, seq.ints,
                                 mot.lens[i], min.scores[i], max.scores[i],
                                 mot.names[i], seq.names, score.mats.rc[[i]],
                                 "-", seq.lens, use.freq)
      if (progress) update_pb(i, res.rc.len)
    }

    res <- c(res, res.rc)

  }

#-------------------------------------------------------------------------------
# OUTPUT RESULTS

  if (verbose > 0)
    cat(" * Generating output\n")

  res <- res_list_to_df_cpp(res)

  if (nrow(res) == 0) {
    message(" ! No matches to motifs found")
    return(invisible(NULL))
  }

  if (!alph %in% c("DNA", "RNA")) res <- res[, colnames(res) != "strand"]

  rownames(res) <- NULL

  res

}

#-------------------------------------------------------------------------------

scan_score_motif <- function(seqs, score.1, thresh) {
  lapply(seqs, function(x) scan_seq_internal(x, score.1, thresh))
}

scan_score_motif2 <- function(seqs, score.1, thresh) {
  lapply(seqs, function(x) scan_seq_internal2(x, score.1, thresh))
}

# perhaps implement this in C++?
scan_process_seqs <- function(seq.matrix, seq.aschar, k) {

  for (i in seq_len(k)) {
    to.remove <- k - i
    if (to.remove == 0) {
      seq.matrix[1, ] <- seq.aschar[seq_len(ncol(seq.matrix))]
      next
    }
    seq.i <- seq.aschar[-seq_len(to.remove)]
    seq.matrix[k - i + 1, ] <- seq.i[seq_len(ncol(seq.matrix))]
  }

  seq.matrix

}

#===============================================================================

# For use with enrich_motifs(..., return.scan.results = FALSE).
# Frankly, this isn't really worth it; the get_res_cpp() step really isn't that
# expensive. The scan_seq_internal() step is by far the slowest step in either
# version of scan_sequences.

# INPUT
#   score.mats: a list of PWM matrices
#   sequences: XStringSet object
#   thresholds.int: integer threshold (X1000)
#   k: use.freq
#   alph: character vector of alphabet letters
#   verbose: 0-3
#   mot.names: character vector of motif names
#   RC, PB, BP: TRUE/FALSE
# OUTPUT
#   List of data.frames. One list entry per motif. The data.frames have
#   columns sequence and start, as well as strand if RC=TRUE.
scan_sequences_slim <- function(score.mats, sequences, thresholds.int, k,
                                alph, verbose, mot.names, RC, PB, BP) {

  if (verbose <= 0) PB <- FALSE

  progress <- PB

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))

#-------------------------------------------------------------------------------
# SEQS --> SPLIT SEQS

  if (verbose > 0) cat(" * Processing sequences\n")

  if (verbose > 1) {
    cat("   * Number of sequences:", length(sequences), "\n")
    cat("   * Mean sequence width:", mean(width(sequences)), "\n")
  }

  seqs.aschar <- as.character(sequences)

  if (progress && !BP && verbose > 0)
    cat("   * Splitting up sequences ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Splitting up sequences\n")

  seqs.aschar <- lapply_(seqs.aschar, safeExplode, BP = BP, PB = PB)

#-------------------------------------------------------------------------------
# SEQS --> SEQ MATRICES

  seq.lens <- width(sequences)

  seq.matrices <- lapply(seq.lens, function(x) matrix(ncol = x - k + 1,
                                                      nrow = k))

  if (progress && !BP && verbose > 0)
    cat("   * Creating sequence matrices ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Creating sequence matrices\n")

  seq.matrices <- mapply_(scan_process_seqs, seq.matrices, seqs.aschar,
                          MoreArgs = list(k = k), SIMPLIFY = FALSE,
                          BP = BP, PB = PB)

  alph.int <- as.integer(seq_along(alph))

#-------------------------------------------------------------------------------
# STRINGS --> FACTORS

  if (progress && !BP && verbose > 0)
    cat("   * Converting sequences to factors ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Converting sequences to factors\n")

  seq.matrices <- lapply_(seq.matrices, function(x) string_to_factor(x, alph),
                          BP = BP, PB = PB)

#-------------------------------------------------------------------------------
# CHECK FOR NAs

  ## BUG FIX: can't deal with NAs generated from non-standard DNA letters
  if (progress && !BP && verbose > 0)
    cat("   * Checking for non-standard letters ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Checking for non-standard letters\n")

  na.check <- lapply_(seq.matrices, function(x) any(is.na(x)), BP = BP, PB = PB)
  na.check <- any(as.logical(na.check))

#-------------------------------------------------------------------------------
# FACTORS --> INTEGERS

  if (progress && !BP && verbose > 0)
    cat("   * Converting sequences to integers ...")
  else if ((progress && BP && verbose > 0) || verbose > 1)
    cat("   * Converting sequences to integers\n")

  seq.ints <- lapply_(seq.matrices, function(x) LETTER_to_int(as.integer(x) - 1,
                                                              k, alph.int),
                      BP = BP, PB = PB)

#-------------------------------------------------------------------------------
# SCANNING

  if (progress && !RC && !BP && verbose > 0)
    cat(" * Scanning sequences for motifs ...")
  else if ((progress && RC && verbose > 0) || verbose > 0)
    cat(" * Scanning sequences for motifs\n")

  if (RC && progress && !BP && verbose > 0)
    cat("   * Forward strand ...")
  else if (RC && verbose > 0)
    cat("   * Forward strand\n")

  if (!na.check) {
    to.keep <- lapply_(seq_along(score.mats),
                       function(x) scan_score_motif(seq.ints, score.mats[[x]],
                                                    thresholds.int[x]),
                       BP = BP, PB = PB)
  } else {
    to.keep <- lapply_(seq_along(score.mats),
                       function(x) scan_score_motif2(seq.ints, score.mats[[x]],
                                                     thresholds.int[x]),
                       BP = BP, PB = PB)
  }

  if (verbose > 2) {
    num.matches <- sum(sapply(to.keep, function(x) vapply(x, sum, integer(1))))
    cat("     * Found", num.matches,
        ifelse(num.matches == 1, "match\n", "matches\n"))
  }

  if (RC) {
    score.mats.rc <- lapply(score.mats, function(x) matrix(rev(as.numeric(x)),
                                                           ncol = ncol(x)))

    if (progress && !BP && verbose > 0)
      cat("   * Reverse strand ...")
    else if ((progress && BP && verbose > 0) || verbose > 0)
      cat("   * Reverse strand\n")

    if (!na.check) {
      to.keep.rc <- lapply_(seq_along(score.mats),
                            function(x) scan_score_motif(seq.ints, score.mats.rc[[x]],
                                                         thresholds.int[x]),
                            BP = BP, PB = PB)
    } else {
      to.keep.rc <- lapply_(seq_along(score.mats),
                            function(x) scan_score_motif2(seq.ints, score.mats.rc[[x]],
                                                          thresholds.int[x]),
                            BP = BP, PB = PB)
    }

    if (verbose > 2) {
      num.matches.rc <- sum(sapply(to.keep.rc, function(x) vapply(x, sum, integer(1))))
      cat("     * Found", num.matches.rc,
          ifelse(num.matches.rc == 1, "match\n", "matches\n"))
    }

  }

#-------------------------------------------------------------------------------
# MAKE RESULTS

  if (progress && !RC && !BP && verbose > 0)
    cat(" * Processing results ...")
  else if ((progress && RC && verbose > 0) || verbose > 0)
    cat(" * Processing results\n")

  if (RC && progress && !BP && verbose > 0)
    cat("   * Forward strand ...")
  else if (RC && verbose > 0)
    cat("   * Forward strand\n")

  for (i in seq_along(to.keep)) {
    to.keep[[i]] <- lapply(to.keep[[i]], res_to_index)
  }
  to.keep.lens <- lapply(to.keep, function(x) vapply(x, length, integer(1)))
  to.keep <- mapply_(function(x, y) index_list_to_df_cpp(x, seq.names, y),
                     to.keep, to.keep.lens, SIMPLIFY = FALSE, BP = BP, PB = PB)
  names(to.keep) <- mot.names

  if (RC) {

    if (progress && !BP && verbose > 0)
      cat("   * Reverse strand ...")
    else if (verbose > 0 || (progress && verbose > 0))
      cat("   * Reverse strand\n")

    for (i in seq_along(to.keep.rc)) {
      to.keep.rc[[i]] <- lapply(to.keep.rc[[i]], res_to_index)
    }
    to.keep.rc.lens <- lapply(to.keep.rc, function(x) vapply(x, length, integer(1)))
    to.keep.rc <- mapply_(function(x, y) index_list_to_df_cpp(x, seq.names, y),
                          to.keep.rc, to.keep.rc.lens, SIMPLIFY = FALSE, BP = BP, PB = PB)
    names(to.keep.rc) <- paste0(mot.names)

    # COMBINE FORWARD + REVERSE
    for (i in seq_along(to.keep)) {
      to.keep.rc[[i]]$start <- to.keep.rc[[i]]$start + ncol(score.mats[[i]]) - 1
      # to.keep[[i]]$strand <- rep("+", nrow(to.keep[[i]]))
      # to.keep.rc[[i]]$strand <- rep("-", nrow(to.keep.rc[[i]]))
      to.keep[[i]] <- rbind(to.keep[[i]], to.keep.rc[[i]])  # potential bottleneck?
    }

  }

  to.keep

}
