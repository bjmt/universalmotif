#' Enrich for input motifs in a set of sequences.
#'
#' Given a set of target and background sequences, test if the input motifs
#' are significantly enriched in the targets sequences relative to the
#' background sequences.
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param sequences \code{\link{XStringSet}} Alphabet should match motif.
#' @param bkg.sequences \code{\link{XStringSet}} Optional; if missing,
#'    [shuffle_sequences()] is used to create background sequences from
#'    the input sequences.
#' @param search.mode `character(1)` One of `c('hits', 'positional', 'both')`.
#'    See details.
#' @param max.p `numeric(1)` P-value threshold.
#' @param max.q `numeric(1)` Adjusted P-value threshold. This is only useful
#'    if multiple motifs are being enriched for.
#' @param max.e `numeric(1)`. The E-value is calculated by multiplying the adjusted
#'    P-value with the number of input motifs times two
#'    \insertCite{meme2}{universalmotif}.
#' @param qval.method `character(1)` See [stats::p.adjust()].
#' @param positional.test `character(1)` One of `c('t.test', 'wilcox.test',
#'    'chisq.test', 'shapiro.test')`. If using the Shapiro test for
#'    normality, then only the input sequences are tested for positionality;
#'    the background sequences are ignored. See [stats::t.test()],
#'    [stats::wilcox.test()], [stats::chisq.test()], [stats::shapiro.test()].
#' @param threshold `numeric(1)` Between 1 and 0. See [scan_sequences()].
#' @param threshold.type `character(1)` One of `c('logodds', 'pvalue')`. See
#'    [scan_sequences()].
#' @param verbose `numeric(1)` 0 for no output, 4 for max verbosity.
#' @param RC `logical(1)` Whether to consider the reverse complement of the
#'    sequences. Only available for \code{\link{DNAStringSet}},
#'    \code{\link{RNAStringSet}} sequences.
#' @param use.freq `numeric(1)` If the `multifreq` slot of the motifs are filled,
#'    then they can be used to scan the sequences. See [scan_sequences()].
#' @param shuffle.k `numeric(1)` The k-let size to use when shuffling input
#'    sequences. Only used if no background sequences are input. See
#'    [shuffle_sequences()].
#' @param shuffle.method `character(1)` One of `c('markov', 'linear', 'random')`.
#'    See [shuffle_sequences()]. The `'random'` method is deprecated and 
#'    will be removed in the next minor version.
#' @param shuffle.leftovers `character(1)` One of `c('asis', 'first', 'split', 'discard')`.
#'    Only used if `shuffle.method = 'random'`. See [shuffle_sequences()].
#' @param return.scan.results `logical(1)` Return output from
#'    [scan_sequences()]. For large jobs, leaving this as
#'    `FALSE` can save a small amount time by preventing construction of the complete 
#'    results `data.frame` from [scan_sequences()].
#' @param progress `logical(1)` Show progress. Note recommended if `BP = TRUE`.
#'    Set to `FALSE` if `verbose = 0`
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [enrich_motifs()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for exceptionally large jobs
#'    (be wary of memory usage however, as [enrich_motifs()] does not try and
#'    limit itself in this regard).
#'    Furthermore, the behaviour of `progress = TRUE` is
#'    changed if `BP = TRUE`; the default \pkg{BiocParallel} progress bar will
#'    be shown (which unfortunately is much less informative).
#'
#' @return `data.frame` Motif enrichment results. The resulting 
#'    `data.frame` contains the following columns:
#'     * `motif` Motif name.
#'     * `total.seq.hits` Total number of matches across all target
#'       sequences.
#'     * `num.seqs.hits` Number of target sequences which contain matches.
#'     * `num.seqs.total` Number of target sequences.
#'     * `total.bkg.hits` Total number of matches across all background
#'       sequences.
#'     * `num.bkg.hits` Number of background sequences which contain
#'       matches.
#'     * `num.bkg.total` Number of background sequences.
#'     * `Pval.hits` P-value of enrichment. Only shown if
#'       `search.mode = c('hits', 'both')`.
#'     * `Qval.hits` Q-val of enrichment. Only shown if
#'       `search.mode = c('hits', 'both')`.
#'     * `Eval.hits` E-val of enrichment. Only shown if
#'       `search.mode = c('hits', 'both')`.
#'     * `Pval.pos` P-value of positional comparison. Only
#'       shown if `search.mode = c('positional', 'both')`.
#'     * `Qval.pos` Q-value of positional comparison. Only
#'       shown if `search.mode = c('positional', 'both')`.
#'     * `Eval.pos` E-value of positional comparison. Only
#'       shown if `search.mode = c('positional', 'both')`.
#'
#' @details
#' To find enriched motifs, [scan_sequences()] is run on both 
#' target and background sequences. If `search.mode = 'hits'`,
#' [stats::fisher.test()] is run to test for enrichment. If
#' `search.mode = 'positional'`, then the test as set by
#' `positional.test` is run to check for positional differences
#' between target and background sequences. However if
#' `positional.test = 'shapiro.test'`, then only target sequence
#' hits are considered.
#'
#' @examples
#' data(ArabidopsisPromoters)
#' data(ArabidopsisMotif)
#' enrich_motifs(ArabidopsisMotif, ArabidopsisPromoters, threshold = 0.01)
#'
#' @references
#'    \insertRef{meme2}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay \email{b2tremblay@@uwaterloo.ca}
#' @seealso [scan_sequences()], [shuffle_sequences()],
#'    [add_multifreq()], [motif_pvalue()]
#' @export
enrich_motifs <- function(motifs, sequences, bkg.sequences, search.mode = "hits",
                          max.p = 10e-6, max.q = 10e-6, max.e = 10e-4,
                          qval.method = "fdr", positional.test = "t.test",
                          threshold = 0.001, threshold.type = "pvalue",
                          verbose = 1, RC = FALSE, use.freq = 1,
                          shuffle.k = 2, shuffle.method = "linear",
                          shuffle.leftovers = "asis",
                          return.scan.results = FALSE, progress = TRUE,
                          BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!search.mode %in% c("hits", "positional", "both")) {
    search.mode_check <- paste0(" * Incorrect 'search.mode': expected `hits`,",
                                " `positional` or `both`; got `",
                                search.mode, "`")
    search.mode_check <- wmsg2(search.mode_check, 4, 2)
    all_checks <- c(all_checks, search.mode_check)
  }
  if (!qval.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                          "fdr", "none")) {
    qval.method_check <- paste0(" * Incorrect 'qval.method': expected `holm`, ",
                                "`hochberg`, `hommel`, `bonferroni`, `BH`, `BY`, ",
                                "`fdr` or `none`; got `",
                                qval.method, "`")
    qval.method_check <- wmsg2(qval.method_check, 4, 2)
    all_checks <- c(all_checks, qval.method_check)
  }
  if (!positional.test %in% c("t.test", "wilcox.test", "chisq.test",
                              "shapiro.test")) {
    positional.test_check <- paste0(" * Incorrect 'positional.test': expected ",
                                    "`t.test`, `wilcox.test`, `chisq.test` or ",
                                    "`shapiro.test`; got `",
                                    positional.test, "`")
    positional.test_check <- wmsg2(positional.test_check, 4, 2)
    all_checks <- c(all_checks, positional.test_check)
  }
  if (!threshold.type %in% c("logodds", "pvalue")) {
    threshold.type_check <- paste0(" * Incorrect 'threshold.type': expected ",
                                   "`logodds` or `pvalue`; got `",
                                   threshold.type, "`")
    threshold.type_check <- wmsg2(threshold.type_check, 4, 2)
    all_checks <- c(all_checks, threshold.type_check)
  }
  if (!shuffle.method %in% c("markov", "linear", "random")) {
    shuffle.method_check <- paste0(" * Incorrect 'shuffle.method': expected ",
                                   "`markov`, `linear` or `random`; got `",
                                   shuffle.method, "`")
    shuffle.method_check <- wmsg2(shuffle.method_check, 4, 2)
    all_checks <- c(all_checks, shuffle.method_check)
  }
  if (!shuffle.leftovers %in% c("asis", "first", "split", "discard")) {
    shuffle.leftovers_check <- paste0(" * Incorrect 'shuffle.leftovers': expected ",
                                      "`asis`, `first`, `split` or `discard`; got `",
                                      shuffle.leftovers, "`")
    shuffle.leftovers_check <- wmsg2(shuffle.leftovers_check)
    all_checks <- c(all_checks, shuffle.leftovers_check)
  }
  char_check <- check_fun_params(list(search.mode = args$search.mode,
                                      qval.method = args$qval.method,
                                      positional.test = args$positional.test,
                                      threshold.type = args$threshold.type,
                                      shuffle.method = args$shuffle.method,
                                      shuffle.leftovers = args$shuffle.leftovers),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(max.p = args$max.p, max.q = args$max.q,
                                     max.e = args$max.e,
                                     threshold = args$threshold,
                                     verbose = args$verbose, use.freq = args$use.freq,
                                     shuffle.k = args$shuffle.k),
                                c(1, 1, 1, 0, 1, 1, 1), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(RC = args$RC, progress = args$progress,
                                      return.scan.results = args$return.scan.results,
                                      BP = args$BP),
                                 numeric(), logical(), TYPE_LOGI)
  s4_check <- check_fun_params(list(sequences = args$sequences,
                                    bkg.sequences = args$bkg.sequences),
                               numeric(), c(FALSE, TRUE, FALSE), TYPE_S4)
  all_checks <- c(all_checks, char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (verbose <= 0) progress <- FALSE

  if (verbose > 2) {
    cat(" > Input parameters\n")
    cat("   > motifs:              ", deparse(substitute(motifs)), "\n")
    cat("   > sequences:           ", deparse(substitute(sequences)), "\n")
    cat("   > bkg.sequences:       ", ifelse(!missing(bkg.sequences),
                                        deparse(substitute(bkg.sequences)),
                                        "none"), "\n")
    cat("   > search.mode:         ", search.mode, "\n")
    cat("   > max.p:               ", max.p, "\n")
    cat("   > max.q:               ", max.q, "\n")
    cat("   > max.e:               ", max.e, "\n")
    cat("   > qval.method:         ", qval.method, "\n")
    cat("   > positional.test:     ", positional.test, "\n")
    cat("   > threshold:           ", threshold, "\n")
    cat("   > threshold.type:      ", threshold.type, "\n")
    cat("   > verbose:             ", verbose, "\n")
    cat("   > RC:                  ", RC, "\n")
    cat("   > use.freq:            ", use.freq, "\n")
    cat("   > shuffle.k:           ", shuffle.k, "\n")
    cat("   > shuffle.method:      ", shuffle.method, "\n")
    cat("   > return.scan.results: ", return.scan.results, "\n")
  }

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PWM")

  if (missing(bkg.sequences)) {
    if (progress && !BP) cat(" > Shuffling input sequences ...")
    else if ((progress && BP) || verbose > 0) cat(" > Shuffling input sequences\n")
    bkg.sequences <- shuffle_sequences(sequences, shuffle.k, shuffle.method,
                                       shuffle.leftovers, progress = progress,
                                       BP = BP)
  } 

  if (!is.list(motifs)) motifs <- list(motifs)
  motcount <- length(motifs)

  if (use.freq == 1) {
    score.mats <- lapply(motifs, function(x) x@motif)
  } else {
    score.mats <- lapply(motifs,
                         function(x) x@multifreq[[as.character(use.freq)]])
    for (i in seq_along(score.mats)) {
      score.mats[[i]] <- apply(score.mats[[i]], 2, ppm_to_pwmC,
                               nsites = motifs[[i]]@nsites,
                               pseudocount = motifs[[i]]@pseudocount)
    }
  }

  if (threshold.type == "pvalue") {
    if (progress && !BP)
      cat(" > Converting P-values to logodds thresholds ...")
    else if ((progress && BP) || verbose > 0)
      cat(" > Converting P-values to logodds thresholds\n")
    threshold <- motif_pvalue(motifs, pvalue = threshold, use.freq = use.freq,
                              k = 6, progress = progress, BP = BP)
    max.scores <- vapply(score.mats, function(x) sum(apply(x, 2, max)), numeric(1))
    min.scores <- vapply(score.mats, function(x) sum(apply(x, 2, min)), numeric(1))
    for (i in seq_along(threshold)) {
      if (threshold[i] > max.scores[i]) threshold[i] <- max.scores[i]
    }
    if (verbose > 3) {
      mot.names <- vapply(motifs, function(x) x@name, character(1))
      for (i in seq_along(threshold)) {
        cat("   > Motif ", mot.names[i], ": max.score = ", max.scores[i],
            ", threshold.score = ", threshold[i], "\n", sep = "")
      }
    }
    threshold <- (threshold + abs(min.scores)) / (abs(max.scores) + abs(min.scores))
    threshold.type <- "logodds"
  }

  res.all <- enrich_mots2(motifs, sequences, bkg.sequences, threshold,
                          verbose, RC, use.freq, positional.test,
                          search.mode, threshold.type, motcount,
                          return.scan.results, progress = progress, BP = BP)

  if (return.scan.results) {
    res.scan <- res.all[2:3]
    res.all <- res.all[[1]]
    if (nrow(res.scan[[1]]) == 0) res.scan[[1]] <- NULL
    if (nrow(res.scan[[2]]) == 0) res.scan[[2]] <- NULL
  }

  if (nrow(res.all) < 1 || is.null(res.all)) {
    message(" ! No enriched motifs")
    if (return.scan.results) {
      res.all <- c(list(NULL), res.scan)
      names(res.all) <- c("enrichment.report", "input.scan", "bkg.scan")
      return(res.all)
    } else return(invisible(NULL))
  } 

  res.all$Qval.hits <- p.adjust(res.all$Pval.hits, method = qval.method)
  res.all$Qval.pos <- p.adjust(res.all$Pval.pos, method = qval.method)
  res.all$Eval.hits <- res.all$Qval.hits * length(motifs) * 2
  res.all$Eval.pos <- res.all$Qval.pos * length(motifs) * 2

  if (search.mode == "hits") {
    res.all <- res.all[order(res.all$Qval.hits), ]
    res.all <- res.all[res.all$Pval.hits < max.p, ]
    res.all <- res.all[res.all$Qval.hits < max.q, ]
    res.all <- res.all[res.all$Eval.hits < max.e, ]
    res.all <- res.all[, !colnames(res.all) %in% c("seq.pos.mean", "seq.pos.sd",
                                                   "bkg.pos.mean", "bkg.pos.sd",
                                                   "Pval.pos", "Qval.pos",
                                                   "Eval.pos")]
  } else if (search.mode == "positional") {
    res.all <- res.all[order(res.all$Qval.pos), ]
    res.all <- res.all[res.all$Pval.pos < max.p, ]
    res.all <- res.all[res.all$Qval.pos < max.q, ]
    res.all <- res.all[res.all$Eval.pos < max.e, ]
    res.all <- res.all[, !colnames(res.all) %in% c("Pval.hits", "Qval.hits",
                                                    "Eval.hits")]
  } else if (search.mode == "both") {
    res.all <- res.all[order(res.all$Qval.hits), ]
    res.all <- res.all[res.all$Pval.hits < max.p | res.all$Pval.pos < max.p, ]
    res.all <- res.all[res.all$Qval.hits < max.q | res.all$Qval.pos < max.q, ]
    res.all <- res.all[res.all$Eval.hits < max.e | res.all$Eval.pos < max.e, ]
  } else stop("unknown 'search.mode'")


  res.all <- res.all[!is.na(res.all$motif), ]
  rownames(res.all) <- NULL

  if (nrow(res.all) < 1 || is.null(res.all)) {
    message(" ! No enriched motifs")
    if (return.scan.results) {
      res.all <- c(list(NULL), res.scan)
      names(res.all) <- c("enrichment.report", "input.scan", "bkg.scan")
      return(res.all)
    } else return(invisible(NULL))
  } 

  if (return.scan.results) {
    res.all <- c(list(res.all), res.scan)
    names(res.all) <- c("enrichment.report", "input.scan", "bkg.scan")
  }

  res.all

}

enrich_mots2 <- function(motifs, sequences, bkg.sequences, threshold,
                         verbose, RC, use.freq, positional.test,
                         search.mode, threshold.type, motcount,
                         return.scan.results, progress, BP) {

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
  bkg.names <- names(bkg.sequences)
  if (is.null(seq.names)) bkg.names <- seq_len(length(bkg.sequences))

  seq.widths <- width(sequences)
  bkg.widths <- width(bkg.sequences)

  if (return.scan.results) {

    if (verbose > 0) cat(" > Scanning input sequences\n")
    results <- scan_sequences(motifs, sequences, threshold, threshold.type,
                              RC, use.freq, progress = progress,
                              verbose = verbose - 1, BP = BP)

    if (verbose > 0) cat(" > Scanning background sequences\n")
    results.bkg <- scan_sequences(motifs, bkg.sequences, threshold,
                                  threshold.type, RC, use.freq,
                                  verbose = verbose - 1, progress = progress,
                                  BP = BP)

    results2 <- split_by_motif_enrich(motifs, results)
    results.bkg2 <- split_by_motif_enrich(motifs, results.bkg)

  } else {

    seq.alph <- seqtype(sequences)
    bkg.alph <- seqtype(bkg.sequences)
    if (seq.alph != bkg.alph)
      stop("'sequences' and 'bkg.sequences' alphabets do not match")
    mot.alph <- unique(vapply(motifs, function(x) x@alphabet, character(1)))
    if (length(mot.alph) > 1)
      stop("Not all motifs share the same alphabet")
    if (seq.alph != mot.alph && !mot.alph %in% c("DNA", "RNA", "AA"))
      stop("Motif and sequence alphabets differ")
    alph <- switch(mot.alph, "DNA" = DNA_BASES, "RNA" = RNA_BASES,
                   "AA" = AA_STANDARD, safeExplode(mot.alph))

    mot.names <- vapply(motifs, function(x) x@name, character(1))

    inf.check <- vapply(motifs, function(x) any(is.infinite(x@motif)), logical(1))
    if (any(inf.check)) {
      warning("-Inf values found in motifs, normalizing")
      motifs <- mapply(function(x, y) if (x) normalize(y) else y, inf.check, motif,
                       SIMPLIFY = FALSE)
    }
    score.mats <- lapply(motifs, function(x) x@motif)
    score.mats <- lapply(score.mats, numeric_to_integer_matrix)
    max.scores <- vapply(score.mats, function(x) sum(apply(x, 2, max)), integer(1))
    min.scores <- vapply(score.mats, function(x) sum(apply(x, 2, min)), integer(1))
    total.scores <- abs(max.scores) + abs(min.scores)
    threshold <- total.scores * threshold - abs(min.scores)

    if (verbose > 0) cat(" > Scanning input sequences\n")
    results2 <- scan_sequences_slim(score.mats, sequences, threshold,
                                    use.freq, alph, verbose - 1, mot.names,
                                    RC, progress, BP)

    if (verbose > 0) cat(" > Scanning background sequences\n")
    results.bkg2 <- scan_sequences_slim(score.mats, bkg.sequences, threshold,
                                        use.freq, alph, verbose - 1, mot.names,
                                        RC, progress, BP)

  }

  if (length(results2) == 0) {
    return(data.frame())
  }
  if (length(results.bkg2) == 0) results.bkg2 <- data.frame()

  if (progress && !BP) cat(" > Testing motifs for enrichment ...")
  else if ((progress && BP) || verbose > 0) cat(" > Testing motifs for enrichment\n")

  tmp_pb <- progress
  if (verbose > 3) tmp_pb <- FALSE
  results.all <- mapply_(function(x, y, z)
                           enrich_mots2_subworker(x, y, z, seq.widths,
                                                  bkg.widths, search.mode,
                                                  positional.test, sequences,
                                                  RC, bkg.sequences, verbose),
                          results2, results.bkg2, motifs,
                          SIMPLIFY = FALSE, BP = BP, PB = tmp_pb)

  results.all <- do.call(rbind, results.all)

  if (return.scan.results) {
    results.all <- list(results.all, results, results.bkg)
  }

  results.all

}

split_by_motif_enrich <- function(motifs, results) {

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  results.sep <- lapply(mot.names, function(x) results[results$motif == x, ])

  results.sep

}

enrich_mots2_subworker <- function(results, results.bkg, motifs,
                                   seq.widths, bkg.widths, search.mode,
                                   positional.test, sequences, RC,
                                   bkg.sequences, verbose) {

  seq.hits <- results$start
  if (length(seq.hits) == 0 || is.null(seq.hits)) {
    seq.hits.mean <- 0
  } else seq.hits.mean <- mean(seq.hits)

  bkg.hits <- results.bkg$start
  if (length(bkg.hits) == 0 || is.null(bkg.hits)) {
    bkg.hits.mean <- 0
  } else bkg.hits.mean <- mean(bkg.hits)

  if (seq.hits.mean > 0 && bkg.hits.mean == 0) {
    warning(wmsg("Found hits for motif '", motifs@name,
                 "' in target sequences but none in bkg, ",
                 "significance will not be calculated"))
    skip.calc <- TRUE
  } else {
    skip.calc <- FALSE
  }

  if (length(seq.hits) == 0 && length(bkg.hits) == 0) return(NULL)

  seq.total <- (mean(seq.widths) - ncol(motifs@motif) + 1) * length(sequences)
  if (RC) seq.total <- seq.total * 2
  seq.no <- seq.total - length(seq.hits)
  bkg.total <- (mean(bkg.widths) - ncol(motifs@motif) + 1) * length(bkg.sequences)
  if (RC) bkg.total <- bkg.total * 2
  bkg.norm <- seq.total / bkg.total
  bkg.no <- bkg.total - length(bkg.hits)

  results.table <- matrix(c(length(seq.hits), seq.no,
                            length(bkg.hits) * bkg.norm, bkg.no * bkg.norm),
                          nrow = 2, byrow = TRUE)

  if (verbose > 3) cat(" > Testing motif for enrichment:", motifs@name, "\n")

  if (!skip.calc) {
    hits.p <- fisher.test(results.table, alternative = "greater")$p.value
  } else {
    hits.p <- 0
  }

  pos.p <- NA

  if (search.mode %in% c("both", "positional") && length(seq.hits) > 0 && !skip.calc) {
    if (positional.test == "t.test" && length(bkg.hits) > 0) {
      if (verbose > 3) cat("   > Running t.test\n")
      tryCatch({
        pos.p <- t.test(seq.hits / mean(seq.widths),
                        bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("t.test failed"))
    } else if (positional.test == "wilcox.test" && length(bkg.hits) > 0) {
      if (verbose > 3) cat("   > Running wilcox.test\n")
      tryCatch({
        pos.p <- wilcox.test(seq.hits / mean(seq.widths),
                             bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("wilcox.test failed"))
    } else if (positional.test == "chisq.test" && length(bkg.hits) > 0) {
      if (verbose > 3) cat("   > Running chisq.test\n")
      tryCatch({
        pos.p <- chisq.test(seq.hits / mean(seq.widths),
                            bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("chisq.test failed"))
    } else if (positional.test == "shapiro.test") {
      if (verbose > 3) cat("   > Running shapiro.test\n")
      tryCatch({
      pos.p <- shapiro.test(seq.hits)$p.value
      }, error = function(e) warning("shapiro.test failed"))
    }
  } else if (skip.calc) {
    pos.p <- 0
  }

  if (verbose > 3 && search.mode %in% c("hits", "both")) {
    cat("       occurrences p-value:", hits.p, "\n")
  }
  if (verbose > 3 && search.mode %in% c("positional", "both")) {
    cat("       positional bias p-value:", pos.p, "\n")
  }

  results <- data.frame(motif = motifs@name,
                        total.seq.hits = length(seq.hits),
                        num.seqs.hit = length(unique(results$sequence)),
                        num.seqs.total = length(sequences),
                        seq.pos.mean = seq.hits.mean / mean(seq.widths),
                        seq.pos.sd = sd(seq.hits) / mean(seq.widths),
                        total.bkg.hits = length(bkg.hits),
                        num.bkg.hit = length(unique(results.bkg$sequence)),
                        num.bkg.total = length(bkg.sequences),
                        bkg.pos.mean = bkg.hits.mean / mean(bkg.widths),
                        bkg.pos.sd = sd(bkg.hits) / mean(bkg.widths),
                        Pval.hits = hits.p, Pval.pos = pos.p,
                        stringsAsFactors = FALSE)

  results

}
