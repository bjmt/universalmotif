#' Enrich for input motifs in a set of sequences.
#'
#' Given a set of target and background sequences, test if the input motifs
#' are significantly enriched in the targets sequences relative to the
#' background sequences.
#'
#' @param motifs See \code{\link{convert_motifs}} for acceptable motif formats.
#' @param sequences \code{XStringSet} Alphabet should match motif.
#' @param bkg.sequences \code{XStringSet} Optional; if missing,
#'    \code{\link{shuffle_sequences}} is used to create background sequences from
#'    the input sequences.
#' @param search.mode \code{character(1)} One of \code{c('hits', 'positional', 'both')}.
#'    See details.
#' @param max.p \code{numeric(1)} P-value threshold.
#' @param max.q \code{numeric(1)} Adjusted P-value threshold. This is only useful
#'    if multiple motifs are being enriched for.
#' @param max.e \code{numeric(1)}. The E-value is calculated by multiplying the adjusted
#'    P-value with the number of input motifs \insertCite{meme2}{universalmotif}.
#' @param qval.method \code{character(1)} See \code{\link[stats]{p.adjust}}.
#' @param positional.test \code{character(1)} One of \code{c('t.test', 'wilcox.test',
#'    'chisq.test', 'shapiro.test')}. If using the Shapiro test for
#'    normality, then only the input sequences are tested for positionality;
#'    the background sequences are ignored. See \code{\link[stats]{t.test}},
#'    \code{\link[stats]{wilcox.test}}, \code{\link[stats]{chisq.test}},
#'    \code{\link[stats]{shapiro.test}}.
#' @param threshold \code{numeric(1)} Between 1 and 0. See \code{\link{scan_sequences}}.
#' @param threshold.type \code{character(1)} One of \code{c('logodds', 'pvalue')}. See
#'    \code{\link{scan_sequences}}.
#' @param verbose \code{numeric(1)} 0 for no output, 2 for max verbosity.
#' @param RC \code{logical(1)} Whether to consider the reverse complement of the
#'    sequences. Only available for \code{DNAStringSet}, \code{RNAStringSet}
#'    sequences.
#' @param use.freq \code{numeric(1)} If the \code{multifreq} slot of the motifs are filled,
#'    then they can be used to scan the sequences. See
#'    \code{\link{scan_sequences}}.
#' @param shuffle.k \code{numeric(1)} The k-let size to use when shuffling input
#'    sequences. Only used if no background sequences are input. See
#'    \code{\link{shuffle_sequences}}.
#' @param shuffle.method \code{character(1)} See \code{\link{shuffle_sequences}}.
#' @param shuffle.leftovers \code{character(1)} See \code{\link{shuffle_sequences}}.
#' @param progress_bar \code{logical(1)} Show progress bar from
#'    \code{\link{scans_sequences}}.
#' @param return.scan.results \code{logical(1)} Return output from
#'    \code{\link{scan_sequences}}.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return \code{data.frame} Motif enrichment results. The resulting 
#'    \code{data.frame} contains the following columns:
#'    \itemize{
#'       \item \code{motif} Motif name.
#'       \item \code{total.seq.hits} Total number of matches accross all target
#'          sequences.
#'       \item \code{num.seqs.hits} Number of target sequences which contain matches.
#'       \item \code{num.seqs.total} Number of target sequences.
#'       \item \code{total.bkg.hits} Total number of matches accross all background
#'          sequences.
#'       \item \code{num.bkg.hits} Number of background sequences which contain
#'          matches.
#'       \item \code{num.bkg.total} Number of background sequences.
#'       \item \code{Pval.hits} P-value of enrichment. Only shown if
#'          \code{search.mode = c('hits', 'both')}.
#'       \item \code{Qval.hits} Q-val of enrichment. Only shown if
#'          \code{search.mode = c('hits', 'both')}.
#'       \item \code{Eval.hits} E-val of enrichment. Only shown if
#'          \code{search.mode = c('hits', 'both')}.
#'       \item \code{Pval.pos} P-value of positional comparison. Only
#'          shown if \code{search.mode = c('positional', 'both'}).
#'       \item \code{Qval.pos} Q-value of positional comparison. Only
#'          shown if \code{search.mode = c('positional', 'both'}).
#'       \item \code{Eval.pos} E-value of positional comparison. Only
#'          shown if \code{search.mode = c('positional', 'both'}).
#'    }
#'
#' @details
#' To find enriched motifs, \code{\link{scan_sequences}} is run on both 
#' target and background sequences. If \code{search.mode = 'hits'},
#' \code{\link[stats]{fisher.test}} is run to test for enrichment. If
#' \code{search.mode = 'positional'}, then the test as set by
#' \code{positional.test} is run to check for positional differences
#' between target and background sequences. However if
#' \code{positional.test = 'shapiro.test'}, then only target sequence
#' hits are considered.
#'
#' @examples
#' target.sequences <- create_sequences(monofreqs = c(0.7, 0.1, 0.1, 0.1))
#' bkg.sequences <- create_sequences()
#' motif <- create_motif(bkg = c(0.7, 0.1, 0.1, 0.1))
#' enrich_motifs(motif, target.sequences, bkg.sequences)
#'
#' @references
#'    \insertRef{meme2}{universalmotif}
#'
#' @author Benjamin Tremblay \email{b2tremblay@@uwaterloo.ca}
#' @seealso \code{\link{scan_sequences}}, \code{\link{shuffle_sequences}},
#'    \code{\link{add_multifreq}}, \code{\link{motif_pvalue}}
#' @export
enrich_motifs <- function(motifs, sequences, bkg.sequences, search.mode = "hits",
                          max.p = 10e-6, max.q = 10e-6, max.e = 10e-4,
                          qval.method = "fdr",
                          positional.test = "t.test", threshold = 0.001,
                          threshold.type = "pvalue",
                          verbose = 1, RC = FALSE, use.freq = 1,
                          shuffle.k = 2, shuffle.method = "linear",
                          shuffle.leftovers = "asis",
                          progress_bar = FALSE,
                          return.scan.results = FALSE,
                          BPPARAM = SerialParam()) {

  v1 <- FALSE; v2 <- FALSE
  if (verbose == 1) v1 <- TRUE else if (verbose == 2) v2 <- TRUE

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)

  if (missing(bkg.sequences)) {
    if (verbose > 0) cat(" > Shuffling input sequences\n")
    bkg.sequences <- shuffle_sequences(sequences, shuffle.k, shuffle.method,
                                       shuffle.leftovers, BPPARAM)
  } 

  if (!is.list(motifs)) motifs <- list(motifs)
  motcount <- length(motifs)

  if (threshold.type == "pvalue") {
    if (verbose > 0) cat(" > Converting P-values to logodds thresholds\n")
    threshold <- motif_pvalue(motifs, pvalue = threshold, use.freq = use.freq,
                              progress_bar = progress_bar, BPPARAM = BPPARAM,
                              k = 5)
    if (use.freq == 1) {
      max.scores <- vapply(motifs, function(x) sum(apply(x["motif"], 2, max)),
                           numeric(1))
    } else {
      max.scores <- vapply(motifs,
            function(x) sum(apply(x["multifreq"][[as.character(use.freq)]])),
                           numeric(1))
    }
    threshold <- threshold / max.scores
    threshold.type <- "logodds"
  }

  res.all <- .enrich_mots2(motifs, sequences, bkg.sequences, threshold,
                           verbose, RC, use.freq, positional.test,
                           BPPARAM, search.mode, threshold.type,
                           motcount, v1, v2, progress_bar,
                           return.scan.results)

  if (return.scan.results) {
    res.scan <- res.all[2:3]
    res.all <- res.all[[1]]
    if (nrow(res.scan[[1]]) == 0) res.scan[[1]] <- NULL
    if (nrow(res.scan[[2]]) == 0) res.scan[[2]] <- NULL
  }

  if (nrow(res.all) < 1) {
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

  if (nrow(res.all) < 1) {
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

.enrich_mots2 <- function(motifs, sequences, bkg.sequences, threshold,
                          verbose, RC, use.freq, positional.test,
                          BPPARAM, search.mode, threshold.type,
                          motcount, v1, v2, progress_bar,
                          return.scan.results) {

  if (v1 || v2) cat(" > Scanning input sequences\n")
  results <- scan_sequences(motifs, sequences, threshold, threshold.type,
                            RC, use.freq, progress_bar = progress_bar,
                            BPPARAM = BPPARAM, verbose = v2)

  if (v1 || v2) cat(" > Scanning background sequences\n")
  results.bkg <- scan_sequences(motifs, bkg.sequences, threshold,
                                threshold.type, RC, use.freq,
                                progress_bar = progress_bar,
                                verbose = v2, BPPARAM = BPPARAM)

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
  bkg.names <- names(bkg.sequences)
  if (is.null(seq.names)) bkg.names <- seq_len(length(bkg.sequences))

  seq.widths <- width(sequences)
  bkg.widths <- width(bkg.sequences)

  results2 <- .split_by_motif(motifs, results)
  results.bkg2 <- .split_by_motif(motifs, results.bkg)
  if (length(results2) == 0) stop("no matches to motifs found")
  if (length(results.bkg2) == 0) results.bkg2 <- data.frame()

  if (v1) cat(" > Testing motifs for enrichment\n")

  results.all <- bpmapply(function(x, y, z)
                          .enrich_mots2_subworker(x, y, z, seq.widths,
                                                  bkg.widths, search.mode,
                                                  positional.test, sequences,
                                                  RC, bkg.sequences, v2, v1),
                          results2, results.bkg2, motifs,
                          SIMPLIFY = FALSE, BPPARAM = BPPARAM)

  results.all <- do.call(rbind, results.all)

  if (return.scan.results) {
    results.all <- list(results.all, results, results.bkg)
  }

  results.all

}

.split_by_motif <- function(motifs, results) {

  mot.names <- vapply(motifs, function(x) x["name"], character(1))
  results.sep <- list(length = length(mot.names))

  for (i in seq_along(mot.names)) {
    results.sep[[i]] <- results[results$motif == mot.names[i], ]
  }

  results.sep

}

.enrich_mots2_subworker <- function(results, results.bkg, motifs,
                                    seq.widths, bkg.widths, search.mode,
                                    positional.test, sequences, RC,
                                    bkg.sequences, v2, v1) {

  seq.hits <- results$start
  if (length(seq.hits) == 0) {
    seq.hits.mean <- 0
  } else seq.hits.mean <- mean(seq.hits)

  bkg.hits <- results.bkg$start
  if (length(bkg.hits) == 0) {
    bkg.hits.mean <- 0
  } else bkg.hits.mean <- mean(bkg.hits)

  if (length(seq.hits) == 0 && length(bkg.hits) == 0) return(NULL)

  seq.total <- (mean(seq.widths) - ncol(motifs["motif"]) + 1) *
               length(sequences)
  if (RC) seq.total <- seq.total * 2
  seq.no <- seq.total - length(seq.hits)
  bkg.total <- (mean(bkg.widths) - ncol(motifs["motif"]) + 1) *
               length(bkg.sequences)
  if (RC) bkg.total <- bkg.total * 2
  bkg.norm <- seq.total / bkg.total
  bkg.no <- bkg.total - length(bkg.hits)

  results.table <- matrix(c(length(seq.hits), seq.no,
                            length(bkg.hits) * bkg.norm, bkg.no * bkg.norm),
                          nrow = 2, byrow = TRUE)

  if (v2) cat(" > Testing motif for enrichment:", motifs["name"], "\n")

  hits.p <- fisher.test(results.table, alternative = "greater")$p.value

  pos.p <- NA

  if (search.mode %in% c("both", "positional") && length(seq.hits) > 0) {
    if (positional.test == "t.test" && length(bkg.hits) > 0) {
      tryCatch({
        pos.p <- t.test(seq.hits / mean(seq.widths),
                        bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("t.test failed"))
    } else if (positional.test == "wilcox.test" && length(bkg.hits) > 0) {
      tryCatch({
        pos.p <- wilcox.test(seq.hits / mean(seq.widths),
                             bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("wilcox.test failed"))
    } else if (positional.test == "chisq.test" && length(bkg.hits) > 0) {
      tryCatch({
        pos.p <- chisq.test(seq.hits / mean(seq.widths),
                            bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("chisq.test failed"))
    } else if (positional.test == "shapiro.test") {
      tryCatch({
      pos.p <- shapiro.test(seq.hits)$p.value
      }, error = function(e) warning("shapiro.test failed"))
    } 
  }

  if (v2 && search.mode %in% c("hits", "both")) {
    cat("       occurrences p-value:", hits.p, "\n")
  }
  if (v2 && search.mode %in% c("positional", "both")) {
    cat("       positional bias p-value:", pos.p, "\n")
  } 

  results <- data.frame(motif = motifs["name"],
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
