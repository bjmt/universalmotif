#' Enrich for input motifs in a set of sequences.
#'
#' @param motifs Motifs.
#' @param sequences DNAStringSet.
#' @param bkg.sequences DNAStringSet.
#' @param mode Character.
#' @param max.p Numeric.
#' @param max.q Numeric.
#' @param qval.method Numeric.
#' @param positional.test Character.
#' @param threshold Numeric.
#' @param RC Logical.
#' @param HMMorder Numeric.
#' @param shuffle.k Numeric.
#' @param shuffle.method Character.
#' @param shuffle.leftovers Character.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return A list of results.
#'
#' @author Benjamin Tremblay \email{b2tremblay@@uwaterloo.ca}
#' @export
enrich_motifs <- function(motifs, sequences, bkg.sequences, mode = "hits",
                          max.p = 0.001, max.q = 0.001, qval.method = "fdr",
                          positional.test = "t.test", threshold = 0.6,
                          verbose = TRUE, RC = FALSE, HMMorder = 0,
                          shuffle.k = 1, shuffle.method = "linear",
                          shuffle.leftovers = "asis", plot = FALSE,
                          plot.top = 5, BPPARAM = bpparam(), ...) {

  if (missing(bkg.sequences)) {
    if (verbose) cat("Shuffling input sequences\n")
    bkg.sequences <- shuffle_sequences(sequences, shuffle.k, shuffle.method,
                                       shuffle.leftovers, BPPARAM)
  }

  if (!is.list(motifs)) motifs <- list(motifs)

  res.all <- lapply(motifs, function(x) .enrich_mots(x, sequences,
                                         bkg.sequences,
                                         threshold, verbose, RC, HMMorder,
                                         positional.test, BPPARAM))

  res.all <- do.call(rbind, res.all)

  if (nrow(res.all) < 1) {
    message("no enriched motifs")
    return(NULL)
  } 

  res.all$Qval.hits <- p.adjust(res.all$Pval.hits, method = qval.method)
  res.all$Qval.pos <- p.adjust(res.all$Pval.pos, method = qval.method)

  if (mode == "hits") {
    res.all <- res.all[order(res.all$Qval.hits), ]
    res.all <- res.all[res.all$Pval.hits < max.p, ]
    res.all <- res.all[res.all$Qval.hits < max.q, ]
  } else if (mode == "positional") {
    res.all <- res.all[order(res.all$Qval.pos), ]
    res.all <- res.all[res.all$Pval.pos < max.p, ]
    res.all <- res.all[res.all$Qval.pos < max.q, ]
  } else if (mode == "both") {
    res.all <- res.all[order(res.all$Qval.hits), ]
    res.all <- res.all[res.all$Pval.hits < max.p |
                       res.all$Pval.pos < max.p, ]
    res.all <- res.all[res.all$Qval.hits < max.q |
                       res.all$Qval.pos < max.q, ]
  } else stop("unknown 'mode'")

  rownames(res.all) <- NULL

  if (nrow(res.all) < 1) {
    message("no enriched motifs")
    return(NULL)
  } 

  return(res.all)

}

.enrich_mots <- function(motifs, sequences, bkg.sequences, threshold, verbose,
                         RC, HMMorder, positional.test, BPPARAM = BPPARAM) {

 if (verbose) cat("Enriching sequences for motif:", motifs["name"], "\n")

  if (verbose) cat("  scanning input sequences\n")
  results <- suppressMessages(scan_sequences(motifs, sequences, threshold, RC,
                                             HMMorder, BPPARAM))
  if (verbose) cat("  scanning background sequences\n")
  results.bkg <- suppressMessages(scan_sequences(motifs, bkg.sequences, threshold,
                                                 RC, HMMorder, BPPARAM))

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
  bkg.names <- names(bkg.sequences)
  if (is.null(bkg.names)) bkg.names <- seq_len(length(bkg.sequences))

  seq.widths <- width(sequences)
  bkg.widths <- width(bkg.sequences)

  seq.hits <- results$start
  bkg.hits <- results.bkg$start

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

  if (verbose) cat("  testing for enrichment\n")

  hits.p <- fisher.test(results.table, alternative = "greater")$p.value

  if (positional.test == "t.test") {
    pos.p <- t.test(seq.hits / mean(seq.widths),
                    bkg.hits / mean(bkg.widths))$p.value
  } else if (positional.test == "wilcox.test") {
    pos.p <- wilcox.test(seq.hits / mean(seq.widths),
                         bkg.hits / mean(bkg.widths))$p.value
  } else if (positional.test == "chisq.test") {
    pos.p <- chisq.test(seq.hits / mean(seq.widths),
                        bkg.hits / mean(bkg.widths))$p.value
  } else if (positional.test == "shapiro.test") {
    pos.p <- shapiro.test(seq.hits)$p.value
  } else stop("unknown 'positional.test'")


  if (verbose && mode %in% c("hits", "both")) {
    cat("    occurrences p-value:", hits.p, "\n")
  }
  if (verbose && mode %in% c("positional", "both")) {
    cat("    positional bias p-value:", pos.p, "\n")
  } 

  results <- data.frame(motif = motifs["name"], sequence.hits = length(seq.hits),
                        num.sequences = length(sequences),
                        seq.pos.mean = mean(seq.hits) / mean(seq.widths),
                        seq.pos.sd = sd(seq.hits) / mean(seq.widths),
                        bkg.hits = length(bkg.hits),
                        num.bkg = length(bkg.sequences),
                        bkg.pos.mean = mean(bkg.hits) / mean(bkg.widths),
                        bkg.pos.sd = sd(bkg.hits) / mean(bkg.widths),
                        Pval.hits = fisher.p, Pval.pos = pos.p)

  results

}
