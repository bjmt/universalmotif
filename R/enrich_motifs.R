#' Enrich for input motifs in a set of sequences.
#'
#' DNA/RNA only.
#'
#' @param motifs Motifs.
#' @param sequences DNAStringSet.
#' @param bkg.sequences DNAStringSet.
#' @param search.mode Character.
#' @param max.p Numeric.
#' @param max.q Numeric.
#' @param qval.method Numeric.
#' @param positional.test Character.
#' @param threshold Numeric.
#' @param threshold.type Character. One of 'logodds' and 'pvalue'.
#' @param verbose Logical.
#' @param RC Logical.
#' @param use.freq Numeric.
#' @param shuffle.k Numeric.
#' @param shuffle.method Character.
#' @param shuffle.leftovers Character.
#' @param BPPARAM See \code{\link[BiocParallel]{SerialParam}}.
#'
#' @return A list of results.
#'
#' @author Benjamin Tremblay \email{b2tremblay@@uwaterloo.ca}
#' @export
enrich_motifs <- function(motifs, sequences, bkg.sequences, search.mode = "hits",
                          max.p = 0.001, max.q = 0.001, qval.method = "fdr",
                          positional.test = "t.test", threshold = 0.6,
                          threshold.type = "logodds",
                          verbose = TRUE, RC = TRUE, use.freq = 1,
                          shuffle.k = 1, shuffle.method = "linear",
                          shuffle.leftovers = "asis", BPPARAM = SerialParam()) {

  # sequences <- DNAStringSet(sequences)

  if (missing(bkg.sequences)) {
    if (verbose) cat("Shuffling input sequences\n")
    bkg.sequences <- shuffle_sequences(sequences, shuffle.k, shuffle.method,
                                       shuffle.leftovers, BPPARAM)
  } #else bkg.sequences <- DNAStringSet(bkg.sequences)

  if (!is.list(motifs)) motifs <- list(motifs)

  res.all <- lapply(motifs, function(x) .enrich_mots(x, sequences,
                                         bkg.sequences,
                                         threshold, verbose, RC, use.freq,
                                         positional.test, BPPARAM,
                                         search.mode, threshold.type))

  res.all <- do.call(rbind, res.all)

  if (nrow(res.all) < 1) {
    message("no enriched motifs")
    return(NULL)
  } 

  res.all$Qval.hits <- p.adjust(res.all$Pval.hits, method = qval.method)
  res.all$Qval.pos <- p.adjust(res.all$Pval.pos, method = qval.method)

  if (search.mode == "hits") {
    res.all <- res.all[order(res.all$Qval.hits), ]
    res.all <- res.all[res.all$Pval.hits < max.p, ]
    res.all <- res.all[res.all$Qval.hits < max.q, ]
    res.all <- res.all[, -c(4, 5, 8, 9, 11, 13)]
  } else if (search.mode == "positional") {
    res.all <- res.all[order(res.all$Qval.pos), ]
    res.all <- res.all[res.all$Pval.pos < max.p, ]
    res.all <- res.all[res.all$Qval.pos < max.q, ]
    res.all <- res.all[, -c(2, 3, 6, 7, 10, 12)]
  } else if (search.mode == "both") {
    res.all <- res.all[order(res.all$Qval.hits), ]
    res.all <- res.all[res.all$Pval.hits < max.p |
                       res.all$Pval.pos < max.p, ]
    res.all <- res.all[res.all$Qval.hits < max.q |
                       res.all$Qval.pos < max.q, ]
  } else stop("unknown 'search.mode'")

  rownames(res.all) <- NULL

  if (nrow(res.all) < 1) {
    message("no enriched motifs")
    return(NULL)
  } 

  return(res.all)

}

.enrich_mots <- function(motifs, sequences, bkg.sequences, threshold, verbose,
                         RC, use.freq, positional.test, 
                         BPPARAM = BPPARAM, search.mode, threshold.type) {

  if (verbose) cat("Enriching sequences for motif:", motifs["name"], "\n")

  if (verbose) cat("  scanning input sequences\n")
  results <- suppressMessages(scan_sequences(motifs, sequences, threshold,
                                             threshold.type, RC,
                                             use.freq, BPPARAM))
  if (verbose) cat("  scanning background sequences\n")
  results.bkg <- suppressMessages(scan_sequences(motifs, bkg.sequences, threshold,
                                                 threshold.type, RC, use.freq,
                                                 BPPARAM))

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))
  bkg.names <- names(bkg.sequences)
  if (is.null(bkg.names)) bkg.names <- seq_len(length(bkg.sequences))

  seq.widths <- width(sequences)
  bkg.widths <- width(bkg.sequences)

  seq.hits <- results$start
  if (length(seq.hits) == 0) {
    seq.hits.mean <- 0
  } else seq.hits.mean <- mean(seq.hits)
  bkg.hits <- results.bkg$start
  if (length(bkg.hits) == 0) {
    bkg.hits.mean <- 0
  } else bkg.hits.mean <- mean(bkg.hits)

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

  pos.p <- NA
  if (search.mode %in% c("both", "positional")) {
    if (positional.test == "t.test") {
      tryCatch({
        pos.p <- t.test(seq.hits / mean(seq.widths),
                        bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("t.test failed"))
    } else if (positional.test == "wilcox.test") {
      tryCatch({
        pos.p <- wilcox.test(seq.hits / mean(seq.widths),
                             bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("wilcox.test failed"))
    } else if (positional.test == "chisq.test") {
      tryCatch({
        pos.p <- chisq.test(seq.hits / mean(seq.widths),
                            bkg.hits / mean(bkg.widths))$p.value
      }, error = function(e) warning("chisq.test failed"))
    } else if (positional.test == "shapiro.test") {
      tryCatch({
      pos.p <- shapiro.test(seq.hits)$p.value
      }, error = function(e) warning("shapiro.test failed"))
    } else stop("unknown 'positional.test'")
  }


  if (verbose && search.mode %in% c("hits", "both")) {
    cat("    occurrences p-value:", hits.p, "\n")
  }
  if (verbose && search.mode %in% c("positional", "both")) {
    cat("    positional bias p-value:", pos.p, "\n")
  } 

  results <- data.frame(motif = motifs["name"], sequence.hits = length(seq.hits),
                        num.sequences = length(sequences),
                        seq.pos.mean = seq.hits.mean / mean(seq.widths),
                        seq.pos.sd = sd(seq.hits) / mean(seq.widths),
                        bkg.hits = length(bkg.hits),
                        num.bkg = length(bkg.sequences),
                        bkg.pos.mean = bkg.hits.mean / mean(bkg.widths),
                        bkg.pos.sd = sd(bkg.hits) / mean(bkg.widths),
                        Pval.hits = hits.p, Pval.pos = pos.p)

  results

}
