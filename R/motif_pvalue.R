#' Motif P-value and scoring utility
#'
#' For calculating p-values/scores for any number of motifs. Vectorized;
#' arguments are recycled.
#'
#' @param motif Single motif or list of motifs.
#' @param score Numeric. 
#' @param pvalue Numeric.
#' @param bkg.probs Numeric. If supplying individual background probabilities
#'    for each motif, a list.
#' @param use.freq Numeric.
#' @param k Numeric.
#' @param progress_bar Logical.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Numeric.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_pvalue <- function(motifs, score, pvalue, bkg.probs, use.freq = 1, k = 6,
                         progress_bar = FALSE, BPPARAM = SerialParam()) {

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, "PWM")

  if (!is.list(motifs)) motifs <- list(motifs)
  if (use.freq == 1) {
    motifs <- lapply(motifs, function(x) x["motif"])
  } else {
    motifs <- lapply(motifs, function(x) x["multifreq"][[as.character(use.freq)]])
  }

  if (missing(bkg.probs)) {
    bkg.probs <- lapply(motifs, function(x) rep( 1 / nrow(x), nrow(x)))
  } else if (!is.list(bkg.probs)) bkg.probs <- list(bkg.probs)

  if (progress_bar) {
    pb_prev <- BPPARAM$progressbar
    BPPARAM$progressbar <- TRUE
  }

  if (!missing(score) && missing(pvalue)) {
    out <- bpmapply(motif_pval, motifs, score, bkg.probs, k, BPPARAM = BPPARAM)
  } else if (missing(score) && !missing(pvalue)) {
    out <- bpmapply(motif_score, motifs, pvalue, bkg.probs, k,
                    BPPARAM = BPPARAM)
  } else stop("only one of 'score' and 'pvalue' can be used at a time")

  if (progress_bar) BPPARAM$progressbar <- pb_prev

  out

}

motif_pval <- function(score.mat, score, bkg.probs, k = 6, num2int = TRUE) {

  if (num2int) {
    score <- as.integer(score * 1000)
    score.mat <- matrix(as.integer(score.mat * 1000), nrow = nrow(score.mat))
  }
  total.max <- sum(apply(score.mat, 2, max))
  total.min <- sum(apply(score.mat, 2, min))

  if (missing(bkg.probs)) bkg.probs <- rep(1 / nrow(score.mat), nrow(score.mat))

  score.mat <- score.mat[, order(apply(score.mat, 2, max), decreasing = TRUE)]
  alph.sort <- apply(score.mat, 2, order, decreasing = TRUE)
  for (i in seq_len(ncol(score.mat))) {
    score.mat[, i] <- score.mat[alph.sort[, i], i]
  }
  mot.len <- ncol(score.mat)
  alph.len <- nrow(score.mat)

  if (mot.len > k) {

    times.tosplit <- mot.len %/% k
    leftover.split <- mot.len %% k
    mot.split <- vector("list", times.tosplit + ifelse(leftover.split > 0, 1, 0))
    alph.sort.split <- mot.split
    mot.split[[1]] <- score.mat[, seq_len(k)]
    alph.sort.split[[1]] <- alph.sort[, seq_len(k)]

    if (times.tosplit > 1) {
      for (i in seq_len(times.tosplit - 1)) {
        mot.split[[i + 1]] <- score.mat[, (i * k + 1):(i * k + k)]
        alph.sort.split[[i + 1]] <- alph.sort[, (i * k + 1):(i * k + k)]
      }
    }

    if (leftover.split > 0) {
      mot.split[[length(mot.split)]] <- score.mat[, (mot.len - leftover.split + 1):
                                                     mot.len]
      alph.sort.split[[length(mot.split)]] <- alph.sort[, (mot.len - leftover.split + 1):
                                                        mot.len]
      if (!is.matrix(mot.split[[length(mot.split)]])) {
        mot.split[[length(mot.split)]] <- matrix(mot.split[[length(mot.split)]])
        alph.sort.split[[length(mot.split)]] <- matrix(alph.sort.split[[length(mot.split)]])
      }

    }

  } else {
    mot.split <- list(score.mat)
    alph.sort.split <- list(alph.sort)
  }

  split.max <- vapply(mot.split, function(x) sum(apply(x, 2, max)), integer(1))

  split.min <- vector("integer", length(split.max))
  for (i in seq_along(split.max)) {
    split.min[i] <- score - sum(split.max[-i])
  }

  all.paths <- vector("list", length(mot.split))
  for (i in seq_along(all.paths)) {
    all.paths[[i]] <- .branch_and_bound_kmers(mot.split[[i]], split.min[i])
    # all.paths[[i]] <- branch_and_bound_cpp(mot.split[[i]], split.min[i])
  }

  all.scores <- vector("list", length(mot.split))
  for (i in seq_along(all.scores)) {
    all.scores[[i]] <- calc_scores_cpp(all.paths[[i]], mot.split[[i]])
  }

  all.probs <- vector("list", length(mot.split))
  for (i in seq_along(all.probs)) {
    all.probs[[i]] <- kmer_mat_to_probs_k1_cpp(all.paths[[i]], bkg.probs,
                                               alph.sort.split[[i]])
  }

  max.scores5 <- vapply(all.scores, max, numeric(1))

  if (length(mot.split) > 2) {
  
    times.toloop <- length(mot.split) - 2

    for (i in seq_len(times.toloop)) {
    
      for (j in seq_along(all.probs[[i + 1]])) {
        all.probs[[i + 1]][j] <- all.probs[[i + 1]][j] *
          sum(all.probs[[i + 2]][all.scores[[i + 2]] >
              score - all.scores[[i + 1]][i] - sum(max.scores5[seq_len(i)])])
      }
    
    }
  
  }

  if (length(mot.split) > 1) {
    final.probs <- vector("numeric", length(all.scores[[1]]))
    for (i in seq_along(final.probs)) {
      final.probs[i] <- all.probs[[1]][i] *
        sum(all.probs[[2]][all.scores[[2]] > score - all.scores[[1]][i]])
    }
    # final.probs <- calc_final_probs_cpp(all.probs, all.scores, score)
  }

  if (length(mot.split) == 1) {
    final.probs <- all.probs[[1]]
  }

  sum(final.probs)

}

motif_score <- function(score.mat, pval, bkg.probs, k = 6, tolerance = 0.75) {

  score.mat <- matrix(as.integer(score.mat * 1000), nrow = nrow(score.mat))
  max.score <- sum(apply(score.mat, 2, max)) - 1
  min.score <- sum(apply(score.mat, 2, min)) + 1
  if (missing(bkg.probs)) bkg.probs <- rep(1 / nrow(score.mat), nrow(score.mat))

  pv.refine <- motif_pval(score.mat, 0, bkg.probs, k, num2int = FALSE)

  if (pv.refine > pval) score <- 100 else score <- -100

  if (score < 0) {

    repeat {

      pv.old <- pv.refine
      score.old <- score
    
      if (pv.refine >= pval * tolerance && pv.refine <= pval) break

      pv.factor <- pval / pv.refine 
      score <- as.integer(score * pv.factor + 10)

      if (score < min.score) score <- min.score
      pv.refine <- motif_pval(score.mat, score, bkg.probs, k, num2int = FALSE)

      if (pv.old > pv.refine && pv.refine <= pval) break
      if (score == score.old) break

    }
  
  } else if (score > 0) {

    repeat {

      pv.old <- pv.refine
      score.old <- score
    
      if (pv.refine >= pval * tolerance && pv.refine <= pval) break

      pv.factor <- pv.refine / pval
      score <- as.integer(score * pv.factor + 10)

      if (score > max.score) score <- max.score
      pv.refine <- motif_pval(score.mat, score, bkg.probs, k, num2int = FALSE)
    
      if (pv.old > pv.refine && pv.refine <= pval) break
      if (score == score.old) break

    }
  
  }

  score / 1000.0

}

.branch_and_bound_kmers <- function(score.mat, min.score) {

  max.scores <- c(rev(cumsum(rev(apply(score.mat, 2, max)))), 0)

  if (min.score > max.scores[1]) stop("input score '", min.score / 1000.0,
                                  "' is higher than max possible score: '",
                                  max.scores[1] / 1000.0, "'")

  mot_len <- ncol(score.mat)

  paths <- init_paths_cpp(score.mat, min.score, max.scores[2])
  if (mot_len == 1) return(paths)

  for (i in seq_len(mot_len - 1) + 1) {
    paths <- calc_next_path_cpp(score.mat, paths, min.score, max.scores[i + 1])
    # paths <- do.call(rbind, paths)
  }

  paths

}
