#' Motif P-value and scoring utility
#'
#' @param motif Single motif.
#' @param score Numeric.
#' @param pvalue Numeric.
#' @param bkg.probs Numeric.
#' @param k Numeric.
#' @param tolerance Numeric.
#'
#' @return Numeric.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_pvalue <- function(motif, score, pvalue, bkg.probs, k = 6,
                         tolerance = 0.75) {

  motif <- convert_motifs(motif)
  if (is.list(motif)) motif <- motif[[1]]
  motif <- convert_type(motif, "PWM")["motif"]

  if (missing(bkg.probs)) bkg.probs <- rep(1 / nrow(motif), nrow(motif))

  if (!missing(score) && missing(pvalue)) {
    return(motif_pval(score.mat = motif, score = score,
                       bkg.probs = bkg.probs, k = k))
  } else if (missing(score) && !missing(pvalue)) {
    return(motif_score(score.mat = motif, pval = pvalue,
                       bkg.probs = bkg.probs, k = k, tolerance = tolerance))
  } else stop("only one of 'score' and 'pvalue' can be used at a time")

}

motif_pval <- function(score.mat, score, bkg.probs, k = 6) {

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

  split.max <- vapply(mot.split, function(x) sum(apply(x, 2, max)), numeric(1))

  split.min <- vector("numeric", length(split.max))
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

  max.score <- sum(apply(score.mat, 2, max)) - 0.0000001
  min.score <- sum(apply(score.mat, 2, min)) + 0.0000001
  if (missing(bkg.probs)) bkg.probs <- rep(1 / nrow(score.mat), nrow(score.mat))

  pv.refine <- motif_pval(score.mat, 0, bkg.probs, k)

  if (pv.refine > pval) score <- 5 else score <- -5
  pv.refine <- motif_pval(score.mat, score, bkg.probs, k)

  if (score < 0) {

    repeat {

      pv.old <- pv.refine
      score.old <- score
    
      if (pv.refine >= pval * tolerance && pv.refine <= pval) break

      pv.factor <- pval / pv.refine 
      score <- score * pv.factor + 0.1

      if (score < min.score) score <- min.score
      pv.refine <- motif_pval(score.mat, score, bkg.probs, k)

      if (pv.old > pv.refine && pv.refine <= pval) break
      if (score == score.old) break

    }
  
  } else if (score > 0) {

    repeat {

      pv.old <- pv.refine
      score.old <- score
    
      if (pv.refine >= pval * tolerance && pv.refine <= pval) break

      pv.factor <- pv.refine / pval
      score <- score * pv.factor + 0.1

      if (score > max.score) score <- max.score
      pv.refine <- motif_pval(score.mat, score, bkg.probs, k)
    
      if (pv.old > pv.refine && pv.refine <= pval) break
      if (score == score.old) break

    }
  
  }

  score

}

.branch_and_bound_kmers <- function(score.mat, min.score) {

  max.scores <- c(rev(cumsum(rev(apply(score.mat, 2, max)))), 0)

  if (min.score > max.scores[1]) stop("input score '", min.score,
                                  "' is higher than max possible score: ",
                                  max.scores[1])

  mot_len <- ncol(score.mat)

  paths <- init_paths_cpp(score.mat, min.score, max.scores[2])
  if (mot_len == 1) return(paths)

  for (i in seq_len(mot_len - 1) + 1) {
    paths <- calc_next_path_cpp(score.mat, paths, min.score, max.scores[i + 1])
    # paths <- do.call(rbind, paths)
  }

  paths

}
