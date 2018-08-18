#' Motif P-value and scoring utility
#'
#' For calculating p-values/logodds scores for any number of motifs.
#'
#' @param motifs See \code{\link{convert_motifs}} for acceptable motif formats.
#' @param score \code{numeric} Get a p-value for a motif from a logodds score.
#' @param pvalue \code{numeric} Get a logodds score for a motif from a 
#'    p-value.
#' @param bkg.probs \code{numeric, list} If supplying individual background
#'    probabilities for each motif, a list. If missing, assumes a uniform
#'    background. Currently does not supported if \code{use.freq > 1}.
#' @param use.freq \code{numeric(1)} By default uses the regular motif matrix;
#'    otherwise uses the corresponding \code{multifreq} matrix.
#' @param k \code{numeric(1)} For speed, scores/p-values can be approximated after 
#'    subsetting the motif every \code{k} columns. If \code{k} is a value
#'    equal or higher to the size of input motif(s), then the calculations
#'    are exact.
#' @param progress_bar \code{logical(1)} If given multiple motifs, show
#'    progress.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return \code{numeric} A vector of scores/p-values.
#'
#' @references
#'    \insertRef{pvalues}{universalmotif}
#'
#' @details
#' Calculating p-values for motifs can be very computationally intensive. This
#' is due to how p-values must be calculated: for a given score, all possible
#' sequences which score equal or higher must be found, and the probability for
#' each of these sequences (based on background probabilities) summed. For a DNA
#' motif of length 10, the number of possible unique sequences is 4^10 = 1,048,576.
#' Finding all possible sequences higher than a given score can be done
#' very efficiently and quickly with a branch-and-bound algorithm, but as the
#' motif length increases this calculation becomes impractical. To get
#' around this, the p-value calculation can be approximated.
#'
#' In order to calculate p-values for longer motifs, this function uses the
#' approximation proposed by \insertCite{pvalues;textual}{universalmotif}, where
#' the motif is subset, p-values calculated for the subsets, and finally
#' combined for a total p-value. The smaller the size of the subsets, the
#' faster the calculation; but also, the bigger the approximation. This can be
#' controlled by setting \code{k}. In fact, for smaller motifs (< 13 positions)
#' calculating exact p-values can be done in reasonable time by setting
#' \code{k = 12}.
#'
#' To calculate a score based on a given p-value, the function simply guesses
#' different scores until it finds one which when used to calculate a p-value,
#' returns a p-value reasonably close to the given p-value.
#'
#' Note that as \code{k} increases (and thus the approximation increases) the
#' resulting p-values increase; meaning the p-values will always be on the
#' conservative side.
#'
#' @examples
#' data(examplemotif)
#'
#' ## p-value/score calculations are performed using the PWM version of the
#' ## motif; these calculations do not work if any -Inf values are present
#' examplemotif["pseudocount"] <- 1
#'
#' ## get a minimum score based on a p-value
#' motif_pvalue(examplemotif, pvalue = 0.001)
#'
#' ## get the probability of a particular sequence hit
#' motif_pvalue(examplemotif, score = 0)
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_pvalue <- function(motifs, score, pvalue, bkg.probs, use.freq = 1, k = 6,
                         progress_bar = FALSE, BPPARAM = SerialParam()) {

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list(score = args$score, pvalue = args$pvalue,
                                     use.freq = args$use.freq, k = args$k),
                                c(0, 0, 1, 1), c(TRUE, TRUE, FALSE, FALSE),
                                "numeric")
  bkg_check <- character()
  if (!missing(bkg.probs)) {
    if (!is.list(bkg.probs) && !is.numeric(bkg.probs)) {
      bkg_check <- paste0(" * Incorrect type for 'bkg.probs': ",
                          "expected 'list' or 'numeric'; got `",
                          class(bkg.probs), "`")
    }
  }
  logi_check <- check_fun_params(list(progress_bar = args$progress_bar),
                                 1, FALSE, "logical")
  s4_check <- check_fun_params(list(BPPARAM = args$BPPARAM), numeric(), FALSE, "S4")
  all_checks <- c(num_check, logi_check, s4_check, bkg_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, "PWM")
  if (!is.list(motifs)) motifs <- list(motifs)
  anyinf <- vapply(motifs, function(x) any(is.infinite(x["motif"])), logical(1))
  if (any(anyinf)) {
    warning("Found -Inf values in motif PWM, adding a pseudocount of 1")
    for (i in which(anyinf)) {
      motifs[[i]] <- convert_type(motifs[[i]], "PPM")
      motifs[[i]]["pseudocount"] <- 1
      motifs[[i]] <- convert_type(motifs[[i]], "PWM")
    }
  }

  if (use.freq == 1) {
    motifs <- lapply(motifs, function(x) x["motif"])
  } else {
    if (!missing(bkg.probs)) warning("'bkg.probs' not supported for 'use.freq' > 1")
    mots <- lapply(motifs, function(x) x["multifreq"][[as.character(use.freq)]])
    motifs <- bpmapply(function(x, y) apply(x, 2, ppm_to_pwmC, 
                                            pseudocount = y["pseudocount"],
                                            nsites = y["nsites"]),
                       mots, motifs, BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  }

  if (missing(bkg.probs) || use.freq > 1) {
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
                    tolerance = 0.75, BPPARAM = BPPARAM)
  } else stop("only one of 'score' and 'pvalue' can be used at a time")

  if (progress_bar) BPPARAM$progressbar <- pb_prev

  out

}

motif_pval <- function(score.mat, score, bkg.probs, k = 6, num2int = TRUE,
                       return_scores = FALSE) {

  if (num2int) {
    if (!return_scores) score <- as.integer(score * 1000)
    score.mat <- matrix(as.integer(score.mat * 1000), nrow = nrow(score.mat))
  }
  total.max <- sum(apply(score.mat, 2, max))
  total.min <- sum(apply(score.mat, 2, min))
  if (return_scores) score <- total.min

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

  if (!return_scores) sum(final.probs) else all.scores

}

fast_fourier_transform <- function(all.scores, pvalue) {

  # https://stats.stackexchange.com/questions/291549/calculate-all-possible-combinations-and-obtain-overall-distribution
  # currently can't find a way to make this work for calculating p-values..

  counts <- vapply(all.scores, length, integer(1))
  n.bins <- 58000
  range.sum <- rowSums(ranges <- vapply(all.scores, range, integer(2)))
  dx <- diff(range.sum) / n.bins

  x.hat <- rep(1, n.bins)

  for (x in all.scores) {
    i <- 1 + round((x - min(x) / dx))
    y <- tabulate(i, nbins = n.bins)
    x.hat <- fft(y / sum(y)) * x.hat
  }

  total.sum.ft <- zapsmall(Re(fft(x.hat, inverse = TRUE))) / n.bins

  i.values <- (1:n.bins - 1/2) * dx + range.sum[1]

  answer <- which(cumsum(total.sum.ft) >= 1 - pvalue)[1]
  i.values[answer]

}

motif_score <- function(score.mat, pval, bkg.probs, k = 6, tolerance = 0.75) {

  score.mat <- matrix(as.integer(score.mat * 1000), nrow = nrow(score.mat))
  max.score <- sum(apply(score.mat, 2, max)) - 1L
  min.score <- sum(apply(score.mat, 2, min)) + 1L
  if (missing(bkg.probs)) bkg.probs <- rep(1 / nrow(score.mat), nrow(score.mat))

  if (ncol(score.mat) <= k) {  # for smaller motifs: exact calculation
    all.scores <- motif_pval(score.mat, min.score, bkg.probs, k,
                             num2int = FALSE, return_scores = TRUE)
    return(quantile(all.scores[[1]], 1 - pval, names = FALSE) / 1000.0)
  }

  pv.refine <- motif_pval(score.mat, 0, bkg.probs, k, num2int = FALSE)

  if (pv.refine > pval) score <- 100L else score <- -100L

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

  max.scores <- c(rev(cumsum(rev(apply(score.mat, 2, max)))), 0L)

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
