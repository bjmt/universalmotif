#' Motif P-value and scoring utility
#'
#' For calculating p-values/logodds scores for any number of motifs.
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param score `numeric` Get a p-value for a motif from a logodds score.
#' @param pvalue `numeric` Get a logodds score for a motif from a
#'    p-value.
#' @param bkg.probs `numeric`, `list` If supplying individual background
#'    probabilities for each motif, a list. If missing, retrieves the
#'    background from the motif `bkg` slot. Note that this only influences
#'    calculating p-values from an input score; calculating a score from an
#'    input p-value currently assumes a uniform background.
#' @param use.freq `numeric(1)` By default uses the regular motif matrix;
#'    otherwise uses the corresponding `multifreq` matrix. Max is 3.
#' @param k `numeric(1)` For speed, scores/p-values can be approximated after
#'    subsetting the motif every `k` columns. If `k` is a value
#'    equal or higher to the size of input motif(s), then the calculations
#'    are (nearly) exact.
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [motif_pvalue()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for exceptionally large
#'    jobs. Furthermore, the behaviour of `progress = TRUE` is changed
#'    if `BP = TRUE`; the default \pkg{BiocParallel} progress bar will be
#'    shown (which unfortunately is much less informative).
#'
#' @return `numeric` A vector of scores/p-values.
#'
#' @references
#'    \insertRef{pvalues}{universalmotif}
#'
#' @details
#'
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
#' controlled by setting `k`. In fact, for smaller motifs (< 13 positions)
#' calculating exact p-values can be done individually in reasonable time by
#' setting `k = 12`.
#'
#' To calculate a score based on a given p-value, the means and variances of
#' each motif subsets are combined to estimate the distribution of all
#' possible scores using [stats::qnorm()]:
#'
#' `qnorm(pvalue, mean = sum(subset.means), sd = sqrt(sum(subset.vars)))`
#'
#' For calculating exact scores, [stats::ecdf()] and [stats::quantile()] are
#' used:
#'
#' `quantile(ecdf(scores), probs = pvalue)`
#'
#' It is important to keep in mind that both approximate and exact score
#' calculations assume uniform backgrounds, so do not use this function for
#' motifs with extremely imbalanced backgrounds. To get all possible scores for
#' each subset, [expand.grid()] is used instead of the branch-and-bound
#' algorithm used for calculating p-values. Keep this in mind for determining
#' the best `k` value for motifs with alphabets longer than those of DNA/RNA
#' motifs.
#'
#' @examples
#' data(examplemotif)
#'
#' ## p-value/score calculations are performed using the PWM version of the
#' ## motif; these calculations do not work if any -Inf values are present
#' examplemotif["pseudocount"] <- 1
#' # or
#' examplemotif <- BiocGenerics::normalize(examplemotif)
#'
#' ## get a minimum score based on a p-value
#' motif_pvalue(examplemotif, pvalue = 0.001)
#'
#' ## get the probability of a particular sequence hit
#' motif_pvalue(examplemotif, score = 0)
#'
#' ## the calculations can be performed for multiple motifs
#' motif_pvalue(list(examplemotif, examplemotif), pvalue = c(0.001, 0.0001))
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_pvalue <- function(motifs, score, pvalue, bkg.probs, use.freq = 1,
                         k = 6, progress = ifelse(length(motifs) > 1, TRUE, FALSE),
                         BP = FALSE) {

  # TODO: Number precision is terrible. Can't get P-values < 1e-8.

  # NOTE: The calculated P-value is the chance of getting a certain score at
  #       one position. To get a P-value from scanning a 2000 bp stretch for
  #       example, the P-value is multiplied by the number of possible positions
  #       the motif can find itself at. For example:
  #
  #       R> motif_pvalue(ArabidopsisMotif, 15)
  #       [1] 9.779e-08
  #
  #       For a 2000 bp sequence (and motif length 15):
  #
  #       R> 9.779e-08 * (2000 - 15 + 1)
  #       [1] 0.0001942
  #
  #       This number is the probability of finding the motif with this score
  #       or higher once in a 2000 bp sequence. Let's say it is found three
  #       times with this exact score:
  #
  #       R> 0.0001942^3
  #       [1] 7.325e-12

  # Previously removed from examples section:
  #
  # ## get motif site p-values after using scan_sequences()
  # data(ArabidopsisMotif)
  # data(ArabidopsisPromoters)
  # res <- scan_sequences(ArabidopsisMotif, ArabidopsisPromoters, RC = FALSE,
  #                       progress = FALSE, verbose = 0, threshold = 0,
  #                       threshold.type = "logodds")[1:100, ]
  # res$pvalue <- motif_pvalue(ArabidopsisMotif, score = res$score)
  #

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list(score = args$score, pvalue = args$pvalue,
                                     use.freq = args$use.freq, k = args$k),
                                c(0, 0, 1, 1), c(TRUE, TRUE, FALSE, FALSE),
                                "numeric")
  logi_check <- check_fun_params(list(progress = args$progress, BP = args$BP),
                                 numeric(), logical(), "logical")
  bkg_check <- character()
  if (!missing(bkg.probs)) {
    if (!is.list(bkg.probs) && !is.numeric(bkg.probs)) {
      bkg_check <- paste0(" * Incorrect type for 'bkg.probs': ",
                          "expected 'list' or 'numeric'; got `",
                          class(bkg.probs), "`")
    }
  }
  use.freq_check <- character()
  if (use.freq > 3) {
    use.freq_check <- " * Incorrect 'use.freq': maximum allowed is 3"
  }
  all_checks <- c(num_check, bkg_check, logi_check, use.freq_check)
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

  motifs2 <- motifs

  if (use.freq == 1) {
    motifs <- lapply(motifs, function(x) x["motif"])
  } else {
    mots <- lapply(motifs, function(x) x["multifreq"][[as.character(use.freq)]])
    motifs <- mapply(function(x, y) apply(x, 2, ppm_to_pwmC,
                                            pseudocount = y["pseudocount"],
                                            nsites = y["nsites"]),
                       mots, motifs, SIMPLIFY = FALSE)
  }

  if (!missing(score) && missing(pvalue)) {

    if (!missing(bkg.probs)) {

      if (!is.list(bkg.probs)) bkg.probs <- list(bkg.probs)
      bkg.probs.len <- lapply(bkg.probs, length)
      motif.nrow <- lapply(motifs, nrow)
      alph.len.check <- mapply(function(x, y) x != y^use.freq,
                               bkg.probs.len, motif.nrow, SIMPLIFY = TRUE)
      if (any(alph.len.check))
        stop("length(bkg.probs) must match nrow(motif['motif'])^use.freq")

    } else bkg.probs <- rep(list(NULL), length(motifs))

    bkg.probs <- mapply(motif_pvalue_bkg, motifs2, bkg.probs,
                        MoreArgs = list(use.freq = use.freq),
                        SIMPLIFY = FALSE)

    out <- mapply_(motif_pval, motifs, score, bkg.probs, k, PB = progress,
                   BP = BP)

  } else if (missing(score) && !missing(pvalue)) {

    out <- mapply_(motif_score, motifs, pvalue, k, PB = progress, BP = BP)

  } else stop("only one of 'score' and 'pvalue' can be used at a time")

  out

}

motif_pvalue_bkg <- function(motif, bkg.probs, use.freq) {

  lets1 <- rownames(motif@motif)
  if (use.freq > 1)
    lets2 <- sort(collapse_rows_df(expand.grid(rep(list(lets1), use.freq),
                                               stringsAsFactors = FALSE)))
  if (is.null(bkg.probs)) {

    if (use.freq == 1) out <- motif@bkg[lets1]
    else {
      out <- rep(1 / length(lets2), length(lets2))
      names(out) <- lets2
    }

  } else {

    if (use.freq == 1) out <- bkg.probs[lets1]
    else {
      out <- bkg.probs[lets2]
      if (any(is.na(out))) {
        message(wmsg("Could not find higher order background probabilities from",
                     " motif object, assuming uniform background"))
        out <- rep(1 / length(lets2), length(lets2))
        names(out) <- lets2
      }
    }

  }

  out

}

motif_pval <- function(score.mat, score, bkg.probs, k = 6, num2int = TRUE,
                       return_scores = FALSE) {

  # lets1 <- rownames(score.mat)
  # use.order <- "0"
  # bkg_len <- length(bkg.probs)
  # alph_len <- nrow(score.mat)
  # if (bkg_len > alph_len) {
    # if (bkg_len >= alph_len + alph_len^2) {
      # lets2 <- expand.grid(rep(list(lets1), 2), stringsAsFactors = FALSE)
      # lets2 <- sort(collapse_rows_df(lets2))
      # lets2 <- bkg.probs[lets2]
      # use.order <- "1"
    # }
    # if (bkg_len >= alph_len + alph_len^2 + alph_len^3) {
      # lets3 <- expand.grid(rep(list(lets1), 3), stringsAsFactors = FALSE)
      # lets3 <- sort(collapse_rows_df(lets3))
      # lets3 <- bkg.probs[lets3]
      # use.order <- "2"
    # }
  # }

  if (num2int) {
    if (!return_scores) score <- as.integer(score * 1000)
    score.mat <- matrix(as.integer(score.mat * 1000), nrow = nrow(score.mat))
  }
  total.max <- sum(apply(score.mat, 2, max))
  if (score > total.max) stop(wmsg("input score '", score / 1000,
                                   "' is higher than max possible score: '",
                                   total.max / 1000, "'"))
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
    all.paths[[i]] <-  branch_and_bound_kmers(mot.split[[i]], split.min[i])
    # all.paths[[i]] <- branch_and_bound_cpp(mot.split[[i]], split.min[i])
  }

  all.scores <- vector("list", length(mot.split))
  for (i in seq_along(all.scores)) {
    all.scores[[i]] <- calc_scores_cpp(all.paths[[i]], mot.split[[i]])
  }

  all.probs <- vector("list", length(mot.split))
  # switch(use.order,
#
    # "0" = {
      for (i in seq_along(all.probs)) {
        all.probs[[i]] <- kmer_mat_to_probs_k1_cpp(all.paths[[i]], bkg.probs,
                                                   alph.sort.split[[i]])
      }
    # },
#
    # "1" = {
      # for (i in seq_along(all.probs)) {
        # all.probs[[i]] <- kmer_mat_to_probs_k2_cpp(all.paths[[i]], bkg.probs[lets1],
                                                   # alph.sort.split[[i]],
                                                   # bkg.probs[lets2])
      # }
    # },
#
    # "2" = {
      # for (i in seq_along(all.probs)) {
        # all.probs[[i]] <- kmer_mat_to_probs_k3_cpp(all.paths[[i]], bkg.probs[lets1],
                                                   # alph.sort.split[[i]],
                                                   # bkg.probs[lets2], bkg.probs[lets3])
      # }
    # }
#
  # )

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

branch_and_bound_kmers <- function(score.mat, min.score) {

  max.scores <- c(rev(cumsum(rev(apply(score.mat, 2, max)))), 0L)

  mot_len <- ncol(score.mat)

  paths <- init_paths_cpp(score.mat, min.score, max.scores[2])
  if (mot_len == 1) return(paths)

  for (i in seq_len(mot_len - 1) + 1) {
    paths <- calc_next_path_cpp(score.mat, paths, min.score, max.scores[i + 1])
  }

  paths

}

motif_score <- function(score.mat, pval, k = 8) {

  # Assumes a _uniform_ background! (or else distribution no longer normal)

  pval <- 1 - pval
  alph.len <- nrow(score.mat)
  mot.len <- ncol(score.mat)

  score.mat <- matrix(score.mat, nrow = alph.len)

  max.score <- sum(apply(score.mat, 2, max))
  min.score <- sum(apply(score.mat, 2, min))

  if (mot.len > k) {

    # Not sure it works correctly if the splits are not the same size?

    score.mat.split <- split_mat(score.mat, k)
    s.split <- lapply(score.mat.split,
                      function(x) rowSums(expand.grid(as.data.frame(x))))

    mean.split <- vapply(s.split, mean, numeric(1))
    var.split <- vapply(s.split, var, numeric(1))

    answer <- qnorm(pval, sum(mean.split), sqrt(sum(var.split)))

    if (answer < min.score) answer <- min.score
    else if (answer > max.score) answer <- max.score

  } else {

    p <- expand.grid(as.data.frame(score.mat))
    s <- rowSums(p)
    e <- ecdf(s)

    answer <- unname(quantile(e, probs = pval))

  }

  answer

}

split_mat <- function(score.mat, k) {

  mot.len <- ncol(score.mat)

  times.tosplit <- mot.len %/% k
  leftover.split <- mot.len %% k
  mot.split <- vector("list", times.tosplit + ifelse(leftover.split > 0, 1, 0))
  mot.split[[1]] <- score.mat[, seq_len(k)]

  if (times.tosplit > 1) {

    for (i in seq_len(times.tosplit - 1)) {
      mot.split[[i + 1]] <- score.mat[, (i * k + 1):(i * k + k)]
    }

  }

  if (leftover.split > 0) {

    mot.split[[length(mot.split)]] <- score.mat[, (mot.len - leftover.split + 1):
                                                   mot.len]

    if (!is.matrix(mot.split[[length(mot.split)]])) {
      mot.split[[length(mot.split)]] <- matrix(mot.split[[length(mot.split)]])
    }

  }

  mot.split

}
