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
#'    are (nearly) exact. The default, 8, is recommended to those looking for
#'    a good tradeoff between speed and accuracy for jobs requiring repeated
#'    calculations.
#' @param progress `logical(1)` Deprecated. Does nothing.
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [motif_pvalue()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for exceptionally large
#'    jobs. Note that this is only used for calculating scores from P-values
#'    (in other words, when the `pvalue` argument is provided).
#' @param nthreads `numeric(1)` Run [motif_pvalue()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads. Note that this is only
#'    used for calculating P-values from scores (in other words, when the `score`
#'    argument is provided).
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
#' each motif subsets are combined to estimate the normal distribution of all
#' possible scores using \code{\link[stats:Normal]{stats::qnorm()}}:
#'
#' `qnorm(pvalue, mean = sum(subset.means), sd = sqrt(sum(subset.vars)))`
#'
#' For calculating exact scores, [stats::ecdf()] and [stats::quantile()] are
#' used:
#'
#' `quantile(ecdf(scores), probs = pvalue)`
#'
#' It is important to keep in mind that the approximate
#' calculation assumes a normal distribution of scores, which is rarely
#' accurate.
#'
#' @examples
#' data(examplemotif)
#'
#' ## p-value/score calculations are performed using the PWM version of the
#' ## motif; these calculations do not work if any -Inf values are present
#' examplemotif["pseudocount"] <- 1
#' # or
#' examplemotif <- normalize(examplemotif)
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
#' ## Compare score thresholds and P-value:
#' scores <- motif_score(examplemotif, c(0.6, 0.7, 0.8, 0.9))
#' motif_pvalue(examplemotif, scores)
#'
#' ## Calculate the probability of getting a certain match or better:
#' TATATAT <- score_match(examplemotif, "TATATAT")
#' TATATAG <- score_match(examplemotif, "TATATAG")
#' motif_pvalue(examplemotif, TATATAT)
#' motif_pvalue(examplemotif, TATATAG)
#'
#' ## Get all possible matches by P-value:
#' get_matches(examplemotif, motif_pvalue(examplemotif, pvalue = 0.0001))
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [motif_score()]
#' @export
motif_pvalue <- function(motifs, score, pvalue, bkg.probs, use.freq = 1,
                         k = 8, progress = FALSE, BP = FALSE, nthreads = 1) {

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
                                     use.freq = args$use.freq, k = args$k,
                                     nthreads = args$nthreads),
                                c(0, 0, 1, 1, 1), c(TRUE, TRUE, FALSE, FALSE, FALSE),
                                TYPE_NUM)
  logi_check <- check_fun_params(list(progress = args$progress, BP = args$BP),
                                 numeric(), logical(), TYPE_LOGI)
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

  if (progress)
    warning("'progress' is deprecated and does nothing")
  progres <- FALSE

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PWM")
  if (!is.list(motifs)) motifs <- list(motifs)
  anyinf <- vapply(motifs, function(x) any(is.infinite(x@motif)), logical(1))
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
    motifs <- lapply(motifs, function(x) x@motif)
  } else {
    mots <- lapply(motifs, function(x) x@multifreq[[as.character(use.freq)]])
    motifs <- mapply(function(x, y) apply(x, 2, ppm_to_pwmC, bkg = numeric(),
                                          pseudocount = y@pseudocount,
                                          nsites = y@nsites),
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
        stop("length(bkg.probs) must match nrow(motif)^use.freq")

    } else bkg.probs <- rep(list(NULL), length(motifs))

    bkg.probs <- mapply(motif_pvalue_bkg, motifs2, bkg.probs,
                        MoreArgs = list(use.freq = use.freq),
                        SIMPLIFY = FALSE)
    
    l1 <- length(motifs)
    l2 <- length(bkg.probs)
    l3 <- length(score)
    lall <- max(c(l1, l2, l3))
    motifs <- rep_len(motifs, lall)
    bkg.probs <- rep_len(bkg.probs, lall)
    score <- rep_len(score, lall)

    out <- motif_pvalue_cpp(motifs, bkg.probs, score, k, nthreads)

  } else if (missing(score) && !missing(pvalue)) {

    out <- mapply_(motif_score_pval, motifs, pvalue, k, PB = progress, BP = BP)

  } else if (missing(score) && missing(pvalue)) {

    stop("Both 'score' and 'pvalue' cannot be missing")

  } else {

    stop("only one of 'score' and 'pvalue' can be used at a time")

  }

  out

}

motif_pvalue_bkg <- function(motif, bkg.probs, use.freq) {

  lets1 <- rownames(motif@motif)
  if (use.freq > 1) lets2 <- get_klets(lets1, use.freq)
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

motif_score_pval <- function(score.mat, pval, k = 8) {

  # Assumes a _uniform_ background! (or else distribution no longer normal)

  # Previously used rowSums(expand.grid(as.data.frame(scores)) instead of
  # expand_scores(scores). The latter function ends up with 1/3 the memory
  # allocations compared to the expand.grid solution.

  pval <- 1 - pval
  alph.len <- nrow(score.mat)
  mot.len <- ncol(score.mat)

  score.mat <- matrix(as.integer(score.mat * 1000), nrow = alph.len)

  max.score <- sum(apply(score.mat, 2, max)) / 1000
  min.score <- sum(apply(score.mat, 2, min)) / 1000

  if (mot.len > k) {

    score.mat.split <- split_mat(score.mat, k)
    s.split <- lapply(score.mat.split,
                      function(x) expand_scores(x) / 1000)

    mean.split <- vapply(s.split, mean, numeric(1))
    var.split <- vapply(s.split, var, numeric(1))

    answer <- qnorm(pval, sum(mean.split), sqrt(sum(var.split)))

    if (answer < min.score) answer <- min.score
    else if (answer > max.score) answer <- max.score

  } else {

    s <- expand_scores(score.mat) / 1000
    e <- ecdf(s)

    answer <- unname(quantile(e, probs = pval))

  }

  answer

}

split_mat <- function(score.mat, k) {

  mot.len <- ncol(score.mat)

  times.tosplit <- mot.len %/% k
  leftover.split <- mot.len %% k
  splitl <- times.tosplit + ifelse(leftover.split > 0, 1, 0)
  mot.split <- vector("list", splitl)
  mot.split[[1]] <- score.mat[, seq_len(k)]

  if (times.tosplit > 1) {

    for (i in seq_len(times.tosplit - 1)) {
      mot.split[[i + 1]] <- score.mat[, (i * k + 1):(i * k + k)]
    }

  }

  if (leftover.split > 0) {

    mot.split[[splitl]] <- score.mat[, (mot.len - leftover.split + 1):mot.len]

    if (!is.matrix(mot.split[[splitl]])) {
      mot.split[[splitl]] <- matrix(mot.split[[splitl]])
    }

  }

  mot.split

}
