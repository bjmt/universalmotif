#' Motif P-value and scoring utility
#'
#' For calculating p-values/logodds scores for any number of motifs.
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param score `numeric`, `list` Get a p-value for a motif from a logodds score.
#' @param pvalue `numeric`, `list` Get a logodds score for a motif from a
#'    p-value.
#' @param bkg.probs `numeric`, `list` A vector background probabilities.
#'    If supplying individual background
#'    probabilities for each motif, a list of such vectors. If missing,
#'    retrieves the background from the motif `bkg` slot. Note that this
#'    option is only used when providing a P-value and returning a score; for
#'    the inverse, the motifs are first converted to PWMs via [convert_type()],
#'    which uses the motif `bkg` slot for background adjustment.
#' @param use.freq `numeric(1)` By default uses the regular motif matrix;
#'    otherwise uses the corresponding `multifreq` matrix. Max is 3.
#' @param k `numeric(1)` For speed, scores/p-values can be approximated after
#'    subsetting the motif every `k` columns. If `k` is a value
#'    equal or higher to the size of input motif(s), then the calculations
#'    are (nearly) exact. The default, 8, is recommended to those looking for
#'    a good tradeoff between speed and accuracy for jobs requiring repeated
#'    calculations.
#' @param nthreads `numeric(1)` Run [motif_pvalue()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#' @param rand.tries `numeric(1)` When `ncol(motif) < k`, an approximation is
#'    used. This involves randomly approximating the overall
#'    motif score distribution. To increase accuracy, the distribution is
#'    approximated `rand.tries` times and the final scores averaged.
#' @param rng.seed `numeric(1)` In order to allow [motif_pvalue()] to perform
#'    C++ level parallelisation, it must work independently from R. This means
#'    it cannot communicate with R to get/set the R RNG state. To get around
#'    this, the RNG seed used by the C++ function can be set with `rng.seed`.
#'    To make sure each thread gets a different seed however, the seed
#'    is multiplied with the iteration count. For example: when working with
#'    two motifs, the second motif gets the following seed: `rng.seed * 2`.
#'     The default is to pick a random
#'    number as chosen by [sample()], which effectively makes [motif_pvalue()]
#'    dependent on the R RNG state.
#' @param allow.nonfinite `logical(1)` If `FALSE`, then apply a pseudocount if
#'    non-finite values are found in the PWM. Note that if the motif has a
#'    pseudocount greater than zero and the motif is not currently of type PWM,
#'    then this parameter has no effect as the pseudocount will be
#'    applied automatically when the motif is converted to a PWM internally. This
#'    value is set to `FALSE` by default in order to stay consistent with
#'    pre-version 1.8.0 behaviour.
#'
#' @return `numeric` A vector of scores/p-values.
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
#' motif length increases even this calculation becomes impractical. To get
#' around this, the p-value calculation can be approximated.
#'
#' In order to calculate p-values for longer motifs, this function uses the
#' approximation proposed by Hartmann et al. (2013), where
#' the motif is subset, p-values calculated for the subsets, and finally
#' combined for a total p-value. The smaller the size of the subsets, the
#' faster the calculation; but also, the bigger the approximation. This can be
#' controlled by setting `k`. In fact, for smaller motifs (< 13 positions)
#' calculating exact p-values can be done individually in reasonable time by
#' setting `k = 12`.
#'
#' To calculate a score from a P-value, all possible scores are calculated
#' and the `(1 - pvalue) * 100` nth percentile score returned.
#' When `k < ncol(motif)`, the complete set of scores is instead approximated
#' by randomly adding up all possible scores from each subset.
#' Note that this approximation
#' can actually be potentially quite expensive at times and even slower than
#' the exact version; for jobs requiring lots of repeat calculations, a bit of
#' benchmarking beforehand can be useful to find the optimal settings.
#'
#' To get an idea as to how the score calculator works (without approximation),
#' try the following code with your motif (be careful with longer motifs):
#'
#' `quantile(get_scores(motif), probs = 0.99)`
#'
#' @references
#'
#' Hartmann H, Guthohrlein EW, Siebert M, Soding SLJ (2013).
#' “P-value-based regulatory motif discovery using positional weight
#' matrices.” *Genome Research*, **23**, 181-194.
#'
#' @examples
#' if (R.Version()$arch != "i386") {
#'
#' ## P-value/score calculations are performed using the PWM version of the
#' ## motif
#' data(examplemotif)
#'
#' ## Get a minimum score based on a p-value
#' motif_pvalue(examplemotif, pvalue = 0.001)
#'
#' ## Get the probability of a particular sequence hit
#' motif_pvalue(examplemotif, score = 0)
#'
#' ## The calculations can be performed for multiple motifs
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
#' }
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @seealso [motif_score()]
#' @export
motif_pvalue <- function(motifs, score, pvalue, bkg.probs, use.freq = 1,
  k = 8, nthreads = 1, rand.tries = 10, rng.seed = sample.int(1e4, 1),
  allow.nonfinite = FALSE, method = c("dynamic", "exhaustive")) {

  # TODO: Need to work on support for use.freq > 1.

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
  #                       verbose = 0, threshold = 0,
  #                       threshold.type = "logodds")[1:100, ]
  # res$pvalue <- motif_pvalue(ArabidopsisMotif, score = res$score)
  #

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list( use.freq = args$use.freq, k = args$k,
                                     nthreads = args$nthreads),
                                c(1, 1, 1), c(FALSE, FALSE, FALSE),
                                TYPE_NUM)
  bkg_check <- character()
  if (!missing(bkg.probs)) {
    if (!is.list(bkg.probs) && !is.numeric(bkg.probs)) {
      bkg_check <- paste0(" * Incorrect type for 'bkg.probs': ",
                          "expected 'list' or 'numeric'; got `",
                          class(bkg.probs)[1], "`")
    }
  }
  use.freq_check <- character()
  if (use.freq > 3) {
    use.freq_check <- " * Incorrect 'use.freq': maximum allowed is 3"
  }
  all_checks <- c(num_check, bkg_check, use.freq_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  if (!missing(score)) {
    if (!is.list(score) && !is.numeric(score))
      stop(wmsg("`score` must be either a numeric vector or a list of ",
          "numeric vectors"), call. = FALSE)
    else if (is.list(score) && !all(vapply(score, is.numeric, logical(1))))
      stop(wmsg("If `score` is a list then it must be a list of ",
          "numeric vectors"), call. = FALSE)
  }
  if (!missing(pvalue)) {
    if (!is.list(pvalue) && !is.numeric(pvalue))
      stop(wmsg("`pvalue` must be either a numeric vector or a list of ",
          "numeric vectors"), call. = FALSE)
    else if (is.list(pvalue) && !all(vapply(pvalue, is.numeric, logical(1))))
      stop(wmsg("If `pvalue` is a list then it must be a list of ",
          "numeric vectors"), call. = FALSE)
  }
  #---------------------------------------------------------

  method <- match.arg(method)
  if (method == "dynamic" && allow.nonfinite)
    stop(wmsg("`method = \"dynamic\"` and `allow.nonfinite = TRUE` are not compatible"),
      call. = FALSE)

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PWM")
  if (!is.list(motifs)) motifs <- list(motifs)
  anyinf <- vapply(motifs, function(x) any(is.infinite(x@motif)), logical(1))
  if (any(anyinf) && !allow.nonfinite) {
    message(wmsg("Note: found -Inf values in motif PWM(s), adding a pseudocount. ",
      "Set `allow.nonfinite = TRUE` to prevent this behaviour when ",
      "`method = \"exhaustive\"`."))
    for (i in which(anyinf)) {
      motifs[[i]] <- suppressMessages(normalize(motifs[[i]]))
    }
  }

  motifs2 <- motifs

  if (use.freq == 1) {
    motifs <- lapply(motifs, function(x) x@motif)
  } else {
    motifs <- lapply(motifs, function(x)
      MATRIX_ppm_to_pwm(x@multifreq[[as.character(use.freq)]],
                        nsites = x@nsites,
                        pseudocount = x@pseudocount,
                        bkg = x@bkg[rownames(x@multifreq[[as.character(use.freq)]])]))
  }

  motnrows <- vapply(motifs, nrow, integer(1))
  motnrows.k <- motnrows^k
  if (method == "exhaustive" && any(motnrows.k > 1e8)) {
    while (any(motnrows.k > 1e8)) {
      k <- k - 1
      motnrows.k <- motnrows^k
    }

    # 1e8 bytes = 200 megabytes
    # 1e9 bytes = 2 gigabytes

    warning(wmsg("Be careful when using motif_pvalue(method=\"exhaustive\") for ",
        "motifs with large ",
        "alphabets or with use.freq > 1 in combination with high k ",
        "values. Currently this function does not allow use cases when ",
        "nrow(motif)^k > 1e8 (or the respective use.freq slot). ",
        "Continuing with k=", k, "."), immediate. = TRUE)
  }

  if (!missing(bkg.probs)) {

    if (!is.list(bkg.probs)) bkg.probs <- list(bkg.probs)
    if (length(bkg.probs) != length(motifs) && length(bkg.probs) > 1) {
      stop(wmsg("`bkg.probs` must be either a single numeric vector set of ",
          "background probabilities ",
          "or a list of background probabilities equal to the number of motifs"),
        call. = FALSE)
    }
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

  if (!missing(score) && missing(pvalue)) {

    wasList <- FALSE
    if (is.list(score)) wasList <- TRUE
    input <- sanitize_input(motifs, bkg.probs, score, method)

    if (method == "exhaustive") {
      out <- motif_pvalue_cpp(input$motifs, input$bkg.probs, input$x,
        k, nthreads, allow.nonfinite)
      if (wasList) out <- restore_list(out, input$nM, input$nX)
    } else if (method == "dynamic") {
      out <- motif_pvalue_dynamic(input$motifs, input$bkg.probs, input$x)
      if (!wasList) out <- out[[1]]
    }

  } else if (missing(score) && !missing(pvalue)) {

    wasList <- FALSE
    if (is.list(pvalue)) wasList <- TRUE
    input <- sanitize_input(motifs, bkg.probs, pvalue, method)

    if (method == "exhaustive") {
      out <- motif_score_cpp(input$motifs, input$x, rng.seed, k, nthreads,
        rand.tries, allow.nonfinite)
      if (allow.nonfinite) {
        out[out <= -min_max_doubles()$max] <- -Inf
      }
      if (wasList) out <- restore_list(out, input$nM, input$nX)
    } else if (method == "dynamic") {
      out <- motif_score_dynamic(input$motifs, input$bkg.probs, input$x)
      if (!wasList) out <- out[[1]]
    }

  } else if (missing(score) && missing(pvalue)) {

    stop("Both 'score' and 'pvalue' cannot be missing")

  } else {

    stop("only one of 'score' and 'pvalue' can be used at a time")

  }

  out

}

restore_list <- function(x, nM, nX) {
  if (length(nX) == 1) {
    i <- rep(seq_len(nM), each = nX)
    unname(tapply(x, i, list))
  } else {
    out <- vector("list", nM)
    nX <- cumsum(nX)
    out[[1]] <- x[1:nX[1]]
    for (i in seq_along(nX)[-1]) {
      out[[i]] <- x[(nX[i - 1] + 1):nX[i]]
    }
    out
  }
}

sanitize_input <- function(mots, bkgs, x, method) {

  # Note: having this function is somewhat insane, but I am trying to stay
  # backwards-compatible.

  if (is.numeric(x)) {
    if (any(is.infinite(x))) {
      stop(wmsg("`score`/`pvalue` cannot be non-finite"), call. = FALSE)
    }
  } else if (is.list(x)) {
    if (any(vapply(x, function(y) any(is.infinite(y)), logical(1)))) {
      stop(wmsg("`score`/`pvalue` cannot be non-finite"), call. = FALSE)
    }
  }

  nM <- NULL
  nX <- NULL

  if (method == "dynamic") {

    if (length(mots) == 1 && is.numeric(x)) {
      x <- list(x)
    } else if (length(mots) == 1 && is.list(x)) {
      if (length(x) != length(mots)) {
        stop(wmsg("If `pvalue`/`score` is a `list`, then it must be a single ",
            "list item or be of equal length to the number of motifs"),
          call. = FALSE)
      }
    } else if (length(mots) > 1 && is.numeric(x)) {
      if (length(x) == 1) {
        x <- as.list(rep_len(x, length(mots)))
      } else if (length(x) == length(mots)) {
        x <- as.list(x)
      } else {
        stop(wmsg("If `pvalue/`score` is a `numeric` vector, then it must be ",
            "of length one or of equal length to the number of motifs"),
          call. = FALSE)
      }
    } else if (length(mots) > 1 && is.list(x)) {
      if (length(x) == 1) {
        x <- rep_len(x, length(mots))
      } else if (length(x) != length(mots)) {
        stop(wmsg("If `pvalue`/`score` is a `list`, then it must be a single ",
            "list item or be of equal length to the number of motifs"),
          call. = FALSE)
      }
    }

  } else if (method == "exhaustive") {

    if (length(mots) == 1 && is.numeric(x)) {
      if (length(x) > 1) {
        mots <- rep_len(mots, length(x))
        bkgs <- rep_len(bkgs, length(x))
      }
    } else if (length(mots) == 1 && is.list(x)) {
      if (length(x) > 1) {
        stop(wmsg("If `pvalue`/`score` is a `list`, then it must be a single ",
            "list item or be of equal length to the number of motifs"),
          call. = FALSE)
      }
      x <- x[[1]]
      nM <- 1
      nX <- length(x)
      mots <- rep_len(mots, length(x))
      bkgs <- rep_len(bkgs, length(x))
    } else if (length(mots) > 1 && is.numeric(x)) {
      if (length(x) == 1) {
        x <- rep_len(x, length(mots))
      } else if (length(x) != length(mots)) {
        stop(wmsg("If `pvalue/`score` is a `numeric` vector, then it must be ",
            "of length one or of equal length to the number of motifs"),
          call. = FALSE)
      }
    } else if (length(mots) > 1 && is.list(x)) {
      if (length(x) == 1) {
        nM <- length(mots)
        nX <- length(x[[1]])
        mots <- rep(mots, each = nX)
        bkgs <- rep(bkgs, each = nX)
        x <- unlist(rep(x, each = nM))
      } else if (length(x) != length(mots)) {
        stop(wmsg("If `pvalue`/`score` is a `list`, then it must be a single ",
            "list item or be of equal length to the number of motifs"),
          call. = FALSE)
      } else {
        nM <- length(mots)
        nX <- vapply(x, length, integer(1))
        mots <- unlist(mapply(function(z1, z2) rep(list(z1), length(z2)),
          mots, x, SIMPLIFY = FALSE), recursive = FALSE)
        bkgs <- unlist(mapply(function(z1, z2) rep(list(z1), length(z2)),
          bkgs, x, SIMPLIFY = FALSE), recursive = FALSE)
        x <- unlist(x)
      }
    }

  }

  list(motifs = mots, bkg.probs = bkgs, x = x, nM = nM, nX = nX)

}

# TODO: multithreaded c++
motif_score_dynamic <- function(motifs, bkg.probs, pvalue) {
  mapply(motif_score_dynamic_single, motifs, bkg.probs, pvalue, SIMPLIFY = FALSE)
}

motif_pvalue_dynamic <- function(motifs, bkg.probs, score) {
  mapply(motif_pvalue_dynamic_single, motifs, bkg.probs, score, SIMPLIFY = FALSE)
}

motif_score_dynamic_single <- function(mot, bkg, p) {
  maxs <- sum(apply(mot, 2, max))
  mins <- abs(sum(as.integer(apply(mot, 2, min) * 1000)))
  mcdf <- motif_cdf(mot, bkg)
  pindex <- vapply(p, function(x) which.max(mcdf < x), integer(1))
  res <- (pindex - mins) / 1000
  res[res > maxs] <- maxs
  res
}

motif_pvalue_dynamic_single <- function(mot, bkg, s) {
  mins <- abs(sum(as.integer(apply(mot, 2, min) * 1000)))
  sint <- as.integer(s * 1000) + mins
  mcdf <- motif_cdf(mot, bkg)
  mcdf[sint]
}

# uh oh, doesnt work with ArabidopsisMotif
motif_cdf <- function(mot, bkg) {
  mot <- mot * 1000
  mot <- apply(mot, 1:2, as.integer)
  motmin <- abs(min(mot))
  mot <- mot + motmin
  motwidth <- ncol(mot)
  alphlen <- nrow(mot)
  maxscore <- max(mot)
  # maxscore <- sum(apply(mot, 2, max))
  cdflen <- motwidth * maxscore + 1
  pdfnew <- rep(1, cdflen)
  for (i in seq_len(motwidth)) {
    maxstep <- (i - 1) * maxscore + 1
    pdfold <- pdfnew
    pdfnew[seq_len(maxstep + maxscore)] <- 0
    for (j in seq_len(alphlen)) {
      s <- mot[j, i]
      for (k in seq_len(maxstep)) {
        if (pdfold[k] != 0) {
          pdfnew[k + s] <- pdfnew[k + s] + pdfold[k] * bkg[j]
        }
      }
    }
  }
  rev(cumsum(rev(pdfnew / sum(pdfnew))))
}

scale_val <- function(x, w, scale, offset) {
  (x - w * offset) * scale
}

scale_mot <- function(mot) {
  offset <- min(mot)
  scale <- (max(mot) * ncol(mot)) / (max(mot) - min(mot))
  apply(mot, 1:2, function(x) scale_val(x, 1, scale, offset))
}

scale_scores <- function(mot, scores) {
  offset <- min(mot)
  scale <- (max(mot) * ncol(mot)) / (max(mot) - min(mot))
  vapply(scores, function(x) scale_val(x, 1, scale, offset), numeric(1))
}

motif_pvalue_dynamic_single3 <- function(mot, bkg, s) {
  s <- scale_scores(mot, s)
  print(s)
  plen <- max(mot) * ncol(mot)
  print(plen)
  smot <- scale_mot(mot)
  print(smot)
  # Ismot <- apply(smot, 1:2, function(x) as.integer(x * 1000))
  Ismot <- apply(smot, 1:2, function(x) as.integer(x))
  print(Ismot)
  # Is <- as.integer(s * 1000)
  Is <- as.integer(s) 
  print(Is)
  # mcdf <- motif_cdf2(Ismot, as.integer(plen * 1000), bkg)
  mcdf <- motif_cdf2(Ismot, max(Ismot), bkg)
  mcdf[Is]
}

motif_pvalue_dynamic_single2 <- function(mot, bkg, s) {
  mot <- matrix(as.integer(mot * 1000), ncol = ncol(mot))

  motmin <- abs(sum(apply(mot, 2, min)))
  motmax <- abs(sum(apply(mot, 2, max)))
  motrange <- length(seq(from = -motmin, to = motmax - 1, by = 1))
  print(motrange)
  print(mot)

  mot <- mot + abs(min(mot))
  # mot <- mot + motmin
  motmax <- sum(apply(mot, 2, max))
  motrange <- length(seq(from = 1, to = motmax, by = 1))
  print(motrange)
  print(mot)

  sint <- as.integer(s * 1000) + motmin
  mcdf <- motif_cdf2(mot, motmax, bkg)
  mcdf[sint]
}
motif_cdf2 <- function(mot, motmax, bkg) {
  # pdfnew <- get_pdf(mot, max(mot), bkg)
  pdfnew <- get_pdf(mot, motmax, bkg)
  rev(cumsum(rev(pdfnew / sum(pdfnew))))
}

# motif_cdf <- function(mot, bkg) {
#   mot <- mot * 1000
#   mot <- apply(mot, 1:2, as.integer)
#   motmin <- abs(min(mot))
#   mot <- mot + motmin
#   motwidth <- ncol(mot)
#   alphlen <- nrow(mot)
#   maxscore <- max(mot)
#   cdflen <- motwidth * maxscore + 1
#   pdfnew <- rep(1, cdflen)
#   for (i in seq_len(motwidth)) {
#     maxstep <- (i - 1) * maxscore + 1
#     pdfold <- pdfnew
#     pdfnew[seq_len(maxstep + maxscore)] <- 0
#     for (j in seq_len(alphlen)) {
#       s <- mot[j, i]
#       for (k in seq_len(maxstep)) {
#         if (pdfold[k] != 0) {
#           pdfnew[k + s] <- pdfnew[k + s] + pdfold[k] * bkg[j]
#         }
#       }
#     }
#   }
#   rev(cumsum(rev(pdfnew / sum(pdfnew))))
# }

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
