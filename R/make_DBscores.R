#' Create P-value databases.
#'
#' Generate data used by [compare_motifs()] for P-value calculations. By default,
#' [compare_motifs()] uses an internal database based on the JASPAR2018 core motifs
#' \insertCite{jaspar}{universalmotif}. Parameters for distributions are
#' are estimated for every combination of motif `widths`.
#'
#' @param db.motifs `list` Database motifs.
#' @param shuffle.db `logical(1)` Deprecated. Does nothing.
#'    generate random motifs with [create_motif()].
#' @param shuffle.k `numeric(1)` See [shuffle_motifs()].
#' @param shuffle.method `character(1)` See [shuffle_motifs()].
#' @param rand.tries `numeric(1)` Approximate number of comparisons 
#'    to perform for every combination of `widths`.
#' @param widths `numeric` Motif widths to use in P-value database calculation.
#' @param progress `logical(1)` Show progress.
#' @param BP `logical(1)` Deprecated. See `nthreads`.
#' @param nthreads `numeric(1)` Run [compare_motifs()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'
#' @return A `DataFrame` with score distributions for the
#'    input database. If more than one [make_DBscores()] run occurs (i.e. args
#'    `method`, `normalise.scores` or `score.strat` are longer than 1), then
#'    a `list` is returned with a `DataFrame` of scores and a `list` of
#'    function parameters.
#'
#' @details
#' See [compare_motifs()] for more info on comparison parameters.
#'
#' To replicate the internal \pkg{universalmotif} DB scores, run
#' [make_DBscores()] with the default settings. Note that this will be
#' a slow process.
#'
#' Arguments `widths`, `method`, `normalise.scores` and `score.strat` are
#' vectorized; all combinations will be attempted.
#'
#' Randomly generated scores are used to estimate parameters for logistic
#' distributions. This will occasionally fail for `method = "IS"`, which
#' can sometimes fit too poorly to estimate parameters. In such cases,
#' the mean and standard deviation will instead be used. A warning is
#' given in [compare_motifs()] when P-values are asked for when using
#' `method = "IS"`.
#'
#' @examples
#' \dontrun{
#' library(MotifDb)
#' motifs <- convert_motifs(MotifDb[1:100])
#' scores <- make_DBscores(motifs, method = "PCC")
#' compare_motifs(motifs, 1:100, db.scores = scores)
#' }
#'
#' @references
#'    \insertRef{jaspar}{universalmotif}
#'
#' @seealso [compare_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams compare_motifs
#' @export
make_DBscores <- function(db.motifs,
                          method = c("PCC", "EUCL", "SW", "KL", "ALLR", "BHAT",
                                     "HELL", "IS", "SEUCL", "MAN", "ALLR_LL"),
                          shuffle.db = TRUE,
                          shuffle.k = 3, shuffle.method = "linear",
                          rand.tries = 1000, widths = 5:30,
                          min.position.ic = 0,
                          normalise.scores = c(FALSE, TRUE), min.overlap = 6,
                          min.mean.ic = 0.25, progress = TRUE, BP = FALSE,
                          nthreads = 1, tryRC = TRUE,
                          score.strat = c("sum", "a.mean", "g.mean", "median")) {

  args <- as.list(environment())

  if (length(method) > 1 || length(normalise.scores) > 1 || length(score.strat) > 1) {

    # Need to vectorize through all three args and combine into data.frames, one per
    # method. Finally, all data.frames are kept in a list with an additional entry
    # for function args.

    out <- vector("list", length(method) + 1)
    names(out) <- c(method, "args")
    mc <- 1
    total <- length(method) * length(normalise.scores) * length(score.strat)

    if (progress) t1 <- Sys.time()
    
    for (m in method) {

      out[[m]] <- DataFrame(subject = integer(), target = integer(),
                             paramA = numeric(), paramB = numeric(),
                             method = character(), normalised = logical(),
                             strat = character(), distribution = character())

      for (norm in normalise.scores) {
        for (strat in score.strat) {

          if (progress) message("[", mc, "/", total, "]", " method=\"", m,
                                "\" normalise.scores=", norm, " score.strat=\"",
                                strat, "\" ", appendLF = FALSE)
          mc <- mc + 1
          if (strat == "g.mean" && m %in% c("ALLR", "ALLR_LL", "PCC")) {
            if (progress)
              message("\n > Skipping: g.mean not allowed with ALLR/ALLR_LL/PCC\n")
            next
          }
          if (progress) start <- Sys.time()

          tmp <- make_DBscores(db.motifs, m, shuffle.db, shuffle.k, shuffle.method,
                               rand.tries, widths, min.position.ic,
                               norm, min.overlap, min.mean.ic,
                               progress, BP, nthreads, tryRC, strat)
          out[[m]] <- rbind(out[[m]], tmp)

          if (progress) stop <- Sys.time()
          if (progress) message(" > ", format(difftime(stop, start)), "\n")


        }
      }

    }

    if (progress) t2 <- Sys.time()
    if (progress) message(" *** Total runtime: ", format(difftime(t2, t1)), " ***")

    out <- list(scores = do.call(rbind, out))
    for (i in colnames(out$scores)) {
      out$scores[, i] <- Rle(out$scores[, i])
    }
    rownames(out$scores) <- NULL

    out$args <- args[-1]
    return(out)

  }

  # param check --------------------------------------------
  method <- match.arg(method, COMPARE_METRICS)
  char_check <- check_fun_params(list(method = args$method,
                                      shuffle.method = args$shuffle.method,
                                      score.strat = args$score.strat),
                                 c(0, 1, 1), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(shuffle.k = args$shuffle.k,
                                     rand.tries = args$rand.tries,
                                     min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     nthreads = args$nthreads,
                                     min.position.ic = args$min.position.ic),
                                numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(shuffle.db = args$shuffle.db,
                                      progress = args$progress, BP = args$BP,
                                      normalise.scores = args$normalise.scores,
                                      tryRC = args$tryRC),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  # having min.mean.ic > 0 can really mess with scores sometimes

  if (BP) warning("'BP' is deprecated; use 'nthreads' instead", immediate. = TRUE)

  rand.mots <- make_DBscores_motifs(db.motifs, widths, rand.tries, shuffle.k)

  comps <- get_comparisons(widths)
  total <- nrow(comps)
  totry <- DataFrame(subject = comps[, 1], target = comps[, 2],
                      paramA = rep(NA_real_, total),
                      paramB = rep(NA_real_, total),
                      method = rep(method, total),
                      normalised = rep(normalise.scores, total),
                      strat = rep(score.strat, total),
                      distribution = rep("logistic", total))

  if (progress) {
    print_pb(0)
  }
  counter <- 1

  for (i in seq_len(total)) {

    subject <- rand.mots[[as.character(totry$subject[i])]]
    target <- rand.mots[[as.character(totry$target[i])]]

    tmpall <- c(subject, target)

    mtmp <- lapply(tmpall, function(x) x@motif)

    comps <- get_comp_indices(seq_along(subject), length(tmpall))
    comps <- comps[!comps[, 2] %in% seq_along(subject), ]

    res <- compare_motifs_cpp(mtmp, comps[, 1] - 1, comps[, 2] - 1,
                              method, min.overlap, tryRC,
                              get_bkgs(tmpall), 1, FALSE, min.mean.ic,
                              normalise.scores, nthreads,
                              min.position.ic, get_nsites(tmpall),
                              score.strat)

    if (length(unique(res)) == 1)
      stop(wmsg("failed to estimate distribution due to uniform random scores ",
                "at comparison: ", as.vector(totry$subject[i]), " - ",
                as.vector(totry$target[i]), " ; ",
                "perhaps too few reference motifs"))

    res[res == min_max_doubles()$min | res == min_max_doubles()$max] <- NA_real_
    res <- res[!is.na(res)]
    if (length(res) <= 1)
      stop(wmsg("failed to obtain scores due to low motif IC at comparison: ",
                as.vector(totry$subject[i]), " - ", as.vector(totry$target[i]),
                " ; perhaps too few reference motifs or too many low IC motifs"))
    
    d <- find_best_dist(res)
    totry$paramA[i] <- d$paramA
    totry$paramB[i] <- d$paramB
    totry$distribution[i] <- d$dist

    if (progress) update_pb(counter, total)
    counter <- counter + 1

  }

  for (i in colnames(totry)) {
    totry[, i] <- Rle(totry[, i])
  }
  rownames(totry) <- NULL
  totry

}

find_best_dist <- function(x) {

  d <- c("normal", "logistic", "weibull")

  pvals <- numeric(length(d))

  # NORMAL
  n <- tryCatch(suppressWarnings(fitdistr(x, "normal")), error = function(e) FALSE)
  if (isFALSE(n)) {
    pvals[1] <- -1
  } else {
    pvals[1] <- suppressWarnings(ks.test(x, "pnorm", n$estimate["mean"],
                                         n$estimate["sd"])$p)
  }

  # LOGISTIC
  l <- tryCatch(suppressWarnings(fitdistr(x, "logistic")), error = function(e) FALSE)
  if (isFALSE(l)) {
    pvals[2] <- -1
  } else {
    pvals[2] <- suppressWarnings(ks.test(x, "plogis", l$estimate["location"],
                                         l$estimate["scale"])$p)
  }

  # WEIBULL
  w <- tryCatch(suppressWarnings(fitdistr(x, "weibull")), error = function(e) FALSE)
  if (isFALSE(w)) {
    pvals[3] <- -1
  } else {
    pvals[3] <- suppressWarnings(ks.test(x, "pweibull", w$estimate["shape"],
                                         w$estimate["scale"])$p)
  }

  dist <- d[which.max(pvals)[1]]

  switch(dist,
         "normal" = list(dist = "normal", paramA = n$estimate["mean"],
                         paramB = n$estimate["sd"]),
         "logistic" = list(dist = "logistic", paramA = l$estimate["location"],
                           paramB = l$estimate["scale"]),
         "weibull" = list(dist = "weibull", paramA = w$estimate["shape"],
                          paramB = w$estimate["scale"]))

}

get_comparisons <- function(widths) {
  out <- do.call(rbind, combn(widths, 2, simplify = FALSE))
  rbind(out, matrix(c(widths, widths), ncol = 2))
}

make_DBscores_motifs <- function(motifs, widths, rand.tries, shuffle.k) {

  mot.mats <- lapply(motifs, function(x) x@motif)
  mot.bkg <- get_bkgs(motifs)
  mot.nsites <- get_nsites(motifs)
  alph <- unique(vapply(motifs, function(x) x@alphabet, character(1)))
  if (length(alph) > 1) stop("all motifs must share the same alphabet")

  comb.mats <- do.call(cbind, mot.mats)
  if (ncol(comb.mats) < max(widths)) {
    while (ncol(comb.mats) < max(widths)) {
      comb.mats <- cbind(comb.mats, comb.mats)
    }
  }

  if (shuffle.k == 1) comb.mats <- comb.mats[, sample(seq_len(ncol(comb.mats)))]
  else {
    old.order <- as.character(seq_len(ncol(comb.mats)))
    new.order <- shuffle_linear(old.order, shuffle.k, mode = 2)
    comb.mats <- comb.mats[, as.numeric(new.order)]
  }

  nmots <- round(sqrt(rand.tries))
  out <- vector("list", length(widths))
  names(out) <- as.character(widths)

  for (i in seq_along(widths)) {
    j <- sample.int(length(mot.bkg), 1)
    out[[i]] <- lapply(seq_len(nmots),
                       function(x) sample_motif_(comb.mats, widths[i],
                                                 mot.bkg[[j]], mot.nsites[[j]],
                                                 alph))
  }

  out

}

sample_motif_ <- function(pool, size, bkg, nsites, alph) {

  n <- ncol(pool)
  n <- sample.int(n - size + 1, 1)
  mat <- pool[, seq(n, n + size - 1), drop = FALSE]

  create_motif(mat, alphabet = alph, bkg = bkg, nsites = nsites)

}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

make_DBscores_v1 <- function(db.motifs, method, shuffle.db = TRUE,
                             shuffle.k = 3, shuffle.method = "linear",
                             shuffle.leftovers = "asis", rand.tries = 1000,
                             normalise.scores = TRUE, min.overlap = 6,
                             min.mean.ic = 0, progress = TRUE, BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(method = args$method,
                                      shuffle.method = args$shuffle.method,
                                      shuffle.leftovers = args$shuffle.leftovers),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(shuffle.k = args$shuffle.k,
                                     rand.tries = args$rand.tries,
                                     min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic),
                                numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(shuffle.db = args$shuffle.db,
                                      progress = args$progress, BP = args$BP,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  db.motifs <- convert_motifs(db.motifs)
  db.ncols <- vapply(db.motifs, function(x) ncol(x@motif), numeric(1))

  if (shuffle.db) {
    rand.mots <- shuffle_motifs(db.motifs, k = shuffle.k,
                                method = shuffle.method) 
    if (length(rand.mots) != rand.tries) {
      if (length(rand.mots) < rand.tries) {
        while (length(rand.mots) < rand.tries) {
          more.rand.mots <- shuffle_motifs(db.motifs, k = shuffle.k,
                                           method = shuffle.method) 
          rand.mots <- c(rand.mots, more.rand.mots)
        }
      }
      if (length(rand.mots) > rand.tries) {
        rand.mots <- rand.mots[sample(seq_along(rand.mots), rand.tries)]
      }
    }
  } else {
    rand.mots <- lapply(seq_len(rand.tries),
                        function(x) create_motif(sample.int(26, 1) + 4))
  }
  rand.ncols <- vapply(rand.mots, function(x) ncol(x@motif), numeric(1))

  totry <- expand.grid(list(subject = sort(unique(rand.ncols)),
                            target = sort(unique(db.ncols))))
  totry$mean <- rep(NA, nrow(totry))
  totry$sd <- rep(NA, nrow(totry))

  res <- vector("list", nrow(totry))

  if (progress) print_pb(0)

  for (i in seq_len(nrow(totry))) {

    tmp1 <- db.motifs[totry[i, 2] == db.ncols]
    tmp2 <- rand.mots[totry[i, 1] == rand.ncols]

    res[[i]] <- compare_motifs(c(tmp2, tmp1), seq_along(tmp2), method = method,
                               min.overlap = min.overlap, min.mean.ic = min.mean.ic,
                               max.e = Inf, max.p = Inf, BP = BP, progress = FALSE,
                               normalise.scores = normalise.scores)$score

    totry$mean[i] <- mean(res[[i]])
    totry$sd[i] <- sd(res[[i]])

    if (progress) update_pb(i, nrow(totry))

  }

  totry$method <- rep(method, nrow(totry))
  totry$normalised <- rep(normalise.scores, nrow(totry))
  totry

}
