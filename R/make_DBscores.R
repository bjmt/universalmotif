#' Create P-value databases.
#'
#' Generate data used by [compare_motifs()] for P-value calculations. By default,
#' [compare_motifs()] uses an internal database based on the JASPAR2018 core motifs
#' (Khan et al. 2018). Parameters for distributions are
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
#' @param nthreads `numeric(1)` Run [compare_motifs()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'
#' @return A `DataFrame` with score distributions for the
#'    input database. If more than one [make_DBscores()] run occurs (i.e. args
#'    `method`, `normalise.scores` or `score.strat` are longer than 1), then
#'    the function args are included in the `metadata` slot.
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
#' @examples
#' \dontrun{
#' library(MotifDb)
#' motifs <- convert_motifs(MotifDb[1:100])
#' scores <- make_DBscores(motifs, method = "PCC")
#' compare_motifs(motifs, 1:100, db.scores = scores)
#' }
#'
#' @references
#'
#' Khan A, Fornes O, Stigliani A, Gheorghe M, Castro-Mondragon JA,
#' van der Lee R, Bessy A, Cheneby J, Kulkarni SR, Tan G, Baranasic
#' D, Arenillas DJ, Sandelin A, Vandepoele K, Lenhard B, Ballester B,
#' Wasserman WW, Parcy F, Mathelier A (2018). “JASPAR 2018: update of
#' the open-access database of transcription factor binding profiles
#' and its web framework.” *Nucleic Acids Research*, **46**, D260-D266.
#'
#' @seealso [compare_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @inheritParams compare_motifs
#' @export
make_DBscores <- function(db.motifs,
                          method = c("PCC", "EUCL", "SW", "KL", "WEUCL",
                                     "ALLR", "BHAT", "HELL", "WPCC", "SEUCL",
                                     "MAN", "ALLR_LL"),
                          shuffle.db = TRUE,
                          shuffle.k = 3, shuffle.method = "linear",
                          rand.tries = 1000, widths = 5:30,
                          min.position.ic = 0,
                          normalise.scores = c(FALSE, TRUE), min.overlap = 6,
                          min.mean.ic = 0.25, progress = TRUE,
                          nthreads = 1, tryRC = TRUE,
                          score.strat = c("sum", "a.mean", "g.mean", "median",
                                          "wa.mean", "wg.mean", "fzt")) {

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
          if (strat %in% c("g.mean", "wg.mean") && m %in% c("ALLR", "ALLR_LL", "PCC")) {
            if (progress)
              message("\n > Skipping: (w)g.mean not allowed with ALLR/ALLR_LL/PCC\n")
            next
          }
          if (progress) start <- Sys.time()

          tmp <- make_DBscores(db.motifs, m, shuffle.db, shuffle.k, shuffle.method,
                               rand.tries, widths, min.position.ic,
                               norm, min.overlap, min.mean.ic,
                               progress, nthreads, tryRC, strat)
          out[[m]] <- rbind(out[[m]], tmp)

          if (progress) stop <- Sys.time()
          if (progress) message(" > ", format(difftime(stop, start)), "\n")


        }
      }

    }

    if (progress) t2 <- Sys.time()
    if (progress) message(" *** Total runtime: ", format(difftime(t2, t1)), " ***")

    out <- do.call(rbind, out)
    out <- out[order(as.vector(out$method), as.vector(out$normalised),
                     as.vector(out$strat), as.vector(out$distribution),
                     as.vector(out$subject), as.vector(out$target)), ]
    out[, 1] <- Rle(out[, 1])
    out[, 5] <- Rle(out[, 5])
    out[, 6] <- Rle(out[, 6])
    out[, 7] <- Rle(out[, 7])
    out[, 8] <- Rle(out[, 8])

    out@metadata <- args[-1]
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
                                      progress = args$progress,
                                      normalise.scores = args$normalise.scores,
                                      tryRC = args$tryRC),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  # having min.mean.ic > 0 can really mess with scores sometimes

  if (score.strat %in% c("g.mean", "wg.mean") && method %in%
      c("ALLR", "ALLR_LL", "PCC"))
    stop(wmsg("'g.mean'/'wg.mean' is not allowed for methods which can generate negative ",
              "values: ALLR, ALLR_LL, PCC"))

  if (!score.strat %in% c("sum", "a.mean", "g.mean", "median", "wa.mean",
                          "wg.mean", "fzt"))
    stop("'score.strat' must be one of 'sum', 'a.mean', 'g.mean', 'median', ",
         "'wa.mean', 'wg.mean', 'fzt'")

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

  totry <- totry[order(totry$method, totry$normalised, totry$strat,
                       totry$distribution, totry$subject, totry$target), ]
  totry[, 1] <- Rle(totry[, 1])
  totry[, 5] <- Rle(totry[, 5])
  totry[, 6] <- Rle(totry[, 6])
  totry[, 7] <- Rle(totry[, 7])
  totry[, 8] <- Rle(totry[, 8])
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
                               max.e = Inf, max.p = Inf,
                               normalise.scores = normalise.scores)$score

    totry$mean[i] <- mean(res[[i]])
    totry$sd[i] <- sd(res[[i]])

    if (progress) update_pb(i, nrow(totry))

  }

  totry$method <- rep(method, nrow(totry))
  totry$normalised <- rep(normalise.scores, nrow(totry))
  totry

}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

make_DBscores_v3 <- function(db.motifs, rand.tries = 10000, widths = 1:100,
                             progress = TRUE, BP = FALSE) {

  mots <- lapply(db.motifs, function(x) x@motif)
  cols <- do.call(cbind, mots)

  if (progress) message("Calculating scores")

  scores <- lapply_(COMPARE_METRICS, function(x) calc_scores(cols, x),
                    BP = BP, PB = progress)
  names(scores) <- COMPARE_METRICS

  if (progress) message("Calculating distribution parameters")

  params <- lapply_(scores, function(x) get_params(x, rand.tries, widths),
                    BP = BP, PB = progress)

}

get_params <- function(scores, n, widths) {
  means <- numeric(length(widths))
  sds <- numeric(length(widths))
  for (i in seq_along(widths)) {
    ans <- vapply(seq_len(n),
                  function(x) sum(sample(scores, widths[i])),
                  numeric(1))
    means[i] <- mean(ans)
    sds[i] <- sd(ans)
  }
  names(means) <- widths
  names(sds) <- widths
  list(mean = means, sd = sds)
}

calc_scores <- function(cols, method) {
  columns <- switch(method,
               KL = columns + 0.01,
               ALLR = columns + 0.01,
               ALLR_LL = columns + 0.01,
               columns
             )
  out <- numeric(ncol(columns)^2 / 2 - ncol(columns) / 2)
  counter <- 0
  for (i in seq_len(ncol(columns))) {
    for (j in seq_len(ncol(columns))[-seq_len(i)]) {
      counter <- counter + 1
      out[counter] <- compare_columns(columns[, i], columns[, j], method)
    }
  }
  out
}
