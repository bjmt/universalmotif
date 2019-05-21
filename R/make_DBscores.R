#' Create P-value databases.
#'
#' Generate data used by [compare_motifs()] for P-value calculations. By default,
#' [compare_motifs()] uses an internal database based on the JASPAR2018 core motifs
#' \insertCite{jaspar}{universalmotif}. Parameters for a logistic distribution are
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
#' @return A `data.frame` with score distributions for the
#'    input database, or a `list` with a `data.frame` for each method and
#'    an additional `list` entry logging function parameters if more than one
#'    method is provided.
#'
#' @details
#' See [compare_motifs()] for more info on comparison parameters.
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
                          method = c("PCC", "MPCC", "EUCL", "MEUCL", "SW", "MSW",
                                     "KL", "MKL", "ALLR", "MALLR", "BHAT",
                                     "MBHAT", "HELL", "MHELL", "IS", "MIS",
                                     "SEUCL", "MSEUCL", "MAN", "MMAN"),
                          shuffle.db = TRUE,
                          shuffle.k = 3, shuffle.method = "linear",
                          rand.tries = 1000, widths = 5:30,
                          min.position.ic = 0,
                          normalise.scores = FALSE, min.overlap = 1/3,
                          min.mean.ic = 0, progress = TRUE, BP = FALSE,
                          nthreads = 1, tryRC = TRUE) {

  # add a use.freq option?

  args <- as.list(environment())

  if (length(method) > 1) {
    out <- vector("list", length(method) + 1)
    names(out) <- c(method, "args")
    mc <- 1
    for (m in method) {
      if (progress) cat("Method:", paste0("[", mc, "/", length(method), "]"), m, "")
      if (progress) start <- Sys.time()
      out[[m]] <- make_DBscores(db.motifs, m, shuffle.db, shuffle.k, shuffle.method,
                                rand.tries, widths, min.position.ic,
                                normalise.scores, min.overlap, min.mean.ic,
                                progress, BP, nthreads, tryRC)
      if (progress) stop <- Sys.time()
      if (progress) cat(" >", format(difftime(stop, start)), "\n")
      mc <- mc + 1
    }
    out$args <- args[-1]
    return(out)
  }

  # param check --------------------------------------------
  method <- match.arg(method, COMPARE_METRICS)
  char_check <- check_fun_params(list(method = args$method,
                                      shuffle.method = args$shuffle.method),
                                 c(0, 1), logical(), TYPE_CHAR)
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
  totry <- data.frame(subject = comps[, 1], target = comps[, 2],
                      location = rep(NA_real_, total),
                      scale = rep(NA_real_, total),
                      method = rep(method, total),
                      normalised = rep(normalise.scores, total),
                      stringsAsFactors = FALSE)

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
                              min.position.ic, get_nsites(tmpall))

    if (length(unique(res)) == 1)
      stop(wmsg("failed to estimate logistic distribution due to uniform random scores ",
                "at comparison: ", totry$subject[i], " - ", totry$target[i], " ; ",
                "perhaps too few reference motifs"))
    a <- suppressWarnings(fitdistr(res, "logistic"))
    totry$location[i] <- a$estimate["location"]
    totry$scale[i] <- a$estimate["scale"]

    if (progress) update_pb(counter, total)
    counter <- counter + 1

  }

  totry

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
