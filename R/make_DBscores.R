#' Create P-value databases.
#'
#' Generate data used by [compare_motifs()] for P-value calculations. By default,
#' [compare_motifs()] uses an internal database based on the JASPAR2018 core motifs
#' \insertCite{jaspar}{universalmotif}.
#'
#' @param db.motifs `list` Database motifs.
#' @param method `character` Any of `c('PCC', 'MPCC', 'EUCL', 'MEUCL',
#'    'SW', 'MSW', 'KL', 'MKL')`. See [compare_motifs()]. If multiple methods
#'    are provided, each result will be in its own list.
#' @param shuffle.db `logical(1)` Shuffle `db.motifs` rather than
#'    generate random motifs with [create_motif()].
#' @param shuffle.k `numeric(1)` See [shuffle_motifs()].
#' @param shuffle.method `character(1)` See [shuffle_motifs()].
#' @param rand.tries `numeric(1)` Number of random motifs to create for
#'    P-value computation.
#' @param widths `numeric` Motif widths to use in P-value database calculation.
#' @param min.position.ic `numeric(1)` Minimum information content required between
#'    individual alignment positions for it to be counted in the final alignment
#'    score. It is recommended to use this together with `normalise.scores = TRUE`,
#'    as this will help punish scores resulting from only a fraction of an
#'    alignment.
#' @param normalise.scores `logical(1)` See [compare_motifs()].
#' @param min.overlap `numeric(1)` Minimum required motif overlap. See
#'    [compare_motifs()].
#' @param min.mean.ic `numeric(1)` See [compare_motifs()].
#' @param progress `logical(1)` Show progress.
#' @param BP `logical(1)` Deprecated. See `nthreads`.
#' @param nthreads `numeric(1)` Run [compare_motifs()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#'
#' @return A `data.frame` with score distributions for the
#'    input database, or a `list` with a `data.frame` for each method and
#'    an additional `list` entry logging function parameters.
#'
#' @examples
#' \dontrun{
#' library(MotifDb)
#' motifs <- convert_motifs(MotifDb[1:100])
#' make_DBscores(motifs, method = "PCC")
#' }
#'
#' @references
#'    \insertRef{jaspar}{universalmotif}
#'
#' @seealso [compare_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
make_DBscores <- function(db.motifs, method, shuffle.db = TRUE,
                          shuffle.k = 3, shuffle.method = "linear",
                          rand.tries = 100, widths = 5:30,
                          min.position.ic = 0,
                          normalise.scores = TRUE, min.overlap = 6,
                          min.mean.ic = 0.25, progress = TRUE, BP = FALSE,
                          nthreads = 1) {

  args <- as.list(environment())

  if (length(method) > 1) {
    out <- vector("list", length(method) + 1)
    names(out) <- c(method, "args")
    for (m in method) {
      out[[m]] <- make_DBscores(db.motifs, m, shuffle.db, shuffle.k, shuffle.method,
                                rand.tries, widths, min.position.ic,
                                normalise.scores, min.overlap, min.mean.ic,
                                progress, BP, nthreads)
    }
    out$args <- args[-1]
    return(out)
  }

  # param check --------------------------------------------
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
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  db.motifs <- convert_motifs(db.motifs)
  db.motifs <- convert_type(db.motifs, "PPM")

  db.ncols <- vapply(db.motifs, function(x) ncol(x@motif), numeric(1))

  rand.mots <- vector("list", length(widths))
  names(rand.mots) <- as.character(widths)
  rand.ncols <- rep(FALSE, length(widths))
  names(rand.ncols) <- as.character(widths)
  if (shuffle.db) {

    for (i in widths) {
      if (!any(db.ncols == i)) next;
      tmp <- list()
      while (length(tmp) < rand.tries) {
        tmp <- c(tmp, shuffle_motifs(db.motifs[db.ncols == i], k = shuffle.k,
                                     method = shuffle.method))
      }
      tmp <- tmp[seq_len(rand.tries)]
      rand.mots[[as.character(i)]] <- tmp
      rand.ncols[as.character(i)] <- TRUE
    }

  } else {
    for (i in widths) {
      rand.mots[[as.character(i)]] <- lapply(rand.tries, function(x) create_motif(i))
    }
  }

  totry <- expand.grid(list(target = as.numeric(widths[rand.ncols]),
                            subject = sort(unique(db.ncols))))[, 2:1]
  totry$mean <- rep(NA_real_, nrow(totry))
  totry$sd <- rep(NA_real_, nrow(totry))

  res <- vector("list", length(unique(db.ncols)))
  allcomps <- vector("list", length(unique(db.ncols)))
  names(res) <- as.character(sort(unique(db.ncols)))
  names(allcomps) <- as.character(sort(unique(db.ncols)))

  rand.mots <- unlist(rand.mots, recursive = FALSE)
  randcols <- vapply(rand.mots, function(x) ncol(x@motif), integer(1))

  if (progress) print_pb(0)
  counter <- 1
  total <- length(unique(db.ncols))

  for (i in sort(unique(db.ncols))) {

    tmp <- db.motifs[db.ncols == i]

    tmpall <- c(tmp, rand.mots)
    mtmp <- lapply(tmpall, function(x) x@motif)
    btmp <- lapply(tmpall, function(x) x@bkg)
    comps <- get_comp_indices(seq_along(tmp), length(tmpall))
    comps <- comps[!comps[, 2] %in% seq_along(tmp), ]
    allcomps[[as.character(i)]] <- rep(randcols, length(tmp))
    res[[as.character(i)]] <- compare_motifs_cpp(mtmp, comps[, 1] - 1, comps[, 2] - 1,
                                                 method, min.overlap, FALSE,
                                                 btmp, 1, FALSE, min.mean.ic,
                                                 normalise.scores, nthreads,
                                                 min.position.ic)

    if (progress) update_pb(counter, total)
    counter <- counter + 1

  }

  for (i in seq_len(nrow(totry))) {
    tmp <- res[[as.character(totry[i, 1])]]
    totry$mean[i] <- mean(tmp[allcomps[[as.character(totry[i, 1])]] == totry[i, 2]])
    totry$sd[i] <- sd(tmp[allcomps[[as.character(totry[i, 1])]] == totry[i, 2]])
  }

  totry$method <- rep(method, nrow(totry))
  totry$normalised <- rep(normalise.scores, nrow(totry))
  totry

}

