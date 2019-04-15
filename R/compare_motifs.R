#' Compare motifs.
#'
#' Compare motifs using four available metrics: Pearson correlation coefficient
#' \insertCite{pearson}{universalmotif}, Euclidean distance
#' \insertCite{euclidean}{universalmotif}, Sandelin-Wasserman similarity
#' \insertCite{wasserman}{universalmotif}, and Kullback-Leibler divergence
#' \insertCite{kldiv}{universalmotif}.
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param compare.to `numeric` If missing, compares all motifs to all other motifs.
#'    Otherwise compares all motifs to the specified motif(s).
#' @param db.scores `data.frame` See `details`.
#' @param use.freq `numeric(1)`. For comparing the `multifreq` slot.
#' @param use.type `character(1)` One of `'PPM'` and `'ICM'`.
#'    The latter allows for taking into account the background
#'    frequencies if `relative_entropy = TRUE`.
#' @param method `character(1)` One of `c('PCC', 'MPCC', 'EUCL', 'MEUCL',
#'    'SW', 'MSW', 'KL', 'MKL')`. See details.
#' @param tryRC `logical` Try the reverse complement of the motifs as well,
#'    report the best score.
#' @param min.overlap `numeric(1)` Minimum overlap required when aligning the
#'    motifs. Setting this to a number higher then the width of the motifs
#'    will not allow any overhangs. Can also be a number less than 1,
#'    representing the minimum fraction that the motifs must overlap.
#' @param min.mean.ic `numeric(1)` Minimum mean information content between the
#'    two motifs for an alignment to be scored. This helps prevent scoring
#'    alignments between low information content regions of two motifs.
#' @param relative_entropy `logical(1)` For ICM calculation. See
#'    [convert_type()].
#' @param normalise.scores `logical(1)` Favour alignments which leave fewer
#'    unaligned positions, as well as alignments between motifs of similar length.
#'    Similarity scores are multiplied by the ratio of
#'    aligned positions to the total number of positions in the larger motif,
#'    and the inverse for distance scores.
#' @param max.p `numeric(1)` Maximum P-value allowed in reporting matches.
#'    Only used if `compare.to` is set.
#' @param max.e `numeric(1)` Maximum E-value allowed in reporting matches.
#'    Only used if `compare.to` is set. The E-value is the P-value multiplied
#'    by the number of input motifs times two.
#' @param progress `logical(1)` Show progress. Not recommended if `BP = TRUE`.
#' @param BP `logical(1)` Allows the use of \pkg{BiocParallel} within
#'    [compare_motifs()]. See [BiocParallel::register()] to change the default
#'    backend. Setting `BP = TRUE` is only recommended for comparing large numbers
#'    of motifs (>10,000). Furthermore, the behaviour of `progress = TRUE` is
#'    changed if `BP = TRUE`; the default \pkg{BiocParallel} progress bar will
#'    be shown (which unfortunately is much less informative).
#'
#' @return `matrix` if `compare.to` is missing; `data.frame` otherwise.
#' * PCC: 0 represents complete distance, >0 similarity.
#' * MPCC: 0 represents complete distance, 1 complete similarity.
#' * EUCL: 0 represents complete similarity, >0 distance.
#' * MEUCL: 0 represents complete similarity, sqrt(2) complete distance.
#' * SW: 0 represents complete distance, >0 similarity.
#' * MSW: 0 represents complete distance, 2 complete similarity.
#' * KL: 0 represents complete similarity, >0 distance.
#' * MKL: 0 represents complete similarity, >0 complete distance.
#'
#' @details
#' Comparisons are calculated between two motifs at a time. All possible alignments
#' are scored, and the best score is reported. Scores are calculated per position
#' and summed, unless the 'mean' version of the specific metric is chosen. If using
#' a similarity metric, then the sum of scores will favour comparisons between
#' longer motifs; and for distance metrics, the sum of scores will favour
#' comparisons between short motifs. This can be avoided by using the 'mean' of
#' scores.
#'
#' * PCC: Pearson correlation coefficient
#'
#'    Per position:
#'
#'    `PCC = sum(col1 * col2) / sqrt(sum(col1^2) * sum(col2^2))`
#'
#' * MPCC: Mean PCC
#'
#'    `MPCC = mean(PCC)`
#'
#' * EUCL: Euclidian distance
#'
#'    Per position:
#'
#'    `EUCL = sqrt(sum((col1 - col2)^2)) / sqrt(2)`
#'
#' * MEUCL: Mean EUCL
#'
#'    `MEUCL = sum(EUCL) / ncol(alignment)`
#'
#' * SW: Sandelin-Wasserman similarity
#'
#'    Per position:
#'
#'    `SW = 2 - sum((col1 - col2)^2)`
#'
#' * MSW: Mean SW
#'
#'    `MSW = mean(SW)`
#'
#' * KL: Kullback-Leibler divergence
#'
#'    Per position:
#'
#'    `KL = 0.5 * sum(col1 * log(col1/col2) + col2 * log(col2/col1))`
#'
#' * MKL: Mean Kullback-Leibler divergence
#'
#'    `MKL = mean(KL)`
#'
#' To note regarding p-values: p-values are pre-computed using the
#' `make_DBscores` function. If not given, then uses a set of internal
#' precomputed p-values from the JASPAR2018 CORE motifs. These precalculated
#' scores are dependent on the length of the motifs being compared; this takes
#' into account that comparing small motifs with larger motifs leads to higher
#' scores, since the probability of finding a higher scoring alignment is
#' higher.
#'
#' The default p-values have been precalculated for regular DNA motifs; they
#' are of little use for motifs with a different number of alphabet letters
#' (or even the `multifreq` slot).
#'
#' @examples
#' motif1 <- create_motif()
#' motif2 <- create_motif()
#' motif1vs2 <- compare_motifs(list(motif1, motif2), method = "MPCC")
#' ## to get a dist object:
#' as.dist(1 - motif1vs2)
#'
#' @references
#'    \insertRef{euclidean}{universalmotif}
#'
#'    \insertRef{jaspar}{universalmotif}
#'
#'    \insertRef{pearson}{universalmotif}
#'
#'    \insertRef{kldiv}{universalmotif}
#'
#'    \insertRef{wasserman}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [convert_motifs()], [motif_tree()], [view_motifs()]
#' @export
compare_motifs <- function(motifs, compare.to, db.scores, use.freq = 1,
                           use.type = "PPM", method = "MPCC", tryRC = TRUE,
                           min.overlap = 6, min.mean.ic = 0.5,
                           relative_entropy = FALSE, normalise.scores = FALSE,
                           max.p = 0.01, max.e = 10, progress = TRUE,
                           BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  char_check <- check_fun_params(list(use.type = args$use.type,
                                      method = args$method), c(1, 1),
                                 c(FALSE, FALSE), TYPE_CHAR)
  num_check <- check_fun_params(list(compare.to = args$compare.to,
                                     use.freq = args$use.freq,
                                     min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     max.p = args$max.p, max.e = args$max.e),
                                c(0, rep(1, 5)), c(TRUE, rep(FALSE, 5)),
                                TYPE_NUM)
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores,
                                      progress = args$progress, BP = args$BP),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (!use.type %in% c("PPM", "ICM")) {
    type_check <- paste0(" * Incorrect 'use.type': expected `PPM` or `ICM`; ",
                         "got `", use.type, "`")
    all_checks <- c(all_checks, type_check)
  }
  if (!method %in% c("PCC", "MPCC", "EUCL", "MEUCL", "SW", "MSW", "KL", "MKL")) {
    method_check <- paste0(" * 'method': expected `PCC`, `MPCC`, `EUCL`, `MEUCL`",
                           ", `SW`, `MSW`, `KL`, or `MKL`; got `", method, "`")
    method_check <- wmsg2(method_check, 4, 2)
    all_checks <- c(all_checks, method_check)
  }
  if (!missing(db.scores) && !is.data.frame(db.scores)) {
    dbscores_check <- paste0(" * Incorrect type for 'db.scores: expected ",
                             "`data.frame`; got `", class(db.scores), "`")
    all_checks <- c(all_checks, dbscores_check)
  }
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, use.type, relative_entropy = relative_entropy)

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  mot.dup <- mot.names[duplicated(mot.names)]
  if (length(mot.dup) > 0) {
    mot.dup.suffix <- seq_along(mot.dup)
    for (i in seq_along(mot.dup)) {
      mot.dup[i] <- paste0(mot.dup[i], " [duplicated #",
                           mot.dup.suffix[i], "]")
    }
    mot.names[duplicated(mot.names)] <- mot.dup
  }


  if (use.freq == 1) {
    mot.mats <- lapply(motifs, function(x) x@motif)
  } else {
    mot.mats <- lapply(motifs, function(x) x@multifreq[[as.character(use.freq)]])
  }
  mot.ics <- mapply(function(x, y) .pos_iscscores(x, y, relative_entropy),
                      motifs, mot.mats, SIMPLIFY = FALSE)

  if (missing(compare.to)) {

    comparisons <- lapply_(seq_along(motifs)[-length(motifs)],
                            function(x) .compare(x, mot.mats, method,
                                                 min.overlap, tryRC,
                                                 min.mean.ic, mot.ics,
                                                 normalise.scores),
                            PB = progress, BP = BP)

  } else {

    if (missing(db.scores)) {
      if (!normalise.scores) db.scores <- JASPAR2018_CORE_DBSCORES[[method]]
      else db.scores <- JASPAR2018_CORE_DBSCORES_NORM[[method]]
      if (use.freq != 1)
        warning(wmsg("Using the internal P-value database with `use.freq > 1`",
                     " will result in incorrect P-values"),
                immediate. = TRUE)
      mot.alphs <- vapply(motifs, function(x) x["alphabet"], character(1))
      if (!all(mot.alphs == "DNA"))
        warning(wmsg("Using the internal P-value database with non-DNA motifs ",
                     "will result in incorrect P-values"),
                immediate. = TRUE)
    } else {
      db.scores <- check_db_scores(db.scores, method, normalise.scores)
    }

    comparisons <- vector("list", length(compare.to))
    pvals <- vector("list", length(compare.to))

    if (progress) print_pb(0)
    for (i in compare.to) {
      comparisons[[i]] <- lapply_(seq_along(motifs)[-i],
              function(x) motif_simil_internal(mot.mats[[i]],
                                               mot.mats[[x]], method,
                                               min.overlap, tryRC,
                                               mot.ics[[i]],
                                               mot.ics[[x]], min.mean.ic,
                                               normalise.scores), BP = BP)
      comparisons[[i]] <- do.call(c, comparisons[[i]])
      if (!is.null(db.scores)) {
        pvals[[i]] <- pvals_from_db(motifs[compare.to[i]], motifs[-i],
                               db.scores, comparisons[[i]], method)
      } else {
        pvals[[i]] <- rep(NA, length(motifs[-i]))
      }
      if (progress) update_pb(i, length(compare.to))
    }

  }


  if (missing(compare.to)) {

    comparisons <- list_to_matrix_simil(comparisons, mot.names, method)
    
    if (method == "PCC") {
      comparisons <- fix_pcc_diag(comparisons, mot.mats)
    }

  } else {

    comparisons <- lapply(seq_along(compare.to),
                            function(i) {
                         data.frame(subject = rep(mot.names[compare.to[i]],
                                                  length(comparisons[[i]])),
                                    target = mot.names[-compare.to[i]],
                                    score = comparisons[[i]], Pval = pvals[[i]],
                                    Eval = pvals[[i]] * length(motifs) * 2,
                                    stringsAsFactors = FALSE)
                            })
    comparisons <- do.call(rbind, comparisons)
    if (!is.null(db.scores)) {
      comparisons <- comparisons[order(comparisons$Pval, decreasing = FALSE), ]
      comparisons <- comparisons[comparisons$Pval <= max.p &
                                 comparisons$Eval <= max.e, ]
      comparisons <- comparisons[comparisons$subject != comparisons$target, ]
      comparisons <- comparisons[!is.na(comparisons$subject), ]
    }
    comparisons <- get_rid_of_dupes(comparisons)

    if (nrow(comparisons) == 0) {
      message("No significant hits")
      return(invisible(NULL))
    }
    rownames(comparisons) <- NULL

  }

  comparisons

}

check_db_scores <- function(db.scores, method, normalise.scores) {

  if (!is.data.frame(db.scores)) {
    stop("'db.scores' must be a data.frame")
  }

  db_coltypes <- vapply(db.scores, class, character(1))
  db_coltypes_exp <- c("numeric", "numeric", "numeric", "numeric",
                       "character", "logical")
  if (any(db_coltypes != db_coltypes_exp)) {
    stop("'db.scores' must have column types: ",
         paste(db_coltypes_exp, collapse = ", "))
  }

  db_colnames <- colnames(db.scores)
  db_colnames_exp <- c("subject", "target", "mean", "sd", "method", "normalised")
  if (!all(db_colnames_exp %in% db_colnames)) {
    stop(paste0("'db.scores' must have columns: ",
               paste(db_colnames_exp, collapse = ", ")))
  }

  if (!any(method %in% db.scores$method)) {
    stop("could not find method '", method, "' in 'db.scores$method'")
  } else {
    db.scores <- db.scores[db.scores$method == method, ]
  }

  if (!any(normalise.scores %in% db.scores$normalised)) {
    stop("'db.scores' column 'normalised' does not match 'normalise.scores'")
  } else {
    db.scores <- db.scores[db.scores$normalised == normalise.scores, ]
  }

  db.scores

}

fix_pcc_diag <- function(comparisons, mot.mats) {

  mot.lens <- vapply(mot.mats, ncol, numeric(1))
  diag(comparisons) <- mot.lens^2

  comparisons

}

#' get_rid_of_dupes
#'
#' Get rid of entries where comparison is an inverted repeat of a
#' previous comparison. [potential bottleneck?]
#'
#' @param comparisons Data frame of results.
#'
#' @noRd
get_rid_of_dupes <- function(comparisons) {
  comparisons.sorted <- character(nrow(comparisons))
  # for (i in seq_len(nrow(comparisons))) {
    # comparisons.sorted[i] <- paste(sort(comparisons[i, 1:2]), collapse = " ")
  # }
  comparisons.sorted <- vapply(seq_len(nrow(comparisons)),
                               function(x) paste0(sort(comparisons[x, 1:2]),
                                                  collapse = " "),
                               character(1))
  comparisons[!duplicated(comparisons.sorted), ]
} 

#' compare
#'
#' Compare motifs for matrix output. Avoids repeat comparisons.
#'
#' @param x Index to start comparison.
#' @param mot.mats Motif matrices.
#' @param min.overlap Minimum allowed overlap.
#' @param tryRC Compare RC as well.
#' @param min.mean.ic Minimum mean IC for comparison.
#' @param mot.ics Vector of ICs for each position in motifs.
#' @param normalise.scores Normalise scores, logical.
#'
#' @return 1d vector of scores.
#'
#' @noRd
.compare <- function(x, mot.mats, method, min.overlap, tryRC, min.mean.ic,
                     mot.ics, normalise.scores) {
  if (x < length(mot.mats)) x2 <- x + 1 else x2 <- x
  y <- vector("numeric", length = length(mot.mats) - x2)
  index <- 1
  for (j in seq(x2, length(mot.mats))) {
    y[index] <- motif_simil_internal(mot.mats[[x]], mot.mats[[j]], method,
                                 min.overlap, tryRC, mot.ics[[x]],
                                 mot.ics[[j]], min.mean.ic,
                                 normalise.scores)
    index <- index + 1
  }
  y
}

#' pos_icscores
#'
#' Calculate IC for each postion in the motif.
#'
#' @param motif universalmotif class motif.
#' @param mot.mats Motif matrix.
#' @param relative Calculate IC as KL divergence.
#'
#' @noRd
.pos_iscscores <- function(motif, mot.mats, relative = FALSE) {

  bkg <- motif@bkg[rownames(motif@motif)]
  pseudo <- motif@pseudocount
  nsites <- motif@nsites
  if (length(nsites) == 0) nsites <- 100
  apply(mot.mats, 2, function(x) position_icscoreC(x, bkg, "PPM", pseudo,
                                                   nsites, relative))

}

#' @rdname utilities
#' @export
make_DBscores <- function(db.motifs, method, shuffle.db = TRUE,
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
                                method = shuffle.method,
                                leftovers = shuffle.leftovers) 
    if (length(rand.mots) != rand.tries) {
      if (length(rand.mots) < rand.tries) {
        while (length(rand.mots) < rand.tries) {
          more.rand.mots <- shuffle_motifs(db.motifs, k = shuffle.k,
                                           method = shuffle.method,
                                           leftovers = shuffle.leftovers) 
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

#' pvals_from_db
#'
#' Calculate pval for a match from database scores.
#'
#' @param subject Subject motif.
#' @param target Target motifs.
#' @param db DB scores in a data frame.
#' @param scores Scores between subject motif and target motifs.
#' @param method Comparison metric.
#'
#' @return Vector of pvals.
#'
#' @noRd
pvals_from_db <- function(subject, target, db, scores, method) {

  if (method %in% c("PCC", "MPCC", "SW", "MSW")) ltail <- FALSE else ltail <- TRUE
  pvals <- vector("numeric", length(target))
  subject.ncol <- ncol(subject[[1]]@motif)
  if (subject.ncol < db$subject[1]) subject.ncol <- db$subject[1]
  if (subject.ncol > db$subject[nrow(db)]) subject.ncol <- db$subject[nrow(db)]
  possible.ncols <- sort(unique(db$target))
  possible.sub <- sort(unique(db$subject))
  for (i in seq_along(target)) {
    ncol2 <- ncol(target[[i]]@motif)
    if (ncol2 < db$target[1]) ncol2 <- db$target[1]
    if (ncol2 > db$target[nrow(db)]) ncol2 <- db$target[nrow(db)]
    tmp.mean <- db$mean[db$subject == subject.ncol & db$target == ncol2]
    if (length(tmp.mean) == 0) {
      ncol2 <- sort(possible.ncols[possible.ncols <= ncol2])
      ncol2 <- ncol2[length(ncol2)]
      subject.ncol <- sort(possible.sub[possible.sub <= subject.ncol])
      subject.ncol <- subject.ncol[length(subject.ncol)]
      tmp.mean <- db$mean[db$subject == subject.ncol & db$target == ncol2]
    }
    tmp.sd <- db$sd[db$subject == subject.ncol & db$target == ncol2]
    pvals[i] <- pnorm(scores[i], tmp.mean, tmp.sd, ltail)
  }

  pvals

}
