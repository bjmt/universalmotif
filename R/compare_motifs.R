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
#' @param min.position.ic `numeric(1)` Minimum information content required between
#'    individual alignment positions for it to be counted in the final alignment
#'    score. It is recommended to use this together with `normalise.scores = TRUE`,
#'    as this will help punish scores resulting from only a fraction of an
#'    alignment.
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
#' @param progress `logical(1)` Deprecated. Does nothing.
#' @param BP `logical(1)` Deprecated. See `nthreads`.
#' @param nthreads `numeric(1)` Run [compare_motifs()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#' @param db.version `numeric(1)` Select internal P-value database. `1` is the
#'    database created in v1.0.0 of the package, and `2` is the database
#'    created in v1.4.0 of the package.
#'
#' @return `matrix` if `compare.to` is missing; `data.frame` otherwise.
#' * PCC: The number of positions in the shorter motif represents the max possible
#'   score. Values below that represent dissimilarity.
#' * MPCC: 1 represents complete similarity, <1 dissimilarity.
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
#'    `PCC = (nrow * sum(pos1 * pos2) - sum(pos1) * sum(pos2)) /
#'           sqrt((nrow * sum(pos1)^2 - sum(pos1^2)) *
#'                (nrow * sum(pos2)^2 - sum(pos2^2)))`
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
#' @seealso [convert_motifs()], [motif_tree()], [view_motifs()],
#'    [make_DBscores()]
#' @export
compare_motifs <- function(motifs, compare.to, db.scores, use.freq = 1,
                           use.type = "PPM", method = "MPCC", tryRC = TRUE,
                           min.overlap = 6, min.mean.ic = 0.25,
                           min.position.ic = 0,
                           relative_entropy = FALSE, normalise.scores = FALSE,
                           max.p = 0.01, max.e = 10, progress = FALSE,
                           BP = FALSE, nthreads = 1, db.version = 2) {

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
                                     max.p = args$max.p, max.e = args$max.e,
                                     nthreads = args$nthreads,
                                     min.position.ic = args$min.position.ic,
                                     db.version = args$db.version),
                                c(0, rep(1, 7)), c(TRUE, rep(FALSE, 7)),
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
  if (!method %in% COMPARE_METRICS) {
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

  if (!db.version %in% 1:2)
    stop("db.scores can only be 1 or 2")

  if (progress)
    warning("'progress' is deprecated and does nothing", immediate. = TRUE)
  if (BP)
    warning("'BP' is deprecated; use 'nthreads' instead", immediate. = TRUE)

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

  check.nrow <- vapply(mot.mats, nrow, integer(1))
  if (length(unique(check.nrow)) > 1)
    stop("all motifs must have the same number of rows")

  alph <- unique(vapply(motifs, function(x) x@alphabet, character(1)))
  if (length(alph) > 1) stop("all motifs must have the same alphabet")

  alph <- switch(alph, "DNA" = DNA_BASES, "RNA" = RNA_BASES, "AA" = AA_STANDARD,
                 sort_unique_cpp(safeExplode(alph)))

  mot.type <- switch(use.type, "PPM" = 1, "ICM" = 2)
  mot.bkgs <- lapply(motifs, function(x) x@bkg[seq_along(alph)])

  if (missing(compare.to)) {

    comparisons <- compare_motifs_all(mot.mats, method, min.overlap,
                                      tryRC, min.mean.ic, normalise.scores,
                                      nthreads, mot.bkgs, mot.type,
                                      relative_entropy, mot.names,
                                      min.position.ic)

    if (method == "PCC") comparisons <- fix_pcc_diag(comparisons, mot.mats)

  } else {

    if (missing(db.scores)) {

      if (db.version == 2) {
        if (!normalise.scores)
          db.scores <- JASPAR2018_CORE_DBSCORES_2[[method]]
        else
          db.scores <- JASPAR2018_CORE_DBSCORES_NORM_2[[method]]
      } else if (db.version == 1) {
        if (!normalise.scores)
          db.scores <- JASPAR2018_CORE_DBSCORES[[method]]
        else
          db.scores <- JASPAR2018_CORE_DBSCORES_NORM[[method]]
      }

      if (use.freq != 1)
        warning(wmsg("Using the internal P-value database with `use.freq > 1`",
                     " will result in incorrect P-values"),
                immediate. = TRUE)
      mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
      if (!all(mot.alphs == "DNA"))
        warning(wmsg("Using the internal P-value database with non-DNA motifs ",
                     "will result in incorrect P-values"),
                immediate. = TRUE)

    } else {
      db.scores <- check_db_scores(db.scores, method, normalise.scores)
    }

    comps <- get_comp_indices(compare.to, length(motifs))
    comparisons <- compare_motifs_cpp(mot.mats, comps[, 1] - 1, comps[, 2] - 1,
                                      method, min.overlap, tryRC, mot.bkgs,
                                      mot.type, relative_entropy,
                                      min.mean.ic, normalise.scores, nthreads,
                                      min.position.ic)

    pvals <- vector("list", length(compare.to))
    for (i in seq_along(compare.to)) {
      pvals[[i]] <- pvals_from_db(motifs[compare.to[i]], motifs[-compare.to[i]],
                                  db.scores, comparisons[comps[, 1] == compare.to[i]],
                                  method)
    }
    pvals <- unlist(pvals)

    comparisons <- data.frame(subject = mot.names[comps[, 1]],
                              target = mot.names[comps[, 2]],
                              score = comparisons,
                              logPval = pvals,
                              Pval = exp(1)^pvals,
                              Eval = exp(1)^pvals * length(motifs) * 2,
                              stringsAsFactors = FALSE)

    if (!is.null(db.scores)) {
      comparisons <- comparisons[order(comparisons$logPval, decreasing = FALSE), ]
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
         paste0(db_coltypes_exp, collapse = ", "))
  }

  db_colnames <- colnames(db.scores)
  db_colnames_exp <- c("subject", "target", "mean", "sd", "method", "normalised")
  if (!all(db_colnames_exp %in% db_colnames)) {
    stop(paste0("'db.scores' must have columns: ",
               paste0(db_colnames_exp, collapse = ", ")))
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

get_comp_indices <- function(compare.to, n) {

  out <- vector("list", length(compare.to))
  for (i in seq_along(compare.to)) {
    out[[i]] <- data.frame(rep(compare.to[i], n - 1), seq_len(n)[-compare.to[i]])
  }

  do.call(rbind, out)

}

fix_pcc_diag <- function(comparisons, mot.mats) {

  mot.lens <- vapply(mot.mats, ncol, numeric(1))
  diag(comparisons) <- mot.lens

  comparisons

}

get_rid_of_dupes <- function(comparisons) {
  comparisons.sorted <- character(nrow(comparisons))
  comparisons.sorted <- vapply(seq_len(nrow(comparisons)),
                               function(x) paste0(sort(comparisons[x, 1:2]),
                                                  collapse = " "),
                               character(1))
  comparisons[!duplicated(comparisons.sorted), ]
} 

compare_motifs_all <- function(mot.mats, method, min.overlap, RC, min.mean.ic,
                               normalise.scores, nthreads, mot.bkgs,
                               mot.type, relative, mot.names,
                               min.position.ic) {

  ans <- compare_motifs_all_cpp(mot.mats, method, min.overlap, RC, mot.bkgs,
                                mot.type, relative, min.mean.ic,
                                normalise.scores, nthreads, min.position.ic)

  n <- length(mot.mats)
  comp <- comb2_cpp(n)
  get_comparison_matrix(unlist(ans), comp[[1]], comp[[2]], method, mot.names)

}

pvals_from_db <- function(subject, target, db, scores, method) {

  if (method %in% c("PCC", "MPCC", "SW", "MSW"))
    ltail <- FALSE
  else
    ltail <- TRUE

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
    pvals[i] <- pnorm(scores[i], tmp.mean, tmp.sd, ltail, TRUE)

  }

  pvals

}
