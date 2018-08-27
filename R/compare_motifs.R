#' Compare motifs.
#'
#' Compare motifs using four metrics: Pearson correlation coefficient,
#' Euclidean distance, Sandelin-Wasserman similarity, and Kullback-Leibler
#' divergence.
#'
#' @param motifs See \code{\link{convert_motifs}} for acceptable motif formats.
#' @param compare.to \code{numeric} If missing, compares all motifs to all other motifs.
#'    Otherwise compares all motifs to the specified motif(s).
#' @param db.scores \code{data.frame} See \code{details}.
#' @param use.freq \code{numeric(1)}. For comparing the \code{multifreq} slot.
#' @param use.type \code{character(1)} One of code{'PPM'} and \code{'ICM'}.
#'    The latter allows for taking into account the background
#'    frequencies if \code{relative_entropy = TRUE}.
#' @param method \code{character(1)} One of \code{c('PCC', 'MPCC', 'EUCL', 'MEUCL',
#'    'SW', 'MSW', 'KL', 'MKL')}. See details.
#' @param tryRC \code{logical} Try the reverse complement of the motifs as well,
#'    report the best score.
#' @param min.overlap \code{numeric(1)} Minimum overlap required when aligning the
#'    motifs. Setting this to a number higher then the width of the motifs
#'    will not allow any overhangs. Can also be a number less than 1,
#'    representing the minimum fraction that the motifs must overlap.
#' @param min.mean.ic \code{numeric(1)} Minimum information content between the
#'    two motifs for an alignment to be scored. This helps prevent scoring
#'    alignments between low information content regions of two motifs.
#' @param relative_entropy \code{logical(1)} For ICM calculation. See
#'    \code{\link{convert_type}}.
#' @param normalise.scores \code{logical(1)} Favour alignments which leave fewer
#'    unaligned positions, as well as alignments between motifs of similar length.
#'    Similarity scores are multiplied by the ratio of
#'    aligned positions to the total number of positions in the larger motif,
#'    and the inverse for distance scores.
#' @param max.p \code{numeric(1)} Maximum P-value allowed in reporting matches.
#'    Only used if \code{compare.to} is set.
#' @param max.e \code{numeric(1)} Maximum E-value allowed in reporting matches.
#'    Only used if \code{compare.to} is set. The E-value is the P-value multiplied
#'    by the number of input motifs times two.
#' @param progress \code{logical(1)} Show progress.
#' @param BP \code{logical(1)} Use BiocParallel.
#'
#' @return \code{matrix} if \code{compare.to} is missing; \code{data.frame} otherwise.
#'  \itemize{
#'    \item PCC: 0 represents complete distance, >0 similarity.
#'    \item MPCC: 0 represents complete distance, 1 complete similarity.
#'    \item EUCL: 0 represents complete similarity, >0 distance.
#'    \item MEUCL: 0 represents complete similarity, sqrt(2) complete distance.
#'    \item SW: 0 represents complete distance, >0 similarity.
#'    \item MSW: 0 represents complete distance, 2 complete similarity.
#'    \item KL: 0 represents complete similarity, >0 distance.
#'    \item MKL: 0 represents complete similarity, 4.62 complete distance.
#'  }
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
#' \itemize{
#'   \item PCC: Pearson correlation coefficient
#'   \item MPCC: Mean PCC
#'   \item EUCL: Euclidian distance
#'   \item MEUCL: Mean EUCL
#'   \item SW: Sandelin-Wasserman similarity
#'   \item MSW: Mean SW
#'   \item KL: Kullback-Leibler divergence
#'   \item MKL: Mean Kullback-Leibler divergence
#' }
#'
#' To note regarding p-values: p-values are pre-computed using the
#' \code{make_DBscores} function. If not given, then uses a set of internal
#' precomputed p-values from the JASPAR2018 CORE motifs. These precalculated
#' scores are dependent on the length of the motifs being compared; this takes
#' into account that comparing small motifs with larger motifs leads to higher
#' scores, since the probability of finding a higher scoring alignment is
#' higher.
#'
#' The default p-values have been precalculated for regular DNA motifs; they
#' are of little use for motifs with a different number of alphabet letters
#' (or even the \code{multifreq} slot).
#'
#' @examples
#' motif1 <- create_motif()
#' motif2 <- create_motif()
#' motif1vs2 <- compare_motifs(list(motif1, motif2), method = "MPCC")
#' ## to get a dist object:
#' as.dist(1 - motif1vs2)
#'
#' @references
#'    \insertRef{jaspar}{universalmotif}
#'
#'    \insertRef{tfbstools}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso \code{\link{convert_motifs}}, \code{\link[TFBSTools]{PWMSimilarity}},
#'    \code{\link{motif_tree}}, \code{\link{view_motifs}}
#' @export
compare_motifs <- function(motifs, compare.to, db.scores, use.freq = 1,
                           use.type = "PPM", method = "MPCC", tryRC = TRUE,
                           min.overlap = 6, min.mean.ic = 0.5,
                           relative_entropy = FALSE, normalise.scores = FALSE,
                           max.p = 0.01, max.e = 10, progress = TRUE,
                           BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(use.type = args$use.type,
                                      method = args$method), c(1, 1),
                                 c(FALSE, FALSE), "character")
  num_check <- check_fun_params(list(compare.to = args$compare.to,
                                     use.freq = args$use.freq,
                                     min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     max.p = args$max.p, max.e = args$max.e),
                                c(0, rep(1, 5)), c(TRUE, rep(FALSE, 5)),
                                "numeric")
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores,
                                      progress = args$progress, BP = args$BP),
                                 numeric(), logical(), "logical")
  if (!use.type %in% c("PPM", "ICM")) {
    type_check <- paste0(" * Incorrect 'use.type': expected `PPM` or `ICM`; ",
                         "got `", use.type, "`")
    all_checks <- c(all_checks, type_check)
  }
  if (!method %in% c("PCC", "MPCC", "EUCL", "MEUCL", "SW", "MSW", "KL", "MKL")) {
    method_check <- paste0(" * 'method': expected `PCC`, `MPCC`, `EUCL`, `MEUCL`",
                           ", `SW`, `MSW`, `KL`, or `MKL`; got `", method, "`")
    all_checks <- c(all_checks, method_check)
  }
  all_checks <- c(char_check, num_check, logi_check)
  if (!missing(db.scores) && !is.data.frame(db.scores)) {
    dbscores_check <- paste0(" * Incorrect type for 'db.scores: expected ",
                             "`data.frame`; got `", class(db.scores), "`")
    all_checks <- c(all_checks, dbscores_check)
  }
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, use.type, relative_entropy = relative_entropy)

  mot.names <- vapply(motifs, function(x) x["name"], character(1))
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
    mot.mats <- lapply(motifs, function(x) x["motif"])
  } else {
    mot.mats <- lapply(motifs, function(x) x["multifreq"][[as.character(use.freq)]])
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
  for (i in seq_along(mot.lens)) {
    comparisons[i, i] <- mot.lens[i]^2
  }

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
  for (i in seq_len(nrow(comparisons))) {
    comparisons.sorted[i] <- paste(sort(comparisons[i, 1:2]), collapse = " ")
  }
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

  bkg <- motif["bkg"]
  pseudo <- motif["pseudocount"]
  nsites <- motif["nsites"]
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
                                 numeric(), logical(), "character")
  num_check <- check_fun_params(list(shuffle.k = args$shuffle.k,
                                     rand.tries = args$rand.tries,
                                     min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic),
                                numeric(), logical(), "numeric")
  logi_check <- check_fun_params(list(shuffle.db = args$shuffle.db,
                                      progress = args$progress, BP = args$BP,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  db.motifs <- convert_motifs(db.motifs)
  db.ncols <- vapply(db.motifs, function(x) ncol(x["motif"]), numeric(1))

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
                        function(x) create_motif(sample(5:30, 1)))
  }
  rand.ncols <- vapply(rand.mots, function(x) ncol(x["motif"]), numeric(1))

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
  subject.ncol <- ncol(subject[[1]]["motif"])
  if (subject.ncol < db$subject[1]) subject.ncol <- db$subject[1]
  if (subject.ncol > db$subject[nrow(db)]) subject.ncol <- db$subject[nrow(db)]
  possible.ncols <- sort(unique(db$target))
  possible.sub <- sort(unique(db$subject))
  for (i in seq_along(target)) {
    ncol2 <- ncol(target[[i]]["motif"])
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
