#' Compare motifs.
#'
#' Compare motifs using three metrics: Pearson correlation coefficient,
#' Euclidean distance, and Kullback-Leibler divergence.
#'
#' @param motifs See \code{\link{convert_motifs}} for acceptable motif formats.
#' @param compare.to \code{numeric} If missing, compares all motifs to all other motifs.
#'    Otherwise compares all motifs to the specified motif(s).
#' @param db.scores \code{data.frame} See \code{details}.
#' @param use.freq \code{numeric(1)}. For comparing the \code{multifreq} slot.
#' @param use.type \code{character(1)} One of \code{'PCM'} (Pearson only),
#'    \code{'PPM'} (any method),
#'    \code{'PWM'} (Pearson only), and \code{'ICM'} (any method). The latter
#'    two allow for taking into account the background frequencies 
#'    (for ICM, only if \code{relative_entropy = TRUE}).
#' @param method \code{character(1)} One of \code{'Pearson', 'Euclidean', 'KL'}.
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
#'    unaligned positions.
#' @param max.p \code{numeric(1)} Maximum P-value allowed in reporting matches.
#'    Only used if \code{compare.to} is set.
#' @param max.e \code{numeric(1)} Maximum E-value allowed in reporting matches.
#'    Only used if \code{compare.to} is set. The E-value is the P-value multiplied
#'    by the number of input motifs times two.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return \code{matrix} if \code{compare.to} is missing; \code{data.frame} otherwise.
#'    If \code{method = c('Euclidean', 'KL')} then the resulting scores represent
#'    distances; for \code{method = 'Pearson'}, similarities.
#'
#' @details
#' The comparison metrics are identical to those implemented by 
#' \code{\link[TFBSTools]{PWMSimilarity}}, rewritten in C++.
#'
#' To note regarding p-values: p-values are pre-computed using the
#' \code{make_DBscores} function. If not given, then uses a set of internal
#' precomputed p-values from the JASPAR2018 CORE motifs. Furthermore, the
#' comparison calculation does not take into account the difference in length
#' between motifs; this leads to an increased likelihood in higher scores
#' between small and large motifs. In order to overcome this limitation,
#' the p-values are calculated with this in mind.
#'
#' The default p-values have been precalculated for regular DNA motifs; they
#' are of little use for motifs with a different number of alphabet letters
#' (or even the \code{multifreq} slot).
#'
#' @examples
#' motif1 <- create_motif()
#' motif2 <- create_motif()
#' motif1vs2 <- compare_motifs(list(motif1, motif2), method = "Pearson")
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
compare_motifs <- function(motifs, compare.to, db.scores, use.freq = 1, use.type = "PPM",
                           method = "Pearson", tryRC = TRUE, min.overlap = 6,
                           min.mean.ic = 0.5, relative_entropy = FALSE,
                           normalise.scores = TRUE, max.p = 0.05, max.e = 10,
                           BPPARAM = SerialParam()) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(use.type = args$use.type, method = args$method),
                                 c(1, 1), c(FALSE, FALSE), "character")
  num_check <- check_fun_params(list(compare.to = args$compare.to,
                                     use.freq = args$use.freq,
                                     min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     max.p = args$max.p, max.e = args$max.e),
                                c(0, rep(1, 5)), c(TRUE, rep(FALSE, 5)),
                                "numeric")
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores),
                                 c(1, 1, 1), c(FALSE, FALSE, FALSE), "logical")
  s4_check <- check_fun_params(list(BPPARAM = args$BPPARAM), numeric(), FALSE, "S4")
  all_checks <- c(char_check, num_check, logi_check, s4_check)
  if (!missing(db.scores) && !is.data.frame(db.scores)) {
    dbscores_check <- paste0(" * Incorrect type for 'db.scores: expected ",
                             "`data.frame`; got `", class(db.scores), "`")
    all_checks <- c(all_checks, dbscores_check)
  }
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------


  if (use.type %in% c("PCM", "PWM") && method %in% c("Euclidean", "KL")) {
    stop("Method '", method, "' is not supported for type '", use.type, "'")
  }

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, use.type, relative_entropy = relative_entropy,
                         BPPARAM = BPPARAM)

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
  mot.ics <- bpmapply(function(x, y) .pos_iscscores(x, y, relative_entropy),
                      motifs, mot.mats, BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  if (missing(compare.to)) {

    comparisons <- bplapply(seq_along(motifs)[-length(motifs)],
                            function(x) .compare(x, mot.mats, method,
                                                 min.overlap, tryRC,
                                                 min.mean.ic, mot.ics,
                                                 normalise.scores),
                            BPPARAM = BPPARAM)

  } else {

    if (missing(db.scores)) {
      if (!normalise.scores) db.scores <- JASPAR2018_CORE_DBSCORES[[method]]
      else db.scores <- JASPAR2018_CORE_DBSCORES_NORM[[method]]
    }

    comparisons <- vector("list", length(compare.to))
    pvals <- vector("list", length(compare.to))

    for (i in compare.to) {
      comparisons[[i]] <- bplapply(seq_along(motifs)[-i],
              function(x) motif_simil_internal(mot.mats[[i]],
                                               mot.mats[[x]], method,
                                               min.overlap, tryRC,
                                               mot.ics[[i]],
                                               mot.ics[[x]], min.mean.ic,
                                               normalise.scores),
                                   BPPARAM = BPPARAM)
      comparisons[[i]] <- do.call(c, comparisons[[i]])
      pvals[[i]] <- pvals_from_db(motifs[compare.to[i]], motifs[-i],
                             db.scores, comparisons[[i]], method)
    }

  }

  if (missing(compare.to)) {

    comparisons <- list_to_matrix_simil(comparisons, mot.names, method)

  } else {

    comparisons <- bplapply(seq_along(compare.to),
                            function(i) {
                         data.frame(subject = rep(mot.names[compare.to[i]],
                                                  length(comparisons[[i]])),
                                    target = mot.names[-compare.to[i]],
                                    score = comparisons[[i]], Pval = pvals[[i]],
                                    Eval = pvals[[i]] * length(motifs) * 2,
                                    stringsAsFactors = FALSE)
                            }, BPPARAM = BPPARAM)
    comparisons <- do.call(rbind, comparisons)
    comparisons <- comparisons[order(comparisons$Pval, decreasing = FALSE), ]
    comparisons <- comparisons[comparisons$Pval <= max.p &
                               comparisons$Eval <= max.e, ]
    comparisons <- comparisons[comparisons$subject != comparisons$target, ]
    comparisons <- comparisons[!is.na(comparisons$subject), ]
    comparisons <- get_rid_of_dupes(comparisons)

    if (nrow(comparisons) == 0) {
      message("No significant hits")
      return(invisible(NULL))
    }
    rownames(comparisons) <- NULL

  }

  comparisons

}

# potential bottleneck
get_rid_of_dupes <- function(comparisons) {
  comparisons.sorted <- apply(comparisons[, 1:2], 1,
                              function(x) paste(sort(x), collapse = " "))
  comparisons[!duplicated(comparisons.sorted), ]
} 

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
make_DBscores <- function(db.motifs, method = "Pearson", shuffle.db = TRUE,
                          shuffle.k = 3, shuffle.method = "linear",
                          shuffle.leftovers = "asis", rand.tries = 1000,
                          normalise.scores = TRUE,
                          min.overlap = 6, min.mean.ic = 0, progress_bar = TRUE,
                          BPPARAM = SerialParam()) {

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
                                      progress_bar = args$progress_bar,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), "logical")
  s4_check <- check_fun_params(list(BPPARAM = args$BPPARAM),
                               numeric(), logical(), "S4")
  all_checks <- c(char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  db.ncols <- vapply(db.motifs, function(x) ncol(x["motif"]), numeric(1))

  if (shuffle.db) {
    rand.mots <- shuffle_motifs(db.motifs, k = shuffle.k,
                                method = shuffle.method,
                                leftovers = shuffle.leftovers,
                                BPPARAM = BPPARAM) 
    if (length(rand.mots) != rand.tries) {
      if (length(rand.mots) < rand.tries) {
        while (length(rand.mots) < rand.tries) {
          more.rand.mots <- shuffle_motifs(db.motifs, k = shuffle.k,
                                           method = shuffle.method,
                                           leftovers = shuffle.leftovers,
                                           BPPARAM = BPPARAM) 
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

  if (progress_bar) pb <- txtProgressBar(max = length(res), style = 3)

  for (i in seq_len(nrow(totry))) {

    tmp1 <- db.motifs[totry[i, 2] == db.ncols]
    tmp2 <- rand.mots[totry[i, 1] == rand.ncols]

    res[[i]] <- compare_motifs(c(tmp2, tmp1), seq_along(tmp2), method = method,
                               min.overlap = min.overlap, min.mean.ic = min.mean.ic,
                               BPPARAM = BPPARAM, max.e = Inf, max.p = Inf,
                               normalise.scores = normalise.scores)$score

    totry$mean[i] <- mean(res[[i]])
    totry$sd[i] <- sd(res[[i]])

    if (progress_bar) setTxtProgressBar(pb, i)

  }

  if (progress_bar) close(pb)

  totry$method <- rep(method, nrow(totry))
  totry$db <- rep(deparse(substitute(db.motifs)), nrow(totry))
  totry$normalised <- rep(normalise.scores, nrow(totry))
  totry

}

pvals_from_db <- function(subject, target, db, scores, method) {

  if (method == "Pearson") ltail <- FALSE else ltail <- TRUE
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
