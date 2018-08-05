#' Compare motifs.
#'
#' TFBSTools implementation of motif comparison, rewritten in C++. See
#' \code{\link[TFBSTools]{PWMSimilarity}} for details.
#'
#' @param motifs List of motifs. If not a \linkS4class{universalmotif} object,
#'    they will be converted. See \code{\link{convert_motifs}} for supported
#'    classes.
#' @param compare.to Numeric. If NULL, compares all motifs to all other motifs.
#'    Otherwise compares all motifs to the specified motif.
#' @param db.scores data.frame.
#' @param method One of 'Euclidean', 'Pearson', and 'KL'.
#' @param tryRC Try the reverse complement of the motifs as well, report the
#'    better score.
#' @param min.overlap Numeric. Minimum overlap required when aligning the
#'    motifs. Setting this to a number higher then the width of the motifs
#'    will not allow any overhangs.
#' @param min.mean.ic Numeric.
#' @param relative_entropy Logical.
#' @param max.p Numeric.
#' @param max.e Numeric.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Distance matrix for Euclidean or KL; similarity matrix for Pearson.
#'
#' @details
#'    The implementations of this function are the exact same as that of the
#'    \code{\link[TFBSTools]{PWMSimilarity}} function, except that this function
#'    allows for more than two motifs to be compared as well as providing
#'    significant performance gains.
#'
#'    Each possible motif pairs are compared to generate the final matrix. This
#'    is done by first aligning the two motifs, then performing the
#'    distance/similarity calculation for each position pairs. If the two
#'    motifs are not the same length, then the calculation is performed
#'    repeatedly by moving the smaller motif along the larger motif. Afterwards,
#'    either the smallest distance or the largest similarity is reported.
#'
#' @examples
#' motif1 <- create_motif()
#' motif2 <- create_motif()
#' motif1vs2 <- compare_motifs(list(motif1, motif2))
#' # to get a dist object:
#' as.dist(motif1vs2)
#'
#' @references
#'    \insertRef{tfbstools}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso \code{\link{convert_motifs}}, \code{\link[TFBSTools]{PWMSimilarity}},
#'    \code{\link{motif_tree}}
#' @export
compare_motifs <- function(motifs, compare.to, db.scores, method = "Pearson",
                           tryRC = TRUE, min.overlap = 6, min.mean.ic = 0.5,
                           relative_entropy = FALSE, max.p = 0.05, max.e = 10,
                           BPPARAM = SerialParam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  if (method == "KL") {
    motifs <- lapply(motifs, function(x) {
                             x["pseudocount"] <- x["pseudocount"] + 0.0001
                             x
                           })
  }
  motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)

  mot.names <- vapply(motifs, function(x) x["name"], character(1))

  mot.mats <- lapply(motifs, function(x) x["motif"])
  mot.ics <- bpmapply(function(x, y) .pos_iscscores(x, y, relative_entropy),
                      motifs, mot.mats, BPPARAM = BPPARAM)
  if (missing(compare.to)) {
    comparisons <- bplapply(seq_along(motifs)[-length(motifs)],
                            function(x) .compare(x, mot.mats, method,
                                                 min.overlap, tryRC,
                                                 min.mean.ic, mot.ics),
                            BPPARAM = BPPARAM)
  } else {
    comparisons <- vector("list", length(compare.to))
    pvals <- vector("list", length(compare.to))
    for (i in compare.to) {
      comparisons[[i]] <- bplapply(seq_along(motifs)[-compare.to],
              function(x) motif_simil_internal(mot.mats[[compare.to[i]]],
                                               mot.mats[[x]], method,
                                               min.overlap, tryRC,
                                               mot.ics[[compare.to[i]]],
                                               mot.ics[[x]],
                                               min.mean.ic),
                                   BPPARAM = BPPARAM)
      comparisons[[i]] <- do.call(c, comparisons[[i]])
      if (!missing(db.scores)) {
        pvals[[i]] <- pvals_from_db(motifs[compare.to[i]], motifs[-compare.to],
                               db.scores, comparisons[[i]], method)
      }
    }
  }

  if (missing(compare.to)) {
    comparisons <- list_to_matrix_simil(comparisons, mot.names, method)
  } else {
    names(comparisons) <- mot.names[compare.to]
    if (!missing(db.scores)) {
      comparisons <- bplapply(seq_along(compare.to),
                              function(i) {
                           data.frame(subject = rep(mot.names[compare.to[i]],
                                      length(comparisons[[i]])),
                                      target = mot.names[-compare.to],
                                      score = comparisons[[i]], Pval = pvals[[i]],
                                      Eval = pvals[[i]] * length(motifs) * 2,
                                      stringsAsFactors = FALSE)
                              }, BPPARAM = BPPARAM)
    } else {
      comparisons <- bplapply(seq_along(compare.to),
                              function(x) {
                                data.frame(subject = rep(mot.names[compare.to[x]],
                                                         length(comparisons[[x]])),
                                           target = mot.names[-compare.to],
                                           score = comparisons[[x]],
                                           stringsAsFactors = FALSE)
                              }, BPPARAM = BPPARAM)
    }
    comparisons <- do.call(rbind, comparisons)
    if (!missing(db.scores)) {
      comparisons <- comparisons[order(comparisons$Pval, decreasing = FALSE), ]
      comparisons <- comparisons[comparisons$Pval <= max.p &
                                 comparisons$Eval <= max.e, ]
      comparisons <- comparisons[comparisons$subject != comparisons$target, ]
      if (nrow(comparisons) == 0) {
        message("No significant hits")
        return(invisible(NULL))
      }
    } else {
      if (method == "Pearson") order.logi <- TRUE else order.logi <- FALSE
      comparisons <- comparisons[order(comparisons$score,
                                       decreasing = order.logi), ]
      comparisons <- comparisons[comparisons$subject != comparisons$target, ]
    }
    rownames(comparisons) <- NULL
  }

  comparisons

}

.compare <- function(x, mot.mats, method, min.overlap, tryRC, min.mean.ic,
                     mot.ics) {
  if (x < length(mot.mats)) x2 <- x + 1 else x2 <- x
  y <- vector("numeric", length = length(mot.mats) - x2)
  index <- 1
  for (j in seq(x2, length(mot.mats))) {
    y[index] <- motif_simil_internal(mot.mats[[x]], mot.mats[[j]], method,
                                 min.overlap, tryRC, mot.ics[[x]],
                                 mot.ics[[j]], min.mean.ic)
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
make_DBscores <- function(db.motifs, method = "Pearson", tries = 1000,
                          min.overlap = 6, progress_bar = TRUE,
                          BPPARAM = SerialParam()) {

  db.ncols <- vapply(db.motifs, function(x) ncol(x["motif"]), numeric(1))
  db.dist <- table(db.ncols)

  rand.mots <- lapply(seq_len(tries),
                      function(x) create_motif(sample(5:30, 1)))
  rand.ncols <- vapply(rand.mots, function(x) ncol(x["motif"]), numeric(1))
  rand.dist <- table(rand.ncols)

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
                               min.overlap = min.overlap, min.mean.ic = 0,
                               BPPARAM = BPPARAM)$score

    totry$mean[i] <- mean(res[[i]])
    totry$sd[i] <- sd(res[[i]])

    if (progress_bar) setTxtProgressBar(pb, i)

  }

  if (progress_bar) close(pb)

  totry$method <- rep(method, nrow(totry))
  totry$db <- rep(deparse(substitute(db.motifs)), nrow(totry))
  totry

}

pvals_from_db <- function(subject, target, db, scores, method) {

  if (method == "Pearson") ltail <- FALSE else ltail <- TRUE
  pvals <- vector("numeric", length(target))
  subject.ncol <- ncol(subject[[1]]["motif"])
  possible.ncols <- sort(unique(db$target))
  for (i in seq_along(target)) {
    ncol2 <- ncol(target[[i]]["motif"])
    tmp.mean <- db$mean[db$subject == subject.ncol & db$target == ncol2]
    if (length(tmp.mean) == 0) {
      ncol2 <- possible.ncols[possible.ncols < ncol2]
      ncol2 <- ncol2[length(ncol2)]
      tmp.mean <- db$mean[db$subject == subject.ncol & db$target == ncol2]
    }
    tmp.sd <- db$sd[db$subject == subject.ncol & db$target == ncol2]
    pvals[i] <- pnorm(scores[i], tmp.mean, tmp.sd, ltail)
  }

  pvals

}
