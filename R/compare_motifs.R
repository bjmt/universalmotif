#' Compare motifs.
#'
#' Compare motifs using any of several available metrics. See the
#' "Motif comparisons and P-values" vignette for detailed information.
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param compare.to `numeric` If missing, compares all motifs to all other motifs.
#'    Otherwise compares all motifs to the specified motif(s).
#' @param db.scores `data.frame` or `DataFrame`. See `details`.
#' @param use.freq `numeric(1)`. For comparing the `multifreq` slot.
#' @param use.type `character(1)` One of `'PPM'` and `'ICM'`.
#'    The latter allows for taking into account the background
#'    frequencies if `relative_entropy = TRUE`. Note that `'ICM'` is not
#'    allowed when `method = c(ALLR, ALLR_LL)`.
#' @param method `character(1)` One of PCC, EUCL, SW, KL, ALLR, BHAT, HELL, IS,
#'    SEUCL, MAN, ALLR_LL. See details.
#' @param tryRC `logical(1)` Try the reverse complement of the motifs as well,
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
#' @param relative_entropy `logical(1)` Change the ICM calculation affecting
#'    `min.position.ic` and `min.mean.ic`. See [convert_type()].
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
#' @param score.strat `character(1)` How to handle column scores calculated from
#'    motif alignments. "sum": add up all scores. "a.mean": take the arithmetic
#'    mean. "g.mean": take the geometric mean. "median": take the median.
#'    "wa.mean", "wg.mean": weighted arithmetic/geometric mean. Weights are the
#'    total information content shared between aligned columns. 
#' @param output.report `character(1)` Provide a filename for [compare_motifs()]
#'    to write an html ouput report to. The top matches are shown alongside
#'    figures of the match alignments. This requires the `knitr` and `rmarkdown`
#'    packages.
#' @param output.report.max.print `numeric(1)` Maximum number of top matches to
#'    print.
#'
#' @return `matrix` if `compare.to` is missing; `DataFrame` otherwise. For the
#'    latter, function args are stored in the `metadata` slot.
#'
#' @details
#' The following metrics are available:
#'
#' * Euclidean distance (`EUCL`) \insertCite{euclidean}{universalmotif}
#' * Kullback-Leibler divergence (`KL`) \insertCite{kl,kldiv}{universalmotif}
#' * Hellinger distance (`HELL`) \insertCite{hellinger}{universalmotif}
#' * Itakura-Saito distance (`IS`) \insertCite{ISdist}{universalmotif}
#' * Squared Euclidean distance (`SEUCL`)
#' * Manhattan distance (`MAN`)
#' * Pearson correlation coefficient (`PCC`)
#' * Sandelin-Wasserman similarity (`SW`), or sum of squared distances \insertCite{wasserman}{universalmotif}
#' * Average log-likelihood ratio (`ALLR`) \insertCite{wang}{universalmotif}
#' * Lower limit ALLR (`ALLR_LL`) \insertCite{mahony}{universalmotif}
#' * Bhattacharyya coefficient (`BHAT`) \insertCite{bhatt}{universalmotif}
#'
#' Comparisons are calculated between two motifs at a time. All possible alignments
#' are scored, and the best score is reported. In an alignment scores are calculated
#' individually between columns. How those scores are combined to generate the final
#' alignment scores depends on `score.strat`.
#'
#' See the "Motif comparisons and P-values" vignette for a description of the
#' various metrics. Note that PCC, SW, ALLR, ALLR_LL and BHAT are similarity;
#' higher values mean more similar motifs. For the remaining metrics, values closer
#' to zero represent more similar motifs.
#'
#' Small pseudocounts are automatically added when one of the following methods
#' is used: KL, ALLR, IS, ALLR_LL. This is avoid
#' zeros in the calculations.
#'
#' To note regarding p-values: P-values are pre-computed using the
#' [make_DBscores()] function. If not given, then uses a set of internal
#' precomputed P-values from the JASPAR2018 CORE motifs. These precalculated
#' scores are dependent on the length of the motifs being compared; this takes
#' into account that comparing small motifs with larger motifs leads to higher
#' scores, since the probability of finding a higher scoring alignment is
#' higher.
#'
#' The default P-values have been precalculated for regular DNA motifs; they
#' are of little use for motifs with a different number of alphabet letters
#' (or even the `multifreq` slot).
#'
#' @examples
#' motif1 <- create_motif()
#' motif2 <- create_motif()
#' motif1vs2 <- compare_motifs(list(motif1, motif2), method = "PCC")
#' ## to get a dist object:
#' as.dist(1 - motif1vs2)
#'
#' @references
#'
#'    \insertRef{bhatt}{universalmotif}
#'
#'    \insertRef{euclidean}{universalmotif}
#'
#'    \insertRef{hellinger}{universalmotif}
#'
#'    \insertRef{jaspar}{universalmotif}
#'
#'    \insertRef{kl}{universalmotif}
#'
#'    \insertRef{ISdist}{universalmotif}
#'
#'    \insertRef{mahony}{universalmotif}
#'
#'    \insertRef{pearson}{universalmotif}
#'
#'    \insertRef{kldiv}{universalmotif}
#'
#'    \insertRef{wasserman}{universalmotif}
#'
#'    \insertRef{wang}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [convert_motifs()], [motif_tree()], [view_motifs()],
#'    [make_DBscores()]
#' @export
compare_motifs <- function(motifs, compare.to, db.scores, use.freq = 1,
                           use.type = "PPM", method = "ALLR", tryRC = TRUE,
                           min.overlap = 6, min.mean.ic = 0.25,
                           min.position.ic = 0,
                           relative_entropy = FALSE, normalise.scores = FALSE,
                           max.p = 0.01, max.e = 10, progress = FALSE,
                           BP = FALSE, nthreads = 1,
                           score.strat = "a.mean", output.report,
                           output.report.max.print = 10) {

  # Some timings:
  #
  # ~10,000 comparisons:
  #    compare_motifs(MotifDb, 1, nthreads = 6)          2.488 sec
  #
  # ~100,000 comparisons:
  #    compare_motifs(MotifDb, 1:10, nthreads = 6)       5.558 sec
  #
  # ~1,000,000 comparisons:
  #    compare_motifs(MotifDb, 1:100, nthreads = 6)     48.941 sec
  #
  # ~10,000,000 comparisons:
  #    compare_motifs(MotifDb, 1:1000, nthreads = 6)   642.968 sec 

  fun.call <- match.call()

  # param check --------------------------------------------
  method <- match.arg(method, COMPARE_METRICS)
  args <- as.list(environment())
  all_checks <- character(0)
  char_check <- check_fun_params(list(use.type = args$use.type,
                                      method = args$method,
                                      score.strat = args$score.strat,
                                      output.report = args$output.report),
                                 c(1, 1, 1, 1), c(FALSE, FALSE, FALSE, TRUE),
                                 TYPE_CHAR)
  num_check <- check_fun_params(list(compare.to = args$compare.to,
                                     use.freq = args$use.freq,
                                     min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     max.p = args$max.p,
                                     max.e = args$max.e,
                                     nthreads = args$nthreads,
                                     min.position.ic = args$min.position.ic,
                                     output.report.max.print = args$output.report.max.print),
                                c(0, rep(1, 8)), c(TRUE, rep(FALSE, 7), TRUE),
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
  if (!missing(db.scores) && !is.data.frame(db.scores)
      && !is(db.scores, "DataFrame")) {
    dbscores_check <- paste0(" * Incorrect type for 'db.scores: expected ",
                             "`data.frame` or `DataFrame`; got `",
                             class(db.scores), "`")
    all_checks <- c(all_checks, dbscores_check)
  }
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (use.type == "ICM" && method %in% c("ALLR", "ALLR_LL"))
    stop("'use.type = \"ICM\"' is not allowed for ALLR/ALLR_LL methods")

  if (!score.strat %in% c("sum", "a.mean", "g.mean", "median", "wa.mean", "wg.mean"))
    stop("'score.strat' must be one of 'sum', 'a.mean', 'g.mean', 'median', 'wa.mean', 'wg.mean'")

  if (score.strat %in% c("g.mean", "wg.mean") && method %in% c("ALLR", "ALLR_LL", "PCC"))
    stop(wmsg("'g.mean'/'wg.mean' is not allowed for methods which can generate negative ",
              "values: ALLR, ALLR_LL, PCC"))

  if (!missing(compare.to) && method == "IS")
    warning(wmsg("P-values from `method = \"IS\"` will likely be very inaccurate, ",
                 "as random scores from this method are usually very skewed"),
            immediate. = TRUE)

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
  mot.bkgs <- get_bkgs(motifs, use.freq)
  mot.nsites <- get_nsites(motifs)

  if (missing(compare.to)) {

    comparisons <- compare_motifs_all(mot.mats, method, min.overlap,
                                      tryRC, min.mean.ic, normalise.scores,
                                      nthreads, mot.bkgs, mot.type,
                                      relative_entropy, mot.names,
                                      min.position.ic, mot.nsites, score.strat)

    comparisons[comparisons == min_max_doubles()$min
                | comparisons == min_max_doubles()$max] <- NA_real_

    if (anyNA(comparisons)) 
      warning(wmsg("Some comparisons failed due to low motif IC"),
              immediate. = TRUE)

  } else {

    if (missing(db.scores)) {

      db.scores <- JASPAR2018_CORE_DBSCORES

      if (use.freq != 1)
        warning(wmsg("Using the internal P-value database with `use.freq > 1`",
                     " will likely result in incorrect P-values"),
                immediate. = TRUE)
      mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
      if (!all(mot.alphs == "DNA"))
        warning(wmsg("Using the internal P-value database with non-DNA motifs ",
                     "will likely result in incorrect P-values"),
                immediate. = TRUE)

    } 

    db.scores <- check_db_scores(db.scores, method, normalise.scores, score.strat)

    comps <- get_comp_indices(compare.to, length(motifs))
    comparisons <- compare_motifs_cpp(mot.mats, comps[, 1] - 1, comps[, 2] - 1,
                                      method, min.overlap, tryRC, mot.bkgs,
                                      mot.type, relative_entropy,
                                      min.mean.ic, normalise.scores, nthreads,
                                      min.position.ic, mot.nsites, score.strat)

    mot.ncols <- vapply(mot.mats, ncol, integer(1))
    pvals <- pval_extractor(mot.ncols, comparisons, as.integer(comps[[1]] - 1),
                            as.integer(comps[[2]] - 1), method, as.vector(db.scores$subject),
                            as.vector(db.scores$target), as.vector(db.scores$paramA),
                            as.vector(db.scores$paramB), as.vector(db.scores$distribution),
                            nthreads)

    if (any(abs(comparisons) == min_max_doubles()$max))
      warning(wmsg("Some comparisons failed due to low motif IC"),
              immediate. = TRUE)

    comparisons <- DataFrame(subject = mot.names[as.vector(comps[, 1])],
                             subject.i = as.vector(comps[, 1]),
                             target = mot.names[as.vector(comps[, 2])],
                             target.i = as.vector(comps[, 2]),
                             score = comparisons,
                             logPval = pvals,
                             Pval = exp(1)^pvals,
                             Eval = exp(1)^pvals * length(motifs) * 2)

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

    comparisons@metadata <- args[-1]

    if (!missing(output.report)) {
      passed <- tryCatch(compare_motifs_reporter(fun.call, comparisons,
                                                 output.report,
                                                 output.report.max.print, args),
                         error = function(e) FALSE)
      if (isFALSE(passed))
        warning("! Failed to generate output report", immediate. = TRUE)
    }

  }

  comparisons

}

check_db_scores <- function(db.scores, method, normalise.scores, score.strat) {

  if (!is.data.frame(db.scores) && !is(db.scores, "DataFrame")) {
    stop("'db.scores' must be a data.frame or a DataFrame")
  }

  if (!any(method %in% as.vector(db.scores$method))) {
    stop("could not find method '", method, "' in 'db.scores$method'")
  } else {
    db.scores <- db.scores[db.scores$method == method, ]
  }

  if (!any(normalise.scores %in% as.vector(db.scores$normalised))) {
    stop("'db.scores' column 'normalised' does not match 'normalise.scores'")
  } else {
    db.scores <- db.scores[db.scores$normalised == normalise.scores, ]
  }

  if (!any(score.strat %in% as.vector(db.scores$strat))) {
    stop("'db.scores' column 'strat' does not match 'strat'")
  } else {
    db.scores <- db.scores[db.scores$strat == score.strat, ]
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

get_rid_of_dupes <- function(comparisons) {
  comparisons.sorted <- character(nrow(comparisons))
  comparisons.sorted <- vapply(seq_len(nrow(comparisons)),
                               function(x) paste0(sort(c(as.vector(comparisons[x, 2]),
                                                         as.vector(comparisons[x, 4]))),
                                                  collapse = " "),
                               character(1))
  comparisons[!duplicated(comparisons.sorted), ]
} 

compare_motifs_all <- function(mot.mats, method, min.overlap, RC, min.mean.ic,
                               normalise.scores, nthreads, mot.bkgs,
                               mot.type, relative, mot.names,
                               min.position.ic, mot.nsites, score.strat) {

  ans <- compare_motifs_all_cpp(mot.mats, method, min.overlap, RC, mot.bkgs,
                                mot.type, relative, min.mean.ic,
                                normalise.scores, nthreads, min.position.ic,
                                mot.nsites, score.strat)

  n <- length(mot.mats)
  comp <- comb2_cpp(n)
  get_comparison_matrix(unlist(ans), comp[[1]], comp[[2]], method, mot.names)

}

compare_motifs_reporter <- function(fun.call, res, output, max.print, args) {

  motifs <- args$motifs

  if (max.print > nrow(res)) max.print <- nrow(res)

  which.m <- c(as.vector(res[seq_len(max.print), 2]),
               as.vector(res[seq_len(max.print), 4]))
  motifs <- motifs[which.m]

  f <- tempfile(fileext = ".Rmd")
  on.exit(unlink(f))

  m <- tempfile(fileext = ".RDS")
  on.exit(unlink(m))

  saveRDS(motifs, m)

  out <- c("---",
           "title: universalmotif::compare_motifs() results",
           paste0("date: ", Sys.time()),
           "output: html_document",
           "---",
           "",
           "```{r, echo=FALSE}",
           "library(universalmotif)",
           paste0("motifs <- readRDS('", m, "')"),
           "```",
           "",
           "### Function call",
           "",
           paste0("`", deparse(fun.call), "`"),
           "",
           "```{r, eval=FALSE}",
           strsplit(as.yaml(args[-(1:4)]), "\n", fixed = TRUE)[[1]],
           "```",
           "")

  passed <- TRUE

  if (requireNamespace("knitr", quietly = TRUE)) {

    for (i in seq_len(max.print)) {
      if (motifs[[i]]@name == motifs[[i + max.print]]@name)
        motifs[[i + max.print]]@name <- paste(motifs[[i + max.print]]@name, "(duplicated)")
      out <- c(out,
               "---",
               "",
               paste0("## ", i, ": ", motifs[[i]]@name, " [", as.vector(res[i, 2]),
                      "]", " -- ", motifs[[i + max.print]]@name, " [",
                      as.vector(res[i, 4]), "]"),
               "",
               paste0("**Score:** ", as.vector(res[i, 5])),
               "",
               paste0("**logP-value:** ", as.vector(res[i, 6])),
               "",
               paste0("**P-value:** ", as.vector(res[i, 7])),
               "",
               paste0("**E-value:** ", as.vector(res[i, 8])),
               "",
               "```{r, echo=FALSE, out.height=2}",
               paste0("view_motifs(motifs[c(", i, ",", i + max.print, ")],",
                      "method='", args$method, "',", "tryRC=", args$tryRC, ",",
                      "min.overlap=", args$min.overlap, ",", "min.mean.ic=",
                      args$min.mean.ic, ",relative_entropy=", args$relative_entropy,
                      ",normalise.scores=", args$normalise.scores, ",min.position.ic=",
                      args$min.position.ic, ",score.strat='", args$score.strat, "')"),
               "```",
               "")
    }

  } else {
    warning("knitr is required to generate results report", immediate. = TRUE)
    passed <- FALSE
  }

  writeLines(out, f)

  if (passed) {
    if (requireNamespace("rmarkdown", quietly = TRUE)) {
      rmarkdown::render(f, output_file = basename(output),
                        output_dir = dirname(output), quiet = TRUE)
    } else {
      warning("rmarkdown is required to generate results report", immediate. = TRUE)
      passed <- FALSE
    }
  }

  passed

}
