#' Faster minimalist motif comparison.
#'
#' `compare_motifs_lite()` is a deliberately pared-down counterpart to
#' [compare_motifs()]. It uses the
#' per-column Pearson correlation as its similarity metric, computes
#' significance from a TOMTOM-style null PMF (either empirical, from the
#' target database, or parametric, via a Dirichlet-Multinomial over a K=5
#' simplex grid), applies both Bonferroni and Benjamini-Hochberg adjustment,
#' and returns either a long-format `data.frame` of the significant pairs or
#' a square score / p-value / q-value matrix.
#'
#' Use [compare_motifs()] when you need any of the broader features:
#' alternative similarity metrics (EUCL, ALLR, KL, …), score-aggregation
#' strategies, IC-based filters, multifreq scoring, ICM/PPM dispatch,
#' relative-entropy mode, or the JASPAR-derived `db.scores` lookup-table
#' p-value path.
#'
#' P-values are computed in two modes:
#'
#' \describe{
#'   \item{`null = "empirical"` (default)}{Per-query column, score
#'     against every target-database column, histogram into 51 PCC bins,
#'     convolve to overlap length L. Cost per query is O(N x T) where N
#'     is the number of target motifs and T their mean width. Best for
#'     databases with more than a few hundred total columns.}
#'   \item{`null = "parametric"`}{Enumerate the 56 compositions of A/C/G/T
#'     on a K=5 simplex grid, weight each by a Dirichlet-Multinomial PMF
#'     under alpha = 4 * bkg, score against the query column. Per-query
#'     cost is constant in database size. Best for small databases or
#'     when histogram sparsity would distort the empirical estimate.}
#' }
#'
#' The reported `score` is the mean PCC over the overlap (always in
#' \eqn{[-1, 1]}), not the raw sum.
#'
#' @param motifs See [convert_motifs()] for accepted motif formats. DNA or
#'   RNA only.
#' @param compare.to `integer`, `character`, or `NULL`. When `NULL`
#'   (default), all-pairs mode: every motif is compared against every
#'   other and a square matrix is returned (see `matrix.out`). When
#'   `integer`/`character`, the named/indexed motifs are treated as
#'   queries and the *entire* input is treated as the target database;
#'   a long-format `data.frame` of significant matches is returned.
#' @param qvalue `numeric(1)`. q-value cutoff for the long-format output.
#'   Default `0.1`. Ignored in matrix mode.
#' @param min.overlap `integer(1)`. Minimum overlap length (in columns)
#'   for a candidate alignment. Default `5`.
#' @param RC `logical(1)`. If `TRUE` (default), also score query motifs
#'   in reverse-complement orientation and report the best of the two.
#' @param null `character(1)`. Null-PMF mode: `"empirical"` (default) or
#'   `"parametric"`. See Details.
#' @param bkg `numeric(4)` or `NULL`. Background base frequencies (A, C,
#'   G, T). If `NULL`, uniform `c(.25, .25, .25, .25)` is used. Critical
#'   for `"parametric"` (sets the Dirichlet prior); affects only
#'   pseudocount smoothing for `"empirical"`.
#' @param matrix.out `character(1)`. In matrix mode, which metric to
#'   return: `"score"` (mean PCC; default), `"pvalue"`, or `"qvalue"`.
#'   Ignored in long-format mode.
#' @param nthreads `numeric(1)`. Number of threads. `nthreads = 0` uses
#'   all available threads.
#' @param output.report `character(1)`. Path (ending in `.html`) for a
#'   self-contained HTML report of the matches: a summary header, a sortable
#'   results table (query, target, offset, strand, overlap, score, P/q-value),
#'   and a query-vs-target sequence-logo alignment per top match. Requires
#'   `compare.to` to be set and the `knitr` and `rmarkdown` packages (plus
#'   pandoc). Only meaningful in long-format mode.
#' @param output.report.max.print `numeric(1)`. Maximum number of top matches
#'   to show in the report (table and figures). Default `10`.
#'
#' @return
#' - **Matrix mode** (`compare.to = NULL`): a square numeric matrix
#'   indexed by motif name, holding the value selected by `matrix.out`.
#' - **Long-format mode** (`compare.to` provided): a `data.frame` with
#'   columns
#'   `subject`, `subject.i`, `target`, `target.i`, `offset`, `strand`,
#'   `overlap`, `score`, `subject.consensus`, `target.consensus`,
#'   `pvalue`, `qvalue`, filtered to `qvalue <= qvalue` and sorted by
#'   q-value ascending. Coordinates are 1-based; `offset` is the query
#'   column where target column 1 aligns.
#'
#' @references
#'
#' Gupta S, Stamatoyannopoulos JA, Bailey TL, Noble WS (2007). "Quantifying
#' similarity between motifs." *Genome Biology*, **8**(2), R24.
#'
#' @examples
#' library(universalmotif)
#' set.seed(1)
#' motifs <- lapply(1:6, function(i)
#'   create_motif(paste(sample(c("A","C","G","T"), 8, replace = TRUE),
#'                      collapse = ""), name = paste0("M", i)))
#' ## Matrix of mean PCC scores
#' m <- compare_motifs_lite(motifs)
#' round(m, 2)
#' ## Long-format hits at q <= 0.5 with motif 1 as the query
#' compare_motifs_lite(motifs, compare.to = 1, qvalue = 0.5)
#'
#' @seealso [compare_motifs()], [scan_sequences_lite()], [view_motifs_lite()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @family lite motif functions
#' @export
compare_motifs_lite <- function(motifs,
                            compare.to    = NULL,
                            qvalue        = 0.1,
                            min.overlap   = 5L,
                            RC            = TRUE,
                            null          = c("empirical", "parametric"),
                            bkg           = NULL,
                            matrix.out    = c("score", "pvalue", "qvalue"),
                            nthreads      = 1,
                            output.report,
                            output.report.max.print = 10) {

  fun.call <- match.call()

  ## --- arg validation -------------------------------------------------
  if (missing(motifs)) stop("`motifs` is required", call. = FALSE)
  null       <- match.arg(null)
  matrix.out <- match.arg(matrix.out)

  ## Fail fast on output.report misuse, before any comparison work.
  if (!missing(output.report)) {
    .cr_check_output(output.report)
    if (is.null(compare.to))
      stop("`output.report` requires `compare.to` to be set; for all-vs-all ",
           "output see motif_tree().", call. = FALSE)
    for (pkg in c("knitr", "rmarkdown"))
      if (!requireNamespace(pkg, quietly = TRUE))
        stop("`output.report` requires the '", pkg, "' package.",
             call. = FALSE)
  }

  if (!is.numeric(qvalue) || length(qvalue) != 1L || is.na(qvalue) ||
      qvalue <= 0 || qvalue > 1)
    stop("`qvalue` must be a single numeric in (0, 1]", call. = FALSE)
  if (!is.numeric(min.overlap) || length(min.overlap) != 1L ||
      is.na(min.overlap) || min.overlap < 1)
    stop("`min.overlap` must be a positive integer", call. = FALSE)
  min.overlap <- as.integer(min.overlap)
  if (!isTRUEorFALSE(RC))
    stop("`RC` must be a single logical", call. = FALSE)

  ## Background
  if (is.null(bkg)) {
    bkg <- c(A = .25, C = .25, G = .25, T = .25)
  } else {
    if (!is.numeric(bkg) || length(bkg) != 4L || anyNA(bkg) || any(bkg < 0))
      stop("`bkg` must be a length-4 non-negative numeric", call. = FALSE)
    if (abs(sum(bkg) - 1) > 1e-6)
      stop("`bkg` must sum to 1", call. = FALSE)
  }

  nthreads <- resolve_nthreads(nthreads)

  ## --- normalise motifs ------------------------------------------------
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  if (length(motifs) < 2L && is.null(compare.to))
    stop("need at least 2 motifs for all-pairs mode", call. = FALSE)

  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) != 1L)
    stop("all motifs must share the same alphabet", call. = FALSE)
  mot.alph <- unique(mot.alphs)
  if (!mot.alph %in% c("DNA", "RNA"))
    stop("`compare_motifs_lite()` only supports DNA/RNA motifs; got `",
         mot.alph, "`. Use `compare_motifs()` for other alphabets.",
         call. = FALSE)

  ## Coerce to PPM. PPM matrices are used as-is (no pseudocount smoothing).
  ## Callers who want pre-smoothed motifs should apply their own treatment
  ## before passing the motifs in.
  motifs.ppm <- convert_type_internal(motifs, "PPM")
  mot.names <- vapply(motifs.ppm, function(x) x@name, character(1))
  if (any(duplicated(mot.names))) mot.names <- make.unique(mot.names)

  bkg_v <- unname(bkg)
  mot.mats <- lapply(motifs.ppm, function(x) x@motif)

  n_mot <- length(motifs.ppm)
  matrix.mode <- is.null(compare.to)

  ## --- query indices ---------------------------------------------------
  if (matrix.mode) {
    qix <- seq_len(n_mot)
  } else {
    if (is.character(compare.to)) {
      qix <- match(compare.to, mot.names)
      if (anyNA(qix))
        stop("`compare.to` names not found: ",
             paste(compare.to[is.na(qix)], collapse = ", "),
             call. = FALSE)
    } else if (is.numeric(compare.to)) {
      qix <- as.integer(compare.to)
      if (any(qix < 1L | qix > n_mot))
        stop("`compare.to` indices must be in 1:", n_mot, call. = FALSE)
    } else {
      stop("`compare.to` must be NULL, integer, or character",
           call. = FALSE)
    }
  }

  ## --- pair index: every (query, target) -------------------------------
  tix <- seq_len(n_mot)
  pairs <- expand.grid(qi = qix, ti = tix, KEEP.OUT.ATTRS = FALSE)

  ## --- alignment scan --------------------------------------------------
  al <- compare_motifs_lite_align_cpp(mot.mats,
                                  qi          = as.integer(pairs$qi),
                                  ti          = as.integer(pairs$ti),
                                  min_overlap = min.overlap,
                                  RC          = RC,
                                  nthreads    = nthreads)
  ## attach pair indices for downstream
  al$qi <- pairs$qi
  al$ti <- pairs$ti

  ## --- p-value (skip when caller only wants the score matrix) ----------
  ## The best alignment is chosen by score in compare_motifs_lite_align_cpp
  ## (yamcmp.c:1310-1312); p-values are a post-hoc significance step on
  ## that already-chosen alignment. Matrix-mode `matrix.out = "score"`
  ## never reads `pvals` or `qvals`, so we can skip both calls entirely.
  need_pvals <- !matrix.mode || matrix.out != "score"
  if (need_pvals) {
    null_mode <- if (null == "empirical") 0L else 1L
    pvals <- compare_motifs_lite_pvalue_cpp(
      query_mats  = mot.mats,
      target_mats = mot.mats,
      qi          = as.integer(al$qi),
      qstart      = as.integer(al$q_start_oriented),
      overlap     = as.integer(al$overlap),
      strand      = as.integer(al$strand),
      score       = as.numeric(al$score),
      n_tested    = as.integer(al$n_tested),
      bkg         = bkg_v,
      null_mode   = null_mode,
      nthreads    = nthreads
    )
    qvals <- rep(NA_real_, length(pvals))
    for (qq in unique(pairs$qi)) {
      rows <- which(pairs$qi == qq)
      qvals[rows] <- stats::p.adjust(pvals[rows], method = "BH")
    }
  } else {
    pvals <- rep(NA_real_, nrow(al))
    qvals <- rep(NA_real_, nrow(al))
  }

  ## --- mean PCC for display -------------------------------------------
  mean_pcc <- ifelse(al$overlap > 0L, al$score / al$overlap, NA_real_)
  strand_chr <- ifelse(al$strand == 1L, "-", "+")

  ## --- assemble result -------------------------------------------------
  if (matrix.mode) {
    if (qvalue != 0.1)
      warning("`qvalue` is ignored when `compare.to = NULL` (matrix mode); ",
              "filter the returned matrix yourself.", call. = FALSE)
    val <- switch(matrix.out,
                  "score"  = mean_pcc,
                  "pvalue" = pvals,
                  "qvalue" = qvals)
    out <- matrix(val,
                  nrow = n_mot, ncol = n_mot,
                  dimnames = list(mot.names, mot.names))
    return(out)
  }

  ## Long format: per query, keep the best alignment per (query, target)
  ## pair, filter by qvalue, and decorate with consensus strings.
  L <- al$overlap
  qstart_orig <- ifelse(
    al$strand == 1L,
    ## RC: oriented col q_start corresponds to original col w_q - 1 - (q_start + L - 1)
    vapply(seq_along(L), function(k) {
      wq <- ncol(mot.mats[[al$qi[k]]])
      wq - 1L - (al$q_start_oriented[k] + L[k] - 1L)
    }, integer(1)),
    al$q_start_oriented
  )
  ## target start: derived from offset (d_in_q) and overlap
  ## t_start in original target coords = max(0, -offset)
  t_start <- pmax(0L, -al$offset)

  long <- data.frame(
    subject     = mot.names[al$qi],
    subject.i   = al$qi,
    target      = mot.names[al$ti],
    target.i    = al$ti,
    offset      = al$offset,
    strand      = strand_chr,
    overlap     = L,
    score       = mean_pcc,
    pvalue      = pvals,
    qvalue      = qvals,
    stringsAsFactors = FALSE
  )

  ## attach consensus strings (start is 0-based for the C++ helper)
  long$subject.consensus <- compare_motifs_lite_consensus_cpp(
    mot.mats,
    mot_i = as.integer(long$subject.i),
    start = as.integer(qstart_orig),
    len   = as.integer(L)
  )
  long$target.consensus <- compare_motifs_lite_consensus_cpp(
    mot.mats,
    mot_i = as.integer(long$target.i),
    start = as.integer(t_start),
    len   = as.integer(L)
  )

  ## filter + sort
  keep <- !is.na(long$qvalue) & long$qvalue <= qvalue & long$overlap > 0L
  long <- long[keep, , drop = FALSE]
  long <- long[order(long$qvalue, long$pvalue, long$target.i), , drop = FALSE]
  rownames(long) <- NULL

  ## final column order
  long <- long[, c("subject", "subject.i", "target", "target.i",
                   "offset", "strand", "overlap",
                   "score", "subject.consensus", "target.consensus",
                   "pvalue", "qvalue")]

  if (!missing(output.report)) {
    tryCatch(
      compare_report_from_lite(fun.call, long, motifs, mot.names,
                             list(RC = RC, null = null, qvalue = qvalue,
                                  min.overlap = min.overlap),
                             qix, output.report, output.report.max.print),
      error = function(e)
        warning("Failed to generate output report: ", conditionMessage(e),
                immediate. = TRUE, call. = FALSE))
  }

  long
}
