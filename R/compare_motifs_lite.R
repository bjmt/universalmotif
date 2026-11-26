#' Faster minimalist motif comparison.
#'
#' `compare_motifs_lite()` is a deliberately pared-down counterpart to
#' [compare_motifs()]. It uses the
#' per-column Pearson correlation as its similarity metric (or the Euclidean
#' distance, with `method = "EUCL"`), computes
#' significance from a TOMTOM-style null PMF (either empirical, from the
#' target database, or parametric, via a Dirichlet-Multinomial over a K=5
#' simplex grid), applies both Bonferroni and Benjamini-Hochberg adjustment,
#' and returns either a long-format `data.frame` of the significant pairs or
#' a square score / p-value / q-value matrix.
#'
#' Use [compare_motifs()] when you need any of the broader features:
#' other similarity metrics (ALLR, KL, …), score-aggregation
#' strategies, IC-based filters, multifreq scoring, ICM/PPM dispatch,
#' relative-entropy mode, or the JASPAR-derived `db.scores` lookup-table
#' p-value path. To name motifs after the database motifs they resemble most,
#' see [name_motifs()].
#'
#' P-values are computed in two modes:
#'
#' \describe{
#'   \item{`null = "empirical"` (default)}{Per-query column, score
#'     against every target-database column, histogram into 201 bins,
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
#' @details
#' **Euclidean distances.** With `method = "EUCL"`, columns are compared by
#' the Euclidean distance between their probabilities. The `score` is then
#' the mean distance, from 0 for identical columns up to the square root of
#' 2, so that lower means more alike. Summed distances favour short
#' overlaps, which is why these alignments are chosen by their mean distance
#' unless `complete` or background overhangs make every alignment count the
#' same number of columns. P-values are worked out exactly as for
#' correlations, from either null, with distances counting as negative
#' similarities.
#'
#' **Overhanging columns.** With `overhang = "background"`, the columns of
#' either motif hanging over the ends of the other are compared with `bkg`,
#' so that every column of both motifs is scored against whatever faces it.
#' The alignment with the highest total is kept, `score` is the mean of all
#' these column scores, and the P-value comes from the exact distribution of
#' the total, with the target's overhanging columns drawn from the null as
#' well. A partial overlap thus pays for however much its overhanging
#' columns differ from the background. With correlations and a uniform
#' background the overhanging columns score 0, which leaves the alignments
#' much as they are and only lowers the scores of partial overlaps.
#'
#' **Complete scores.** By default only the overlap counts, so that a short
#' target matching part of a long query can score as well as one covering
#' all of it. With `complete = TRUE` whole queries are scored instead, as
#' TOMTOM does by default. Each query column's scores are measured against
#' the median of its scores with the null columns, which makes a typical
#' column score zero and lets query columns hanging over the ends of the
#' target count as typical ones. The P-value of the best alignment is then
#' the probability that the best of all the alignments tried reaches its
#' score, taking them to be independent, rather than a Bonferroni
#' correction. This favours alignments covering more of the query and makes
#' the results depend on which motif is the query, so it suits searching a
#' database with a query motif (as [name_motifs()] does) better than
#' all-against-all comparisons. The reported `score` is still the mean
#' column score over the overlap. Complete scores cannot be combined with
#' background overhangs.
#'
#' **Gapped motifs.** A gapped motif (see [add_gap()]) is compared as each
#' of the motifs it stands for, one for every combination of its gap sizes,
#' with each gap filled by columns of the background (`bkg`). Every
#' combination of one motif is compared with every combination of the
#' other, the comparison with the lowest P-value is kept, and its P-value
#' is multiplied by the number of combinations tried. Its score is then
#' taken over the columns that face no gap. A motif may have at most 32
#' combinations of gap sizes. With `use.gaps = FALSE` the gaps are ignored
#' and the columns compared as if they followed one another.
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
#'   for `"parametric"` (sets the Dirichlet prior); also used by
#'   `overhang = "background"` and to fill the gaps of gapped motifs.
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
#' @param complete `logical(1)`. Score whole queries rather than the overlap
#'   alone, as described in Details. Default `FALSE`.
#' @param method `character(1)`. How columns are compared: `"PCC"`, by their
#'   Pearson correlation (default), or `"EUCL"`, by their Euclidean
#'   distance.
#' @param overhang `character(1)`. `"ignore"` (default) to score only the
#'   overlapping columns, or `"background"` to compare overhanging columns
#'   with `bkg`, as described in Details.
#' @param use.gaps `logical(1)`. Compare gapped motifs with their gaps, as
#'   described in Details. Default `TRUE`.
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
#'   column where target column 1 aligns. When gapped motifs are compared
#'   with their gaps, two more columns, `subject.gaps` and `target.gaps`,
#'   give the size of each gap in the comparison kept (separated by commas,
#'   and empty for motifs without gaps); the offset, overlap and consensus
#'   sequences are then those of the motifs with their gaps filled to these
#'   sizes.
#'
#' Q-values are adjusted among the comparisons of each query with the other
#' motifs. A query's comparison with itself keeps its P-value as its q-value
#' and is left out of the others' adjustment.
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
#' ## The same, scoring the whole query as TOMTOM does
#' compare_motifs_lite(motifs, compare.to = 1, qvalue = 0.5, complete = TRUE)
#' ## Mean Euclidean distances, lower meaning more alike
#' round(compare_motifs_lite(motifs, method = "EUCL"), 2)
#'
#' @seealso [compare_motifs()], [name_motifs()], [scan_sequences_lite()],
#'   [view_motifs_lite()]
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
                            output.report.max.print = 10,
                            complete      = FALSE,
                            method        = c("PCC", "EUCL"),
                            overhang      = c("ignore", "background"),
                            use.gaps      = TRUE) {

  fun.call <- match.call()

  ## --- arg validation -------------------------------------------------
  if (missing(motifs)) stop("`motifs` is required", call. = FALSE)
  null       <- match.arg(null)
  matrix.out <- match.arg(matrix.out)
  method     <- match.arg(method)
  overhang   <- match.arg(overhang)

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
  if (!isTRUEorFALSE(complete))
    stop("`complete` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(use.gaps))
    stop("`use.gaps` must be a single logical", call. = FALSE)
  if (complete && overhang == "background")
    stop("complete scores and background overhangs cannot be combined",
         call. = FALSE)

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

  opts <- list(min.overlap = min.overlap, RC = RC, null = null, bkg = bkg_v,
               method = method, overhang = overhang, complete = complete,
               nthreads = nthreads)

  ## --- alignments and P-values -----------------------------------------
  ## Matrix-mode `matrix.out = "score"` never reads P-values or q-values, so
  ## they are skipped there, except for gapped motifs, whose comparisons are
  ## chosen by their P-values.
  gapped <- use.gaps && lite_any_gapped(motifs.ppm)
  if (gapped) {
    expansion <- lite_expand(motifs.ppm, bkg_v)
    inst.q <- unlist(expansion$of[qix], use.names = FALSE)
    rows <- lite_compare_pairs(expansion$mats, inst.q, opts, need.pvals = TRUE)
    rows <- lite_collapse(rows, expansion, qix, opts)
    all.mats <- expansion$mats
  } else {
    rows <- lite_compare_pairs(mot.mats, qix, opts,
                               need.pvals = !matrix.mode || matrix.out != "score")
    all.mats <- mot.mats
  }

  ## --- assemble result -------------------------------------------------
  if (matrix.mode) {
    if (qvalue != 0.1)
      warning("`qvalue` is ignored when `compare.to = NULL` (matrix mode); ",
              "filter the returned matrix yourself.", call. = FALSE)
    val <- switch(matrix.out,
                  "score"  = rows$score,
                  "pvalue" = rows$pvalue,
                  "qvalue" = rows$qvalue)
    out <- matrix(NA_real_, nrow = n_mot, ncol = n_mot,
                  dimnames = list(mot.names, mot.names))
    out[cbind(rows$qi, rows$ti)] <- val
    return(out)
  }

  ## Long format: the best alignment per (query, target) pair, filtered by
  ## qvalue and decorated with consensus strings.
  L <- rows$overlap
  wq <- vapply(all.mats, ncol, integer(1))[rows$iq]
  ## RC: oriented col q_start corresponds to original col w_q - 1 - (q_start + L - 1)
  qstart_orig <- ifelse(rows$strand == 1L,
                        wq - 1L - (rows$q_start_oriented + L - 1L),
                        rows$q_start_oriented)
  ## target start in original target coords = max(0, -offset)
  t_start <- pmax(0L, -rows$offset)

  long <- data.frame(
    subject     = mot.names[rows$qi],
    subject.i   = rows$qi,
    target      = mot.names[rows$ti],
    target.i    = rows$ti,
    offset      = rows$offset,
    strand      = ifelse(rows$strand == 1L, "-", "+"),
    overlap     = L,
    score       = rows$score,
    pvalue      = rows$pvalue,
    qvalue      = rows$qvalue,
    stringsAsFactors = FALSE
  )

  ## attach consensus strings (start is 0-based for the C++ helper)
  long$subject.consensus <- compare_motifs_lite_consensus_cpp(
    all.mats,
    mot_i = as.integer(rows$iq),
    start = as.integer(qstart_orig),
    len   = as.integer(L),
    alph  = mot.alph
  )
  long$target.consensus <- compare_motifs_lite_consensus_cpp(
    all.mats,
    mot_i = as.integer(rows$it),
    start = as.integer(t_start),
    len   = as.integer(L),
    alph  = mot.alph
  )
  cols <- c("subject", "subject.i", "target", "target.i",
            "offset", "strand", "overlap",
            "score", "subject.consensus", "target.consensus",
            "pvalue", "qvalue")
  if (gapped) {
    sizes <- vapply(expansion$sizes, paste, character(1), collapse = ",")
    long$subject.gaps <- sizes[rows$iq]
    long$target.gaps <- sizes[rows$it]
    cols <- c(cols, "subject.gaps", "target.gaps")
  }

  ## filter + sort
  keep <- !is.na(long$qvalue) & long$qvalue <= qvalue & long$overlap > 0L
  long <- long[keep, , drop = FALSE]
  long <- long[order(long$qvalue, long$pvalue, long$target.i), , drop = FALSE]
  rownames(long) <- NULL

  ## final column order
  long <- long[, cols]

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

## ---------------------------------------------------------------------------
## Internal helpers shared by the lite comparison, merging and display
## functions.
## ---------------------------------------------------------------------------

## The default comparison options, as compare_motifs_lite() fills them in.
lite_opts <- function(min.overlap = 5L, RC = TRUE, null = "empirical",
                      bkg = rep(0.25, 4), method = "PCC", overhang = "ignore",
                      complete = FALSE, nthreads = 1L) {
  list(min.overlap = as.integer(min.overlap), RC = RC, null = null, bkg = bkg,
       method = method, overhang = overhang, complete = complete,
       nthreads = as.integer(nthreads))
}

## The median score of each column of the wanted motifs against the null
## columns, for complete scores; empty for the other motifs.
lite_medians <- function(mats, wanted, opts) {
  compare_motifs_lite_medians_cpp(mats, wanted = wanted, RC = opts$RC,
    null_mode = if (opts$null == "empirical") 0L else 1L, bkg = opts$bkg,
    method = if (opts$method == "EUCL") 1L else 0L, nthreads = opts$nthreads)
}

## The best alignment of each (qi, ti) pair of PPM matrices.
lite_align <- function(mats, qi, ti, opts, medians = NULL) {
  if (opts$complete && is.null(medians)) {
    wanted <- logical(length(mats))
    wanted[unique(qi)] <- TRUE
    medians <- lite_medians(mats, wanted, opts)
  }
  if (is.null(medians)) medians <- rep(list(numeric(0)), length(mats))
  compare_motifs_lite_align_cpp(mats,
    qi          = as.integer(qi),
    ti          = as.integer(ti),
    min_overlap = as.integer(opts$min.overlap),
    RC          = as.logical(opts$RC),
    nthreads    = as.integer(opts$nthreads),
    method      = if (opts$method == "EUCL") 1L else 0L,
    overhang    = if (opts$overhang == "background") 1L else 0L,
    complete    = opts$complete,
    medians     = medians,
    bkg         = as.numeric(opts$bkg))
}

## Compare each query in `qix` with every matrix in `mats`, itself included.
## One row per pair, with the query varying fastest: the matrix indices
## (qi, ti, and the same as iq, it), the alignment, the reported score and,
## with `need.pvals`, P-values and q-values adjusted within each query,
## leaving out its comparison with itself.
lite_compare_pairs <- function(mats, qix, opts, need.pvals = TRUE) {
  pairs <- expand.grid(qi = qix, ti = seq_along(mats), KEEP.OUT.ATTRS = FALSE)
  medians <- NULL
  if (opts$complete) {
    wanted <- logical(length(mats))
    wanted[unique(qix)] <- TRUE
    medians <- lite_medians(mats, wanted, opts)
  }
  al <- lite_align(mats, pairs$qi, pairs$ti, opts, medians)
  al$qi <- al$iq <- pairs$qi
  al$ti <- al$it <- pairs$ti

  if (need.pvals) {
    if (is.null(medians)) medians <- rep(list(numeric(0)), length(mats))
    pvals <- compare_motifs_lite_pvalue_cpp(
      query_mats  = mats,
      target_mats = mats,
      qi          = as.integer(al$qi),
      ti          = as.integer(al$ti),
      qstart      = as.integer(al$q_start_oriented),
      overlap     = as.integer(al$overlap),
      strand      = as.integer(al$strand),
      offset      = as.integer(al$offset),
      score       = as.numeric(al$score),
      chosen      = as.numeric(al$chosen),
      n_tested    = as.integer(al$n_tested),
      bkg         = as.numeric(opts$bkg),
      null_mode   = if (opts$null == "empirical") 0L else 1L,
      nthreads    = as.integer(opts$nthreads),
      method      = if (opts$method == "EUCL") 1L else 0L,
      overhang    = if (opts$overhang == "background") 1L else 0L,
      complete    = opts$complete,
      medians     = medians,
      RC          = as.logical(opts$RC),
      min_overlap = as.integer(opts$min.overlap)
    )
    ## A query's comparison with itself is no test, and its tiny P-value
    ## would shrink every other q-value, so it keeps its P-value and is left
    ## out of the others' adjustment.
    qvals <- pvals
    for (qq in unique(al$qi)) {
      rows <- which(al$qi == qq & al$ti != qq)
      qvals[rows] <- stats::p.adjust(pvals[rows], method = "BH")
    }
  } else {
    pvals <- rep(NA_real_, nrow(al))
    qvals <- rep(NA_real_, nrow(al))
  }
  al$pvalue <- pvals
  al$qvalue <- qvals

  ## The mean column score, over the overlap or with background overhangs
  ## over every column of both motifs, with distances as they are.
  widths <- vapply(mats, ncol, integer(1))
  mean.score <- if (opts$overhang == "background") {
    al$chosen / (widths[al$qi] + widths[al$ti])
  } else al$score / al$overlap
  if (opts$method == "EUCL") mean.score <- 0 - mean.score
  al$score.sum <- al$score
  al$score <- ifelse(al$overlap > 0L, mean.score, NA_real_)
  al
}

## Whether any motif has gaps.
lite_any_gapped <- function(motifs) {
  any(vapply(motifs, function(m)
    isTRUE(m@gapinfo@isgapped) && length(m@gapinfo@gaploc) > 0, logical(1)))
}

## Combinations of gap sizes beyond which a gapped motif is refused. Every
## combination of one motif is compared with every one of another, so the
## cost of a pair grows with the product.
LITE_MAX_GAP_COMBOS <- 32

## Every combination of the gap sizes of a gapped motif, smallest first, as
## the rows of a matrix with one column per gap.
lite_gap_combos <- function(motif) {
  g <- motif@gapinfo
  lo <- as.integer(g@mingap)
  hi <- as.integer(g@maxgap)
  n <- prod(hi - lo + 1)
  if (n > LITE_MAX_GAP_COMBOS)
    stop(wmsg(motif@name, " has ", n, " combinations of gap sizes, more than ",
              "the ", LITE_MAX_GAP_COMBOS, " that can be compared; narrow its ",
              "gap range, or set `use.gaps = FALSE` to compare it without its ",
              "gaps"), call. = FALSE)
  ## expand.grid() varies its first column fastest, so the gaps are given
  ## in reverse to vary the last gap fastest.
  ranges <- rev(Map(seq.int, lo, hi))
  out <- as.matrix(expand.grid(ranges, KEEP.OUT.ATTRS = FALSE))
  out <- out[, rev(seq_len(ncol(out))), drop = FALSE]
  dimnames(out) <- NULL
  out
}

## The matrix of a gapped motif with gaps of the given sizes filled by blank
## columns (`blank.col`, such as the background), and which of its columns are
## blank.
lite_gap_instance <- function(mat, gaploc, sizes, blank.col) {
  bounds <- c(0, gaploc, ncol(mat))
  pieces <- list()
  blank <- logical(0)
  for (k in seq_len(length(bounds) - 1L)) {
    if (k > 1L && sizes[k - 1L] > 0) {
      pieces <- c(pieces, list(matrix(blank.col, nrow = nrow(mat),
                                      ncol = sizes[k - 1L])))
      blank <- c(blank, rep(TRUE, sizes[k - 1L]))
    }
    cols <- seq_len(bounds[k + 1L] - bounds[k]) + bounds[k]
    pieces <- c(pieces, list(mat[, cols, drop = FALSE]))
    blank <- c(blank, rep(FALSE, length(cols)))
  }
  out <- do.call(cbind, pieces)
  dimnames(out) <- list(rownames(mat), NULL)
  list(mat = out, blank = blank)
}

## Each motif as the matrices of its combinations of gap sizes, or as its own
## matrix when it has no gaps: the matrices, the motif each stands for, its gap
## sizes and blank columns, and the matrices of each motif.
lite_expand <- function(motifs, bkg) {
  mats <- list()
  origin <- integer(0)
  sizes <- list()
  blank <- list()
  of <- vector("list", length(motifs))
  for (i in seq_along(motifs)) {
    m <- motifs[[i]]
    mat <- m@motif
    if (isTRUE(m@gapinfo@isgapped) && length(m@gapinfo@gaploc) > 0) {
      combos <- lite_gap_combos(m)
      for (r in seq_len(nrow(combos))) {
        inst <- lite_gap_instance(mat, m@gapinfo@gaploc, combos[r, ],
                                  bkg / sum(bkg))
        mats <- c(mats, list(inst$mat))
        blank <- c(blank, list(inst$blank))
        sizes <- c(sizes, list(combos[r, ]))
        origin <- c(origin, i)
        of[[i]] <- c(of[[i]], length(mats))
      }
    } else {
      mats <- c(mats, list(mat))
      blank <- c(blank, list(logical(ncol(mat))))
      sizes <- c(sizes, list(integer(0)))
      origin <- c(origin, i)
      of[[i]] <- length(mats)
    }
  }
  list(mats = mats, origin = origin, sizes = sizes, blank = blank, of = of)
}

## The mean column score of a comparison of two gap instances over the
## columns that face each other, leaving out those facing a blank column in
## either, which would otherwise dilute it. Distances are given as they are.
lite_gap_free_score <- function(qm, tm, q.blank, t.blank, minus, offset,
                                overlap, method) {
  if (minus) {
    qm <- rev_comp_mat(qm)
    q.blank <- rev(q.blank)
  }
  k <- seq_len(overlap) - 1L
  a <- max(offset, 0L) + k + 1L
  b <- max(-offset, 0L) + k + 1L
  ok <- a <= ncol(qm) & b <= ncol(tm)
  a <- a[ok]
  b <- b[ok]
  ok <- !q.blank[a] & !t.blank[b]
  a <- a[ok]
  b <- b[ok]
  if (!length(a)) return(NaN)
  qa <- qm[, a, drop = FALSE]
  tb <- tm[, b, drop = FALSE]
  if (method == "EUCL") return(mean(sqrt(colSums((qa - tb)^2))))
  qd <- qa - 0.25
  td <- tb - 0.25
  sqq <- colSums(qd * qd)
  stt <- colSums(td * td)
  sqt <- colSums(qd * td)
  r <- ifelse(sqq < 1e-12 | stt < 1e-12, 0, sqt / sqrt(sqq * stt))
  mean(r)
}

## Comparisons of gap instances reduced to those of the motifs they stand
## for. Each pair of motifs keeps the comparison of its instances with the
## lowest P-value, or the best score among equal ones, and its P-value is
## multiplied by the number of pairs of instances tried. Comparisons between
## gapped motifs are scored over the columns that face no gap (with overhangs
## ignored), and q-values are worked out again within each query, leaving out
## its comparison with itself.
lite_collapse <- function(rows, expansion, qix, opts) {
  n <- length(expansion$of)
  rows$qi <- expansion$origin[rows$iq]
  rows$ti <- expansion$origin[rows$it]
  ## The instances are compared query instance by query instance, and the
  ## first of equally good comparisons is kept.
  seen <- order(match(rows$iq, unique(rows$iq)), rows$it)
  rows <- rows[seen, , drop = FALSE]
  score.key <- if (opts$method == "EUCL") rows$score else -rows$score
  best <- order(rows$qi, rows$ti, is.nan(rows$score) | is.na(rows$score),
                rows$pvalue, score.key, seq_len(nrow(rows)))
  rows <- rows[best, , drop = FALSE]
  rows <- rows[!duplicated(rows[, c("qi", "ti")]), , drop = FALSE]
  key <- paste(rows$qi, rows$ti)

  out <- rows[match(paste(rep(qix, each = n), rep(seq_len(n), length(qix))),
                    key), , drop = FALSE]
  if (opts$overhang == "ignore") {
    for (k in which(lengths(expansion$sizes[out$iq]) > 0 |
                    lengths(expansion$sizes[out$it]) > 0)) {
      if (out$overlap[k] == 0L) next
      out$score[k] <- lite_gap_free_score(
        expansion$mats[[out$iq[k]]], expansion$mats[[out$it[k]]],
        expansion$blank[[out$iq[k]]], expansion$blank[[out$it[k]]],
        out$strand[k] == 1L, out$offset[k], out$overlap[k], opts$method)
    }
  }
  n.inst <- lengths(expansion$of)
  out$pvalue <- pmin(out$pvalue * n.inst[out$qi] * n.inst[out$ti], 1)
  out$qvalue <- out$pvalue
  for (j in seq_along(qix)) {
    r <- (j - 1L) * n + seq_len(n)
    others <- r[out$ti[r] != qix[j]]
    out$qvalue[others] <- stats::p.adjust(out$pvalue[others], method = "BH")
  }
  rownames(out) <- NULL
  out
}
