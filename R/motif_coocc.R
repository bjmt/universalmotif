#' Find significantly co-occurring motif pairs in a set of sequences.
#'
#' [motif_coocc()] tests every motif pair for over-co-occurrence across
#' a set of sequences using a one-sided Fisher's exact test on the 2x2
#' contingency table of per-sequence presence/absence. The result is
#' the long-format pair table, BH-corrected across all tested pairs.
#' Optionally, when `max.distance` is non-`NULL`, the function also
#' reports two descriptive spatial columns, `both.clustered`
#' (number of co-occurring sequences with a within-`max.distance`
#' (A, B) hit pair) and `median.distance` (the median nearest-pair
#' spacing), so users can flag heterodimer-like arrangements. The
#' Fisher p-value itself is always computed on the unfiltered 2x2
#' (the spatial filter is informational, not a test).
#'
#' Two input paths are supported:
#'
#'   1. Internal scan (default): pass `motifs` + `sequences` and
#'      `motif_coocc()` will call [scan_sequences2()] with the given
#'      `pvalue` and `RC` arguments to build the hit table. DNA/RNA
#'      only (the restriction comes from `scan_sequences2()`).
#'   2. Precomputed hits: pass `motifs` + `hits` + `n.sequences`.
#'      `hits` is a `data.frame` or `GRanges` with `motif.i` and
#'      `sequence.i` columns (the native return shape of
#'      [scan_sequences()] / [scan_sequences2()]), optionally `start`
#'      when `max.distance` is set. This path accepts any motif
#'      alphabet; once a hit table exists, co-occurrence is pure
#'      set arithmetic on integer indices and no sequence bytes are
#'      read.
#'
#' @param motifs See [convert_motifs()] for accepted motif formats.
#'   A single motif or a list of motifs. Required in both paths to
#'   provide motif names and to dimension the pair table.
#' @param sequences `XStringSet` or `NULL`. If supplied, the function
#'   scans internally via [scan_sequences2()]; DNA/RNA only.
#' @param hits `data.frame`, `GRanges`, or `NULL`. Precomputed hit
#'   table with `motif.i` and `sequence.i` columns. If supplied,
#'   `sequences` is ignored and `n.sequences` must be provided. Any
#'   alphabet accepted on this path.
#' @param n.sequences `integer(1)` or `NULL`. Total number of host
#'   sequences (needed to fill the "neither" cell of the 2x2 when
#'   using the precomputed-hits path).
#' @param pvalue `numeric(1)`. Per-hit p-value threshold for the
#'   internal scan. Default `1e-4`. Ignored on the hit-table path.
#' @param RC `logical(1)`. Scan the reverse complement too? Default
#'   `TRUE`. Ignored on the hit-table path.
#' @param max.distance `integer(1)` or `NULL`. If non-`NULL`, two
#'   extra descriptive columns are added to the output:
#'   `both.clustered` (subset of co-occurring sequences with a within-
#'   `max.distance` (A, B) hit pair) and `median.distance` (median of
#'   those nearest-pair distances). The Fisher p-value is unchanged;
#'   it always tests "do A and B occur in the same sequences more
#'   often than chance?" on the unfiltered 2x2. Requires a `start`
#'   column on the hit table.
#' @param min.coocc `integer(1)`. Pairs with fewer than `min.coocc`
#'   co-occurring sequences are tagged with NA p- and q-values
#'   (not tested). Default `1L` (only test pairs that co-occur at
#'   least once).
#' @param pseudocount `numeric(1)`. Added to every cell of the 2x2
#'   contingency before the Fisher test (matches the convention from
#'   [enrich_motifs2()]). Default `0L`.
#' @param self.pairs `logical(1)`. If `TRUE`, include `(i, i)` rows
#'   (testing whether motif `i` co-occurs with itself). Usually only
#'   meaningful for motifs that bind as homodimers. Default `FALSE`.
#' @param nthreads `integer(1)`. Threads for the internal scan (when
#'   `sequences` is supplied). The pairwise Fisher loop is single-
#'   threaded R code regardless.
#'
#' @return `data.frame` with columns:
#'   * `motif_a`, `motif_b`: motif names (1-based indexing into the
#'     supplied `motifs` list; deduplicated names if duplicates exist).
#'   * `a_only`, `b_only`, `both`, `neither`: 2x2 contingency cells.
#'   * `odds_ratio`: conditional MLE odds ratio from `fisher.test`.
#'   * `pvalue`: one-sided Fisher's exact p-value
#'     (alternative = "greater"). NA when `both < min.coocc`.
#'   * `qvalue`: BH-corrected q-value across all tested pairs.
#'   * `both.clustered` (spatial mode only): subset of `both`
#'     where a within-`max.distance` (A, B) hit pair exists.
#'   * `median.distance` (spatial mode only): median of those
#'     nearest-pair distances; `NA` if no clustered pair.
#'   Rows are sorted by q-value ascending; NA-pvalue rows go last.
#'
#' @examples
#' \dontrun{
#' library(universalmotif)
#' data(ArabidopsisMotif)
#' seqs <- create_sequences("DNA", seqnum = 200, seqlen = 500)
#' ## Toy: implant two copies of the same motif jointly in all seqs.
#' set.seed(1)
#' imp <- implant_motifs(list(ArabidopsisMotif, ArabidopsisMotif),
#'                       seqs, n.per.seq = 1, return.indices = TRUE)
#' co <- motif_coocc(list(ArabidopsisMotif, ArabidopsisMotif),
#'                   imp$sequences, pvalue = 1e-3)
#' co
#' }
#'
#' @seealso [scan_sequences2()], [enrich_motifs2()], [motif_peaks()]
#' @references
#'
#' Heinz S, Benner C, Spann N, et al. (2010). "Simple combinations of
#' lineage-determining transcription factors prime cis-regulatory
#' elements required for macrophage and B cell identities."
#' *Molecular Cell*, **38**(4):576-589. \doi{10.1016/j.molcel.2010.05.004}.
#' (HOMER's pairwise-Fisher co-occurrence convention.)
#'
#' Whitington T, Frith MC, Johnson J, Bailey TL (2011). "Inferring
#' transcription factor complexes from ChIP-seq data." *Nucleic Acids
#' Research*, **39**(15):e98. \doi{10.1093/nar/gkr341}. (SpaMo, the
#' spatial-distance motif-pair analysis.)
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
motif_coocc <- function(motifs, sequences = NULL,
                        hits         = NULL,
                        n.sequences  = NULL,
                        pvalue       = 1e-4,
                        RC           = TRUE,
                        max.distance = NULL,
                        min.coocc    = 1L,
                        pseudocount  = 0L,
                        self.pairs   = FALSE,
                        nthreads     = 1L) {

  ## --- arg validation --------------------------------------------------
  if (missing(motifs))
    stop("`motifs` is required", call. = FALSE)
  if (is.null(sequences) && is.null(hits))
    stop("provide either `sequences` (to scan internally) or `hits` ",
         "(with `n.sequences`)", call. = FALSE)
  if (!is.null(sequences) && !is.null(hits))
    stop("provide exactly one of `sequences` or `hits`, not both",
         call. = FALSE)
  if (!is.null(max.distance) &&
      (!is.numeric(max.distance) || length(max.distance) != 1L ||
       max.distance < 0L))
    stop("`max.distance` must be a non-negative integer or NULL",
         call. = FALSE)
  if (!is.numeric(min.coocc) || length(min.coocc) != 1L || min.coocc < 0L)
    stop("`min.coocc` must be a non-negative integer", call. = FALSE)
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      pseudocount < 0)
    stop("`pseudocount` must be a non-negative numeric", call. = FALSE)
  if (!isTRUEorFALSE(self.pairs))
    stop("`self.pairs` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(RC))
    stop("`RC` must be a single logical", call. = FALSE)

  ## --- normalise motifs ------------------------------------------------
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  n.motifs <- length(motifs)
  if (n.motifs < 2L && !self.pairs)
    stop("need at least 2 motifs (or `self.pairs = TRUE`)", call. = FALSE)
  mot.names <- vapply(motifs, function(x) x@name, character(1))
  empty <- !nzchar(mot.names)
  if (any(empty)) mot.names[empty] <- sprintf("motif_%d", which(empty))
  if (any(duplicated(mot.names))) mot.names <- make.unique(mot.names)

  ## --- choose input path -----------------------------------------------
  if (is.null(hits)) {
    ## Internal-scan path (DNA/RNA only -- enforced by scan_sequences2).
    hits <- scan_sequences2(motifs, sequences,
                            pvalue = pvalue, RC = RC,
                            nthreads = nthreads, return.granges = FALSE)
    n.sequences <- length(sequences)
    if (inherits(hits, "GRanges"))
      hits <- as.data.frame(hits)  # safety; should be data.frame here
  } else {
    ## Hit-table path -- any alphabet OK.
    if (is.null(n.sequences) ||
        !is.numeric(n.sequences) || length(n.sequences) != 1L ||
        n.sequences < 1L)
      stop("`n.sequences` must be a positive integer when `hits` is ",
           "supplied", call. = FALSE)
    n.sequences <- as.integer(n.sequences)
    hits <- coerce_hits(hits, motifs, mot.names, need.start = !is.null(max.distance))
  }

  if (!is.null(max.distance) && !"start" %in% names(hits))
    stop("spatial mode (`max.distance` non-NULL) requires a `start` ",
         "column in the hit table", call. = FALSE)

  ## --- per-motif sequence sets -----------------------------------------
  presence <- vector("list", n.motifs)
  for (i in seq_len(n.motifs))
    presence[[i]] <- unique(hits$sequence.i[hits$motif.i == i])

  ## --- for spatial mode: cache hit starts per (motif, sequence) --------
  hit_pos <- NULL
  if (!is.null(max.distance)) {
    hit_pos <- vector("list", n.motifs)
    for (i in seq_len(n.motifs)) {
      sub <- hits[hits$motif.i == i, , drop = FALSE]
      if (nrow(sub) > 0L)
        hit_pos[[i]] <- split(as.integer(sub$start), sub$sequence.i)
      else
        hit_pos[[i]] <- list()
    }
  }

  ## --- build pair index ------------------------------------------------
  pair_idx <- if (self.pairs) {
    do.call(rbind, lapply(seq_len(n.motifs),
                          function(i) cbind(i, i:n.motifs)))
  } else {
    do.call(rbind, lapply(seq_len(n.motifs - 1L),
                          function(i) cbind(i, (i + 1L):n.motifs)))
  }
  storage.mode(pair_idx) <- "integer"
  rows <- nrow(pair_idx)

  ## --- output skeleton -------------------------------------------------
  out <- data.frame(
    motif_a    = character(rows),
    motif_b    = character(rows),
    a_only     = integer(rows),
    b_only     = integer(rows),
    both       = integer(rows),
    neither    = integer(rows),
    odds_ratio = numeric(rows),
    pvalue     = numeric(rows),
    qvalue     = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!is.null(max.distance)) {
    out$both.clustered  <- integer(rows)
    out$median.distance <- NA_real_
  }

  ## --- pairwise loop ---------------------------------------------------
  ## The Fisher test always answers the SET question: do A and B co-occur
  ## in the same sequences more often than chance? The 2x2 partition
  ## (a_only, b_only, both, neither) is therefore computed without the
  ## spatial filter -- otherwise the partition stops summing to
  ## n.sequences and the test result becomes meaningless.
  ##
  ## When `max.distance` is supplied, we additionally report:
  ##   - both.clustered:  subset of `both` where at least one
  ##                      (A, B) hit pair lies within `max.distance` bp.
  ##   - median.distance: median of those nearest-pair distances.
  ## These are descriptive of HOW co-occurring motifs are arranged;
  ## a proper spatial significance test (e.g. SpaMo's permutation of
  ## distance distributions) is out of scope here.
  pc <- as.numeric(pseudocount)
  for (k in seq_len(rows)) {
    i <- pair_idx[k, 1L]
    j <- pair_idx[k, 2L]
    Si <- presence[[i]]
    Sj <- presence[[j]]
    coocc_seqs <- if (i == j) Si else intersect(Si, Sj)

    both    <- length(coocc_seqs)
    a_only  <- length(Si) - both
    b_only  <- length(Sj) - both
    if (i == j) {
      ## self-pair: degenerate 2x2 (Si == Sj so a_only == b_only == 0).
      a_only <- 0L; b_only <- 0L
    }
    neither <- n.sequences - a_only - b_only - both
    if (neither < 0L) neither <- 0L  # numerical guard

    out$motif_a[k] <- mot.names[i]
    out$motif_b[k] <- mot.names[j]
    out$a_only[k]  <- as.integer(a_only)
    out$b_only[k]  <- as.integer(b_only)
    out$both[k]    <- as.integer(both)
    out$neither[k] <- as.integer(neither)

    if (both >= as.integer(min.coocc)) {
      tbl <- matrix(c(both, b_only, a_only, neither), nrow = 2L) + pc
      ft  <- suppressWarnings(stats::fisher.test(tbl, alternative = "greater"))
      out$pvalue[k]     <- ft$p.value
      out$odds_ratio[k] <- unname(ft$estimate)
    } else {
      out$pvalue[k]     <- NA_real_
      out$odds_ratio[k] <- NA_real_
    }

    ## Descriptive spatial columns (do NOT feed the Fisher test).
    if (!is.null(max.distance) && both > 0L) {
      kept_dists <- numeric(0)
      for (seqkey in as.character(coocc_seqs)) {
        sa <- hit_pos[[i]][[seqkey]]
        sb <- hit_pos[[j]][[seqkey]]
        if (length(sa) == 0L || length(sb) == 0L) next
        if (i == j && length(sa) < 2L) next  # self-pair: need >= 2 hits
        ds <- if (i == j)
                outer(sa, sa, function(x, y) abs(x - y))[upper.tri(diag(length(sa)))]
              else
                as.numeric(outer(sa, sb, function(x, y) abs(x - y)))
        d_min <- min(ds)
        if (d_min <= max.distance) kept_dists <- c(kept_dists, d_min)
      }
      out$both.clustered[k] <- length(kept_dists)
      if (length(kept_dists) > 0L)
        out$median.distance[k] <- stats::median(kept_dists)
    }
  }

  ## --- BH q-values across tested pairs only ----------------------------
  tested <- !is.na(out$pvalue)
  if (any(tested))
    out$qvalue[tested] <- stats::p.adjust(out$pvalue[tested], method = "BH")

  ## --- sort: tested rows by q then p, NA rows last ---------------------
  ord <- order(is.na(out$qvalue), out$qvalue, out$pvalue,
               method = "radix")
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}

## ---- internal helpers ---------------------------------------------------

## Coerce a hit table (data.frame or GRanges) into a data.frame with
## integer `motif.i`, `sequence.i`, and (if needed) `start`.
coerce_hits <- function(hits, motifs, mot.names, need.start) {
  if (inherits(hits, "GRanges")) {
    mc <- as.data.frame(S4Vectors::mcols(hits))
    if (!"start" %in% names(mc)) mc$start <- BiocGenerics::start(hits)
    hits <- mc
  }
  if (!is.data.frame(hits))
    stop("`hits` must be a data.frame or GRanges", call. = FALSE)
  if (!"sequence.i" %in% names(hits))
    stop("`hits` must contain a `sequence.i` column", call. = FALSE)

  ## Allow motif.i as integer index OR character motif name.
  if (!"motif.i" %in% names(hits)) {
    if ("motif" %in% names(hits))
      hits$motif.i <- match(as.character(hits$motif), mot.names)
    else
      stop("`hits` must contain a `motif.i` or `motif` column",
           call. = FALSE)
  } else if (is.character(hits$motif.i)) {
    hits$motif.i <- match(hits$motif.i, mot.names)
  }
  if (any(is.na(hits$motif.i)))
    stop("some entries in `hits$motif.i` could not be matched to ",
         "`motifs` -- check names and indices", call. = FALSE)
  hits$motif.i    <- as.integer(hits$motif.i)
  hits$sequence.i <- as.integer(hits$sequence.i)
  if (max(hits$motif.i) > length(motifs))
    stop("`hits` references motif indices beyond length(motifs)",
         call. = FALSE)
  if (need.start && !"start" %in% names(hits))
    stop("`hits` must contain a `start` column for spatial mode",
         call. = FALSE)
  hits
}
