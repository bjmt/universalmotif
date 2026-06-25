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
#' spacing), so heterodimer-like arrangements can be flagged. The
#' Fisher p-value itself is always computed on the unfiltered 2x2
#' (the spatial filter is informational, not a test).
#'
#' Two input paths are supported:
#'
#'   1. Internal scan (default): pass `motifs` + `sequences` and
#'      `motif_coocc()` will call [scan_sequences_lite()] with the given
#'      `pvalue` and `RC` arguments to build the hit table. DNA/RNA
#'      only (the restriction comes from `scan_sequences_lite()`).
#'   2. Precomputed hits: pass `motifs` + `hits` + `n.sequences`.
#'      `hits` is a `data.frame` or `GRanges` with `motif.i` and
#'      `sequence.i` columns (the native return shape of
#'      [scan_sequences()] / [scan_sequences_lite()]), optionally `start`
#'      when `max.distance` is set. This path accepts any motif
#'      alphabet; once a hit table exists, co-occurrence is pure
#'      set arithmetic on integer indices and no sequence bytes are
#'      read.
#'
#' @param motifs See [convert_motifs()] for accepted motif formats.
#'   A single motif or a list of motifs. Required in both paths to
#'   provide motif names and to dimension the pair table.
#' @param sequences `XStringSet` or `NULL`. If supplied, the function
#'   scans internally via [scan_sequences_lite()]; DNA/RNA only.
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
#'   `TRUE`. Ignored when a hits table is provided.
#' @param no.overlaps `logical(1)`. Passed to [scan_sequences_lite()] on the
#'   internal-scan path: if `TRUE`, overlapping hits of the *same* motif
#'   are collapsed to one (greedy; see [dedup_hits()]) before
#'   co-occurrence is computed. This does not change the presence-based
#'   Fisher test (a sequence either contains a motif or not, however many
#'   overlapping copies it has), so it only affects the descriptive
#'   spatial columns reported under `max.distance`, where it stops
#'   tandem/overlapping matches of one motif being counted as separate
#'   instances. Overlaps *between* the two motifs of a pair are always
#'   kept, so heterodimer footprints are preserved. Ignored (with a
#'   warning) when a precomputed `hits` table is supplied. Default
#'   `FALSE`.
#' @param no.overlaps.by `character(1)`. Tie-break priority when
#'   `no.overlaps = TRUE`: `"pvalue"` (default) keeps the lowest-p hit,
#'   `"score"` the highest-scoring. Passed to [scan_sequences_lite()].
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
#'   [enrich_motifs_lite()]). Default `0L`.
#' @param self.pairs `logical(1)`. If `TRUE`, include `(i, i)` rows
#'   (testing whether motif `i` co-occurs with itself). Usually only
#'   meaningful for motifs that bind as homodimers. Default `FALSE`.
#' @param nthreads `integer(1)`. Threads for the internal scan (when
#'   `sequences` is supplied). The pairwise co-occurrence test itself is
#'   vectorised R (one `phyper()` call over all pairs) and does not use
#'   threads.
#'
#' @return `data.frame` with columns:
#'   * `motif_a`, `motif_b`: motif names (1-based indexing into the
#'     supplied `motifs` list; deduplicated names if duplicates exist).
#'   * `a_only`, `b_only`, `both`, `neither`: 2x2 contingency cells.
#'   * `odds_ratio`: sample odds ratio of the (pseudocount-adjusted)
#'     2x2, `(both * neither) / (a_only * b_only)`; `Inf` when an
#'     off-diagonal cell is zero. `NA` when `both < min.coocc`.
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
#' @seealso [scan_sequences_lite()], [enrich_motifs_lite()], [motif_peaks()]
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
                        no.overlaps  = FALSE,
                        no.overlaps.by = c("pvalue", "score"),
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
  if (!isTRUEorFALSE(no.overlaps))
    stop("`no.overlaps` must be a single logical", call. = FALSE)
  no.overlaps.by <- match.arg(no.overlaps.by)
  if (no.overlaps && !is.null(hits))
    warning("`no.overlaps` only applies to the internal scan; it is ignored ",
            "when a precomputed `hits` table is supplied.", call. = FALSE)

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
    ## Internal-scan path (DNA/RNA only -- enforced by scan_sequences_lite).
    hits <- scan_sequences_lite(motifs, sequences,
                            pvalue = pvalue, RC = RC,
                            no.overlaps = no.overlaps,
                            no.overlaps.by = no.overlaps.by,
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

  ## --- for spatial mode: build global-coordinate hit lists -------------
  ## The descriptive spatial columns need, for every co-occurring sequence,
  ## the nearest within-sequence (A, B) hit pair. Rather than loop
  ## sequence-by-sequence (and build an outer() distance matrix each time),
  ## map every hit to a single global coordinate
  ##     g = sequence.i * BIG + (start - minstart)
  ## with BIG = 2 * (max(start) - min(start)) + 1. Two hits in the same
  ## sequence then differ by at most `rng = max(start) - min(start)`, while
  ## two hits in different sequences differ by at least BIG - rng > rng. So
  ## the nearest hit by global value is always in the same sequence, which
  ## lets a single sorted findInterval()/diff() per motif pair resolve every
  ## within-sequence nearest pair with no per-sequence loop. Coordinates are
  ## kept as doubles (the products can exceed the 2^31 integer range) but
  ## stay well within 2^53 for realistic inputs.
  g_pos <- g_seq <- NULL
  if (!is.null(max.distance)) {
    max.distance <- as.numeric(max.distance)
    minstart <- as.numeric(min(hits$start))
    rng <- as.numeric(max(hits$start)) - minstart
    BIG <- 2 * rng + 1
    g_pos <- vector("list", n.motifs)   # global coords, sorted ascending
    g_seq <- vector("list", n.motifs)   # matching sequence ids
    for (i in seq_len(n.motifs)) {
      sub <- hits[hits$motif.i == i, , drop = FALSE]
      if (nrow(sub) > 0L) {
        g <- as.numeric(sub$sequence.i) * BIG +
             (as.numeric(sub$start) - minstart)
        o <- order(g)
        g_pos[[i]] <- g[o]
        g_seq[[i]] <- as.integer(sub$sequence.i)[o]
      } else {
        g_pos[[i]] <- numeric(0)
        g_seq[[i]] <- integer(0)
      }
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

  ## --- contingency counts (vectorised over all pairs) ------------------
  ## `both[i,j]` (sequences containing BOTH motifs) is exactly the (i,j)
  ## entry of crossprod(incidence), where `incidence` is the sequence x
  ## motif presence matrix. A single matrix product therefore gives every
  ## pair's count, replacing the old per-pair intersect()/fisher.test()
  ## loop that dominated the runtime once there were many motif
  ## combinations. The 2x2 partition is always the SET question (does the
  ## pair share sequences more than chance?) and is computed without the
  ## spatial filter, so a_only + b_only + both + neither == n.sequences.
  i_idx <- pair_idx[, 1L]
  j_idx <- pair_idx[, 2L]
  cnt <- lengths(presence)                 # sequences containing each motif

  ## Restrict the incidence matrix to sequences that carry at least one
  ## hit: hit-free sequences add nothing to any `both` count, and this
  ## keeps the matrix small when n.sequences far exceeds the number of
  ## sequences actually hit.
  present_ids <- sort(unique(unlist(presence, use.names = FALSE)))
  if (length(present_ids)) {
    incidence <- matrix(0, nrow = length(present_ids), ncol = n.motifs)
    for (i in seq_len(n.motifs))
      incidence[match(presence[[i]], present_ids), i] <- 1
    coocc <- crossprod(incidence)          # M[i,j] = sequences with both
    both_v <- as.integer(coocc[pair_idx])  # M[i,i] == cnt[i] for self-pairs
  } else {
    both_v <- integer(rows)
  }
  ## For self-pairs M[i,i] == cnt[i], so a_only == b_only == 0 falls out
  ## of the same formula with no special-casing.
  a_only_v  <- cnt[i_idx] - both_v
  b_only_v  <- cnt[j_idx] - both_v
  neither_v <- n.sequences - a_only_v - b_only_v - both_v
  neither_v[neither_v < 0L] <- 0L          # numerical guard

  out <- data.frame(
    motif_a    = mot.names[i_idx],
    motif_b    = mot.names[j_idx],
    a_only     = as.integer(a_only_v),
    b_only     = as.integer(b_only_v),
    both       = as.integer(both_v),
    neither    = as.integer(neither_v),
    odds_ratio = NA_real_,
    pvalue     = NA_real_,
    qvalue     = NA_real_,
    stringsAsFactors = FALSE
  )

  ## --- vectorised one-sided Fisher's exact test ------------------------
  ## For a 2x2 table the one-sided ("greater") Fisher exact p-value is
  ## exactly the upper tail of the hypergeometric distribution of the
  ## top-left cell given the margins, so phyper() returns every pair's
  ## p-value in a single call. With an all-integer table (the default
  ## pseudocount = 0) this matches fisher.test() to floating-point
  ## precision.
  pc <- as.numeric(pseudocount)
  a  <- out$both    + pc        # top-left cell: sequences with both motifs
  b  <- out$a_only  + pc
  cc <- out$b_only  + pc
  d  <- out$neither + pc
  tested <- out$both >= as.integer(min.coocc)
  out$pvalue[tested] <- stats::phyper(a[tested] - 1,
                                      a[tested] + b[tested],
                                      cc[tested] + d[tested],
                                      a[tested] + cc[tested],
                                      lower.tail = FALSE)
  ## Sample odds ratio of the (pseudocount-adjusted) 2x2. This is `Inf`
  ## when an off-diagonal cell is zero, matching the boundary behaviour
  ## of the conditional-MLE estimate fisher.test() used to report.
  out$odds_ratio[tested] <- (a[tested] * d[tested]) / (b[tested] * cc[tested])

  ## --- descriptive spatial columns (opt-in; do NOT feed the test) ------
  ## Computed only when `max.distance` is supplied and only over pairs that
  ## co-occur:
  ##   - both.clustered:  number of co-occurring sequences whose nearest
  ##                      (A, B) hit pair is within `max.distance` bp.
  ##   - median.distance: median of those per-sequence nearest distances.
  ## For each pair this is fully vectorised over hits using the global
  ## coordinates built above: findInterval() finds each A-hit's nearest
  ## B-hit (cross-pairs), or diff() of the sorted coordinates gives adjacent
  ## gaps (self-pairs); a sorted group-min then collapses to one distance
  ## per sequence. A proper spatial significance test (e.g. SpaMo's
  ## permutation of the distance distribution) is out of scope here.
  if (!is.null(max.distance)) {
    both.clustered  <- integer(rows)
    median.distance <- rep(NA_real_, rows)
    for (k in which(both_v > 0L)) {
      i <- i_idx[k]; j <- j_idx[k]
      if (i == j) {
        ## Self-pair: minimum adjacent gap among the motif's hits in each
        ## sequence. A sequence with a single hit yields only a cross-
        ## sequence diff (>= BIG), which is dropped, so it needs >= 2 hits
        ## to count -- matching the original behaviour.
        gi <- g_pos[[i]]; si <- g_seq[[i]]
        if (length(gi) < 2L) next
        d    <- diff(gi)
        same <- si[-length(si)] == si[-1L]
        if (!any(same)) next
        seqd <- si[-length(si)][same]
        dd   <- d[same]
      } else {
        ## Cross-pair: each A-hit's nearest B-hit (guaranteed same sequence
        ## by the BIG offset), then the per-sequence minimum.
        coocc <- intersect(presence[[i]], presence[[j]])
        if (!length(coocc)) next
        ina <- g_seq[[i]] %in% coocc
        inb <- g_seq[[j]] %in% coocc
        ga  <- g_pos[[i]][ina]; seqd <- g_seq[[i]][ina]
        gb  <- g_pos[[j]][inb]                        # sorted ascending
        len <- length(gb)
        pos <- findInterval(ga, gb)                   # 0 .. len
        below <- ifelse(pos >= 1L, ga - gb[pmax(pos, 1L)],      Inf)
        above <- ifelse(pos < len, gb[pmin(pos + 1L, len)] - ga, Inf)
        dd   <- pmin(below, above)                    # nearest dist per A-hit
      }
      ## Sorted group-min: one nearest distance per co-occurring sequence.
      o    <- order(seqd, dd)
      dmin <- dd[o][!duplicated(seqd[o])]
      qual <- dmin[dmin <= max.distance]
      both.clustered[k] <- length(qual)
      if (length(qual)) median.distance[k] <- stats::median(qual)
    }
    out$both.clustered  <- both.clustered
    out$median.distance <- median.distance
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
