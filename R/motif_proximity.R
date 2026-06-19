#' Test for motif enrichment near a set of genomic anchors.
#'
#' `motif_proximity()` tests whether a motif's hits sit closer to a set
#' of anchor coordinates than expected by chance, genome-wide. It is the
#' genome-wide sibling of [motif_peaks()]: where [motif_peaks()] asks
#' whether hits cluster in a central window of equal-length sequences,
#' `motif_proximity()` asks whether hits cluster near anchors scattered
#' across whole chromosomes. The default test mirrors the CentriMo-style
#' binomial of [motif_peaks()], but the success probability is the
#' fraction of the scanned genome lying within a distance of an anchor
#' rather than the fraction of a sequence lying in a window.
#'
#' Three nulls are available via `method`:
#' \describe{
#'   \item{`"binomial"`}{(default) For each candidate distance `d`, the
#'     "near zone" is the set of genomic positions within `d` of any
#'     anchor; `f_d` is its fraction of the universe `U`. Under the null
#'     that hits are uniformly distributed over `U`, the count `k` of
#'     hits inside the near zone follows `Binomial(N, f_d)`. The best `d`
#'     per motif is reported, Bonferroni-corrected over the candidate
#'     distances, then BH-adjusted across motifs. Fast, but the
#'     uniform-hit null cannot correct composition bias.}
#'   \item{`"permutation"`}{The anchors (not the hits) are randomised
#'     within `U` `n.perm` times; the binomial z-score is recomputed each
#'     time to build an empirical null. Randomising the anchors keeps the
#'     real hits (and their composition-driven clustering) fixed, so the
#'     test asks whether hits are closer to the real anchors than to
#'     random ones, which absorbs composition and mappability bias. This
#'     is the recommended method when composition bias is a concern.
#'     Uses the global RNG, so call [set.seed()] beforehand for
#'     reproducible p-values (see `nthreads`).}
#'   \item{`"ks"`}{A threshold-free one-sided Kolmogorov-Smirnov test of
#'     the per-hit distance-to-nearest-anchor against its analytical null
#'     CDF (`f` as a function of `d`). P-values are approximate under the
#'     integer-distance ties; `"permutation"` is the exact alternative.
#'     Only available for `count = "hits"` and
#'     `orientation = "undirected"`.}
#' }
#'
#' The counting unit is set by `count`. `"hits"` (default) collapses
#' overlapping same-motif hits to one best-scoring representative (the
#' genome-wide analogue of the best-hit-per-sequence convention of
#' [motif_peaks()]) and counts hits in the near zone. `"anchors"` counts
#' the fraction of anchors that have at least one nearby hit; it has no
#' clean analytical null, so `method = "permutation"` is recommended and
#' `method = "binomial"` falls back to a Poisson occupancy approximation
#' (anti-conservative when hits cluster).
#'
#' @param hits Motif hits from [scan_sequences_lite()] or [scan_sequences()],
#'   as a `GRanges` (preferred) or a `data.frame`, carrying genomic
#'   coordinates plus `motif` and `score`. Supply this *or*
#'   (`motifs` + `sequences`), not both. `data.frame` input also requires
#'   `universe` (its seqlengths are unknown).
#' @param anchors `GRanges`. The anchor coordinates to test proximity to.
#'   Anchor width matters: the near zone is the anchor widened by the
#'   distance on each side, so wide anchors enlarge the within-distance
#'   fraction and collapse the resolution of the test (see Caveats). For a
#'   point hypothesis (summits, TSSs, SNPs, motif centres) use 1 bp
#'   anchors, reducing regions first with e.g.
#'   `GenomicRanges::resize(anchors, width = 1, fix = "center")`.
#' @param motifs,sequences Alternative to `hits`: motifs (any format
#'   accepted by [convert_motifs()]) and the sequences to scan. `sequences`
#'   is an `XStringSet` (the scanned chromosomes) or a `BSgenome`; for a
#'   `BSgenome` only the chromosomes carrying anchors are extracted and
#'   scanned. Scanning is done internally with [scan_sequences_lite()].
#' @param universe `GRanges` or `NULL`. The regions actually scanned (the
#'   null support). `NULL` (default) derives whole chromosomes from
#'   `seqlengths(hits)`. **Pass this explicitly for partial or masked
#'   scans**, otherwise the test is anti-conservative.
#' @param method `character(1)`. `"binomial"` (default), `"permutation"`,
#'   or `"ks"`. See Details.
#' @param count `character(1)`. `"hits"` (default) or `"anchors"`. See
#'   Details.
#' @param orientation `character(1)`. `"undirected"` (default; symmetric
#'   near zone) or `"upstream"` / `"downstream"` for a strand-aware,
#'   directional near zone built with [GenomicRanges::flank()]. Directional
#'   modes require stranded anchors and are not available for
#'   `method = "ks"`.
#' @param scan.pvalue `numeric(1)`. Hit P-value threshold forwarded to
#'   [scan_sequences_lite()] when scanning internally. Default `1e-4`.
#' @param RC `logical(1)`. Scan both strands when scanning internally.
#'   Default `TRUE`.
#' @param dist.thresholds `integer` or `NULL`. Candidate distances `d` to
#'   test. `NULL` (default) builds a log-spaced grid capped at the
#'   smallest chromosome. Ignored by `method = "ks"`.
#' @param n.perm `integer(1)`. Number of permutations for
#'   `method = "permutation"`. Default `1000`.
#' @param perm.per.chrom `logical(1)`. Randomise anchors within their own
#'   chromosome (`TRUE`, default, more conservative) or across the whole
#'   universe.
#' @param ignore.strand `logical(1)`. Match hits to anchors ignoring
#'   strand. Default `TRUE`. The near-zone *geometry* is always unstranded.
#' @param qvalue `numeric(1)`. BH q-value cutoff for the reported motifs.
#'   Default `0.1`. Set to `1` to return every motif.
#' @param nthreads `numeric(1)`. Threads for internal scanning and the
#'   per-motif loop. `0` uses all cores. Default `1`. For
#'   `method = "permutation"`, exact reproducibility needs `nthreads = 1`
#'   together with a prior [set.seed()]; parallel workers draw independent
#'   streams that are not bit-reproducible unless `RNGkind("L'Ecuyer-CMRG")`
#'   is in effect.
#'
#' @return A `data.frame` with one row per motif passing the q-value
#'   cutoff, sorted by `pvalue` ascending. Columns: `motif`, `motif.i`,
#'   `method`, `count`, `orientation`, `nhits`, `best.dist`, `hits.near`,
#'   `hits.far`, `frac.near`, `expected.near`, `enrichment`,
#'   `log2.enrichment`, `pvalue`, `qvalue`, `ks.D` (non-`NA` only for
#'   `method = "ks"`), and `dist.to.nearest` (list-column of per-hit
#'   signed distances to the nearest anchor, consumed by
#'   [plot_motif_proximity()]). For `count = "anchors"`, `nhits`,
#'   `hits.near` and `hits.far` count anchors rather than hits.
#'
#' @section Caveats:
#'
#' The universe derived from `seqlengths(hits)` is correct only for whole,
#' unmasked chromosome scans. For masked or partial scans (or for windowed
#' sequences extracted with [Biostrings::getSeq()]), pass `universe`
#' explicitly; otherwise `f_d` is overstated and the test is
#' anti-conservative.
#'
#' The `"binomial"` null assumes hits are uniformly distributed over the
#' universe and cannot correct for composition or mappability bias. Prefer
#' `method = "permutation"` when that bias is plausible.
#'
#' Anchor width strongly affects power and resolution. The near zone is the
#' anchor widened by `d` on each side, so its width is `width(anchor) + 2d`
#' and the within-distance fraction `f_d` grows with anchor width; a hit
#' inside an anchor region also has a distance of zero. Wide anchors
#' therefore inflate the expected count (shrinking the enrichment toward
#' 1), wash out the distance distribution, and weaken the test, while the
#' biological signal is unchanged. For a point hypothesis (summits, TSSs,
#' SNPs, motif centres) use 1 bp anchors, reducing regions with
#' `GenomicRanges::resize(anchors, width = 1, fix = "center")` (or the
#' called summit). Keep regions only when the question really is
#' enrichment within or around intervals, and then scale `dist.thresholds`
#' to the region width.
#'
#' All inputs must share chromosome names. A mismatch (e.g. `chr1` versus
#' `Chr1` versus `1`) raises an error; reconcile with
#' [GenomeInfoDb::seqlevelsStyle()].
#'
#' `method = "permutation"` rebuilds the near zone for every distance
#' threshold on every permutation, so its cost scales with `n.perm` times
#' `length(dist.thresholds)`; for genome-scale anchor sets a run can take
#' minutes per motif. Lower `n.perm`, shorten `dist.thresholds`, or use
#' `method = "binomial"` / `"ks"` when that is prohibitive.
#'
#' @references
#'
#' Bailey TL, Machanick P (2012). "Inferring direct DNA binding from
#' ChIP-seq." *Nucleic Acids Research*, **40**(17):e128.
#' \doi{10.1093/nar/gks433}.
#'
#' @examples
#' \dontrun{
#' library(universalmotif)
#' library(GenomicRanges)
#' ## Treat three created sequences as small chromosomes.
#' set.seed(1)
#' chrs <- do.call(c, lapply(1:3, function(i)
#'   create_sequences(seqnum = 1, seqlen = 5000, rng.seed = i)))
#' names(chrs) <- paste0("chr", 1:3)
#'
#' ## Plant a motif near a set of anchors on chr1.
#' anchor.pos <- seq(500, 4500, by = 400)
#' planted <- Biostrings::DNAString("TTGACATA")
#' for (p in anchor.pos)
#'   subseq(chrs[["chr1"]], start = p + 5L, width = 8) <- planted
#' anchors <- GRanges("chr1", IRanges(anchor.pos, width = 1))
#'
#' m <- create_motif("TTGACATA", name = "example")
#'
#' ## Either scan first and pass hits ...
#' hits <- scan_sequences_lite(m, chrs, pvalue = 1e-3, return.granges = TRUE)
#' motif_proximity(hits, anchors)
#'
#' ## ... or let motif_proximity() scan internally.
#' motif_proximity(motifs = m, sequences = chrs, anchors = anchors)
#'
#' ## A composition-robust permutation test (set.seed() for reproducibility).
#' set.seed(1)
#' motif_proximity(hits, anchors, method = "permutation", n.perm = 199)
#' }
#'
#' @seealso [motif_peaks()], [scan_sequences_lite()], [plot_motif_proximity()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
motif_proximity <- function(hits = NULL, anchors,
                            motifs = NULL, sequences = NULL,
                            universe = NULL,
                            method = c("binomial", "permutation", "ks"),
                            count  = c("hits", "anchors"),
                            orientation = c("undirected", "upstream", "downstream"),
                            scan.pvalue = 1e-4, RC = TRUE,
                            dist.thresholds = NULL, n.perm = 1000L,
                            perm.per.chrom = TRUE, ignore.strand = TRUE,
                            qvalue = 0.1, nthreads = 1) {

  ## --- arg validation --------------------------------------------------
  method      <- match.arg(method)
  count       <- match.arg(count)
  orientation <- match.arg(orientation)

  if (method == "ks" && count == "anchors")
    stop("method = \"ks\" is not available for count = \"anchors\"; ",
         "use method = \"permutation\".", call. = FALSE)
  if (method == "ks" && orientation != "undirected")
    stop("method = \"ks\" is not available for directional orientation; ",
         "use method = \"binomial\" or \"permutation\".", call. = FALSE)

  if (!is.numeric(qvalue) || length(qvalue) != 1L || is.na(qvalue) ||
      qvalue <= 0 || qvalue > 1)
    stop("`qvalue` must be a single numeric in (0, 1]", call. = FALSE)
  if (!is.numeric(n.perm) || length(n.perm) != 1L || n.perm < 1L)
    stop("`n.perm` must be a positive integer", call. = FALSE)
  n.perm <- as.integer(n.perm)
  if (!is.null(dist.thresholds) &&
      (!is.numeric(dist.thresholds) || anyNA(dist.thresholds) ||
       any(dist.thresholds <= 0)))
    stop("`dist.thresholds` must be positive integers", call. = FALSE)

  for (pkg in c("GenomicRanges", "GenomeInfoDb"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("`motif_proximity()` requires the ", pkg, " package", call. = FALSE)

  if (missing(anchors) || !inherits(anchors, "GRanges"))
    stop("`anchors` must be a GRanges", call. = FALSE)
  if (orientation != "undirected") {
    a.strand <- as.character(GenomicRanges::strand(anchors))
    if (!any(a.strand %in% c("+", "-")))
      stop("orientation = \"", orientation, "\" requires stranded anchors",
           call. = FALSE)
  }

  nthreads <- resolve_nthreads(nthreads)

  ## --- resolve / normalise hits + universe -----------------------------
  hits0   <- resolve_proximity_hits(hits, motifs, sequences, anchors,
                                    scan.pvalue, RC, nthreads)
  hits.gr <- normalize_motif_proximity_input(hits0)
  U       <- build_proximity_universe(universe, hits.gr)

  rr      <- restrict_to_anchored_universe(U, hits.gr, anchors)
  U       <- rr$U
  hits.gr <- rr$hits
  anchors <- rr$anchors
  Utot    <- sum(as.numeric(GenomicRanges::width(U)))

  if (length(hits.gr) == 0L || length(anchors) == 0L || Utot <= 0)
    return(empty_proximity_result())

  if (count == "hits")
    hits.gr <- collapse_hits_per_cluster(hits.gr, ignore.strand)

  if (is.null(dist.thresholds))
    dist.thresholds <- default_dist_thresholds(U)
  else
    dist.thresholds <- sort(unique(as.integer(dist.thresholds)))
  if (!length(dist.thresholds))
    stop("no usable `dist.thresholds`", call. = FALSE)

  if (count == "anchors" && method == "binomial")
    message("count = \"anchors\" with method = \"binomial\" uses a Poisson ",
            "occupancy approximation that assumes non-clustered hits; ",
            "method = \"permutation\" is preferred.")

  ## Pre-compute the null CDF once for the KS test (motif-independent).
  Ffun <- NULL
  if (method == "ks") {
    fc   <- proximity_f_curve(anchors, U, Utot, orientation)
    Ffun <- stats::approxfun(fc$d, fc$f, yleft = 0, yright = 1, rule = 2)
  }

  ## --- per-motif test --------------------------------------------------
  motif_ids <- sort(unique(as.character(S4Vectors::mcols(hits.gr)$motif)))

  one <- function(i) {
    hm <- hits.gr[as.character(S4Vectors::mcols(hits.gr)$motif) == motif_ids[i]]
    per_motif_proximity(hm, anchors, U, Utot, dist.thresholds, method, count,
                        orientation, n.perm, perm.per.chrom, ignore.strand,
                        Ffun)
  }
  per <- if (nthreads > 1L && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(seq_along(motif_ids), one, mc.cores = nthreads)
  } else {
    lapply(seq_along(motif_ids), one)
  }

  assemble_proximity_result(motif_ids, per, method, count, orientation, qvalue)
}

#' Plot the distance-to-anchor distribution of motif hits.
#'
#' Companion to [motif_proximity()]. Returns a `ggplot` faceted by motif,
#' showing a histogram of per-hit distances to the nearest anchor (signed
#' when the anchors are stranded, so upstream and downstream sit on
#' opposite sides of zero) with the most-enriched near zone shaded.
#'
#' @param proximity `data.frame` returned by [motif_proximity()]. Must
#'   carry the `dist.to.nearest` list-column.
#' @param motifs `character` or `NULL`. Subset to these motifs (by name).
#'   `NULL` (default) plots every row of `proximity`.
#' @param ncol `integer(1)`. Number of facet columns. Default `2`.
#' @param bins `integer(1)`. Histogram bin count. Default `50`.
#' @param fill.zone `character(1)`. Colour of the best-zone shade.
#'   Default `"#ffe4a3"` (a pale yellow).
#' @param log.x `logical(1)`. Plot `log10` of the absolute distance
#'   instead of the signed distance. Default `FALSE`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [motif_proximity()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_motif_proximity <- function(proximity, motifs = NULL, ncol = 2L,
                                 bins = 50L, fill.zone = "#ffe4a3",
                                 log.x = FALSE) {

  if (!is.data.frame(proximity))
    stop("`proximity` must be a data.frame from motif_proximity()",
         call. = FALSE)
  if (!"dist.to.nearest" %in% names(proximity))
    stop("`proximity` is missing the `dist.to.nearest` list-column ",
         "(was it produced by motif_proximity()?)", call. = FALSE)
  if (nrow(proximity) == 0L)
    stop("`proximity` has 0 rows -- no enriched motifs to plot",
         call. = FALSE)

  if (!is.null(motifs)) {
    keep <- proximity$motif %in% motifs
    if (!any(keep))
      stop("none of the requested motifs are in `proximity`", call. = FALSE)
    proximity <- proximity[keep, , drop = FALSE]
  }

  long <- do.call(rbind, lapply(seq_len(nrow(proximity)), function(i) {
    d <- proximity$dist.to.nearest[[i]]
    d <- d[is.finite(d)]
    if (log.x) d <- abs(d) + 1
    if (!length(d)) return(NULL)
    data.frame(motif = proximity$motif[i], distance = d,
               stringsAsFactors = FALSE)
  }))
  if (is.null(long) || !nrow(long))
    stop("no finite distances to plot", call. = FALSE)
  long$motif <- factor(long$motif, levels = proximity$motif)

  ## Best-zone rectangles (skip motifs without a best.dist, e.g. ks rows).
  rects <- proximity[!is.na(proximity$best.dist), , drop = FALSE]
  rect.df <- NULL
  if (nrow(rects)) {
    xmax <- if (log.x) log10(rects$best.dist + 1) else rects$best.dist
    xmin <- if (log.x) 0 else -rects$best.dist
    rect.df <- data.frame(
      motif = factor(rects$motif, levels = proximity$motif),
      xmin = xmin, xmax = xmax, stringsAsFactors = FALSE
    )
  }

  strip_labels <- setNames(
    sprintf("%s  (p = %.2g)", proximity$motif, proximity$pvalue),
    proximity$motif
  )

  g <- ggplot2::ggplot(long, ggplot2::aes(x = .data$distance))
  if (!is.null(rect.df))
    g <- g + ggplot2::geom_rect(data = rect.df,
                                ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                             ymin = -Inf, ymax = Inf),
                                fill = fill.zone, alpha = 0.5,
                                inherit.aes = FALSE)
  g <- g +
    ggplot2::geom_histogram(bins = bins, colour = NA, fill = "black") +
    ggplot2::facet_wrap(~ motif, ncol = ncol, scales = "free_y",
                        labeller = ggplot2::as_labeller(strip_labels)) +
    ggplot2::labs(x = if (log.x) "Distance to nearest anchor (bp, log10)"
                      else "Signed distance to nearest anchor (bp)",
                  y = "Count") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background    = ggplot2::element_rect(fill = NA, colour = NA),
      panel.grid          = ggplot2::element_blank(),
      panel.border        = ggplot2::element_blank(),
      axis.line.x.bottom  = ggplot2::element_line(colour = "black"),
      axis.line.y.left    = ggplot2::element_line(colour = "black")
    )
  if (log.x) g <- g + ggplot2::scale_x_log10()
  g
}

## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

empty_proximity_result <- function() {
  data.frame(
    motif           = character(0),
    motif.i         = integer(0),
    method          = character(0),
    count           = character(0),
    orientation     = character(0),
    nhits           = integer(0),
    best.dist       = integer(0),
    hits.near       = integer(0),
    hits.far        = integer(0),
    frac.near       = numeric(0),
    expected.near   = numeric(0),
    enrichment      = numeric(0),
    log2.enrichment = numeric(0),
    pvalue          = numeric(0),
    qvalue          = numeric(0),
    ks.D            = numeric(0),
    dist.to.nearest = I(list()),
    stringsAsFactors = FALSE
  )
}

## Dispatch the two input modes: precomputed `hits`, or scan `motifs`
## against `sequences` (XStringSet or BSgenome) internally.
resolve_proximity_hits <- function(hits, motifs, sequences, anchors,
                                   scan.pvalue, RC, nthreads) {
  have_hits <- !is.null(hits)
  have_scan <- !is.null(motifs) || !is.null(sequences)
  if (have_hits && have_scan)
    stop("supply either `hits` or (`motifs` + `sequences`), not both",
         call. = FALSE)
  if (!have_hits && !have_scan)
    stop("supply either `hits` or (`motifs` + `sequences`)", call. = FALSE)

  if (have_hits) return(hits)

  if (is.null(motifs) || is.null(sequences))
    stop("scanning internally needs both `motifs` and `sequences`",
         call. = FALSE)
  sequences <- proximity_get_sequences(sequences, anchors)
  out <- scan_sequences_lite(motifs, sequences, pvalue = scan.pvalue, RC = RC,
                         nthreads = nthreads, return.granges = TRUE)
  if (!inherits(out, "GRanges"))
    stop("internal scan did not return a GRanges; is GenomicRanges installed?",
         call. = FALSE)
  out
}

## Turn `sequences` into an XStringSet to scan. A BSgenome is reduced to
## the chromosomes carrying anchors (the rest are dropped by U* anyway).
proximity_get_sequences <- function(sequences, anchors) {
  if (methods::is(sequences, "BSgenome")) {
    for (pkg in c("BSgenome", "GenomeInfoDb"))
      if (!requireNamespace(pkg, quietly = TRUE))
        stop("scanning a BSgenome requires the ", pkg, " package",
             call. = FALSE)
    chrs <- intersect(as.character(GenomeInfoDb::seqlevelsInUse(anchors)),
                      GenomeInfoDb::seqnames(sequences))
    if (!length(chrs))
      stop("none of the anchor seqnames are present in the genome",
           call. = FALSE)
    seqs <- Biostrings::getSeq(sequences, names = chrs)
    names(seqs) <- chrs
    seqs
  } else if (methods::is(sequences, "XStringSet")) {
    sequences
  } else {
    stop("`sequences` must be an XStringSet or a BSgenome object",
         call. = FALSE)
  }
}

## Coerce hits (GRanges or data.frame) to a GRanges with motif/score mcols.
normalize_motif_proximity_input <- function(hits) {

  if (inherits(hits, "GRanges")) {
    mc <- S4Vectors::mcols(hits)
    if (!"motif" %in% names(mc))
      stop("GRanges `hits` must have a `motif` mcols column", call. = FALSE)
    if (!"score" %in% names(mc))
      stop("GRanges `hits` must have a `score` mcols column", call. = FALSE)
    gr <- hits
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
      motif = as.character(mc$motif),
      score = as.numeric(mc$score)
    )
    return(gr)
  }

  if (is.data.frame(hits)) {
    if (!"motif" %in% names(hits))
      stop("data.frame `hits` must have a `motif` column", call. = FALSE)
    if (!"score" %in% names(hits))
      stop("data.frame `hits` must have a `score` column", call. = FALSE)
    if (!"start" %in% names(hits))
      stop("data.frame `hits` must have a `start` column", call. = FALSE)
    end_v <- if ("end" %in% names(hits)) hits$end
             else if ("stop" %in% names(hits)) hits$stop
             else stop("data.frame `hits` must have an `end` or `stop` column",
                       call. = FALSE)
    seq_v <- if ("sequence" %in% names(hits)) hits$sequence
             else if ("sequence.i" %in% names(hits)) as.character(hits$sequence.i)
             else stop("data.frame `hits` must have a `sequence` or ",
                       "`sequence.i` column", call. = FALSE)
    strand_v <- if ("strand" %in% names(hits)) as.character(hits$strand)
                else "*"
    gr <- GenomicRanges::GRanges(
      seqnames = as.character(seq_v),
      ranges   = IRanges::IRanges(start = as.integer(hits$start),
                                  end   = as.integer(end_v)),
      strand   = strand_v
    )
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
      motif = as.character(hits$motif),
      score = as.numeric(hits$score)
    )
    return(gr)
  }

  stop("`hits` must be a data.frame or GRanges", call. = FALSE)
}

## Build the reduced universe U. Derive whole chromosomes from
## seqlengths(hits) when `universe` is NULL.
build_proximity_universe <- function(universe, hits.gr) {
  if (!is.null(universe)) {
    if (!inherits(universe, "GRanges"))
      stop("`universe` must be a GRanges or NULL", call. = FALSE)
    return(GenomicRanges::reduce(universe, ignore.strand = TRUE))
  }
  sl <- GenomeInfoDb::seqlengths(hits.gr)
  sl <- sl[!is.na(sl) & sl > 0L]
  if (!length(sl))
    stop("could not derive a universe from seqlengths(hits); pass ",
         "`universe` explicitly (required for data.frame input or for ",
         "hits without seqlengths).", call. = FALSE)
  GenomicRanges::GRanges(
    seqnames = names(sl),
    ranges   = IRanges::IRanges(start = 1L, end = as.integer(sl))
  )
}

## Restrict to U* = chromosomes carrying >=1 anchor; reconcile seqlevels
## across hits, anchors and U; harmonise seqinfo to U*.
restrict_to_anchored_universe <- function(U, hits.gr, anchors.gr) {

  a.levels <- as.character(GenomeInfoDb::seqlevelsInUse(anchors.gr))
  common   <- intersect(a.levels, GenomeInfoDb::seqlevels(U))
  if (!length(common))
    stop("anchors and the scanned universe share no common seqnames.\n",
         "  anchors:  ", paste(utils::head(a.levels, 5L), collapse = ", "), "\n",
         "  universe: ", paste(utils::head(GenomeInfoDb::seqlevels(U), 5L),
                               collapse = ", "), "\n",
         "  Do they use the same chromosome naming? ",
         "See GenomeInfoDb::seqlevelsStyle().", call. = FALSE)

  n.anchor0 <- length(anchors.gr)
  n.hit0    <- length(hits.gr)

  U.star <- GenomeInfoDb::keepSeqlevels(
    U, common, pruning.mode = "coarse")
  U.star <- GenomicRanges::reduce(U.star, ignore.strand = TRUE)

  anchors.star <- GenomeInfoDb::keepSeqlevels(
    anchors.gr, common, pruning.mode = "coarse")

  hits.in <- GenomeInfoDb::keepSeqlevels(
    hits.gr[as.character(GenomicRanges::seqnames(hits.gr)) %in% common],
    common, pruning.mode = "coarse")

  ## Harmonise seqinfo so derived zones share U*'s seqlengths (lets trim()
  ## clamp dilated zones and lets intersect() merge cleanly).
  lv <- GenomeInfoDb::seqlevels(U.star)
  si <- GenomeInfoDb::seqinfo(U.star)
  suppressWarnings({
    GenomeInfoDb::seqlevels(anchors.star) <- lv
    GenomeInfoDb::seqinfo(anchors.star)   <- si
    GenomeInfoDb::seqlevels(hits.in)      <- lv
    GenomeInfoDb::seqinfo(hits.in)        <- si
  })

  hits.star <- IRanges::subsetByOverlaps(hits.in, U.star, ignore.strand = TRUE)

  n.hit.drop    <- n.hit0 - length(hits.star)
  n.anchor.drop <- n.anchor0 - length(anchors.star)
  if (n.hit.drop > 0L)
    message(n.hit.drop, " hit(s) outside the anchor-bearing universe ",
            "were dropped")
  if (n.anchor.drop > 0L)
    message(n.anchor.drop, " anchor(s) on chromosomes absent from the ",
            "universe were dropped")

  list(U = U.star, hits = hits.star, anchors = anchors.star)
}

## Keep one best-scoring hit per connected overlap cluster, per motif.
## Genome-wide analogue of best_hit_per_seq() in motif_peaks.R.
collapse_hits_per_cluster <- function(hits.gr, ignore.strand) {
  if (length(hits.gr) <= 1L) return(hits.gr)
  score <- as.numeric(S4Vectors::mcols(hits.gr)$score)
  motif <- as.character(S4Vectors::mcols(hits.gr)$motif)
  keep  <- logical(length(hits.gr))
  for (m in unique(motif)) {
    idx <- which(motif == m)
    sub <- hits.gr[idx]
    clusters <- GenomicRanges::reduce(sub, ignore.strand = ignore.strand,
                                      min.gapwidth = 0L)
    ov  <- GenomicRanges::findOverlaps(sub, clusters,
                                       ignore.strand = ignore.strand)
    cl  <- integer(length(sub))
    cl[S4Vectors::queryHits(ov)] <- S4Vectors::subjectHits(ov)
    sc  <- score[idx]
    o   <- order(cl, -sc)
    first <- !duplicated(cl[o])
    keep[idx[o[first]]] <- TRUE
  }
  hits.gr[keep]
}

## Default log-spaced distance grid, capped at half the smallest chromosome.
default_dist_thresholds <- function(U) {
  segw <- as.numeric(GenomicRanges::width(U))
  cap  <- max(10L, as.integer(floor(min(segw) / 2)))
  grid <- c(10L, 25L, 50L, 100L, 250L, 500L, 1000L, 2500L, 5000L,
            10000L, 25000L, 50000L)
  grid <- grid[grid <= cap]
  if (!length(grid)) grid <- cap
  as.integer(unique(grid))
}

## Per-anchor d-neighbourhood, strand stripped. Symmetric for
## "undirected"; a strand-aware one-sided flank for directional modes.
dilate_anchors <- function(anchors, d, orientation) {
  if (orientation == "undirected") {
    z <- anchors
    GenomicRanges::strand(z) <- "*"
    GenomicRanges::ranges(z) <- IRanges::IRanges(
      start = GenomicRanges::start(z) - d,
      end   = GenomicRanges::end(z)   + d)
    z <- suppressWarnings(GenomicRanges::trim(z))
  } else {
    z <- suppressWarnings(
      GenomicRanges::flank(anchors, width = d,
                           start = orientation == "upstream"))
    GenomicRanges::strand(z) <- "*"
    z <- suppressWarnings(GenomicRanges::trim(z))
  }
  z
}

## Scan the distance grid for one motif. Picks the most significant window
## (smallest binomial p, as motif_peaks() does), and also returns its
## z-score (the statistic used by the permutation null) and the number of
## thresholds actually evaluated (for the Bonferroni correction). Works for
## both count modes.
proximity_scan <- function(hits.m, anchors, U, Utot, dthr, count, orientation) {
  N  <- if (count == "hits") length(hits.m) else length(anchors)
  nh <- length(hits.m)
  best <- list(z = -Inf, p = 1, d = NA_integer_, k = 0L, N = N,
               expct = NA_real_, frac = NA_real_, enr = NA_real_, n.tested = 0L)
  if (N == 0L || nh == 0L) return(best)
  lambda   <- nh / Utot
  n.tested <- 0L

  for (d in dthr) {
    if (count == "hits") {
      zone <- GenomicRanges::reduce(dilate_anchors(anchors, d, orientation),
                                    ignore.strand = TRUE)
      zone <- GenomicRanges::intersect(zone, U, ignore.strand = TRUE)
      f <- sum(as.numeric(GenomicRanges::width(zone))) / Utot
      if (!is.finite(f) || f <= 0 || f >= 1) next
      k <- sum(IRanges::overlapsAny(hits.m, zone, ignore.strand = TRUE))
    } else {
      za   <- dilate_anchors(anchors, d, orientation)
      nbhd <- as.numeric(GenomicRanges::width(za))
      g    <- 1 - exp(-lambda * nbhd)
      f    <- mean(g)
      if (!is.finite(f) || f <= 0 || f >= 1) next
      k <- sum(IRanges::overlapsAny(za, hits.m, ignore.strand = TRUE))
    }
    n.tested <- n.tested + 1L
    expct <- N * f
    p     <- stats::pbinom(k - 1L, N, f, lower.tail = FALSE)
    z     <- (k - expct) / sqrt(N * f * (1 - f))
    ## Select by smallest p; break ties by larger z.
    if (p < best$p || (p == best$p && is.finite(z) && z > best$z)) {
      best <- list(z = z, p = p, d = as.integer(d), k = as.integer(k), N = N,
                   expct = expct, frac = k / N,
                   enr = if (expct > 0) k / expct else NA_real_,
                   n.tested = n.tested)
    }
  }
  best$n.tested <- n.tested
  best
}

## Null CDF f(d) over a grid, for the KS test.
proximity_f_curve <- function(anchors, U, Utot, orientation, npts = 200L) {
  dmax <- max(as.numeric(GenomicRanges::width(U)))
  grid <- unique(c(0, as.integer(round(
    exp(seq(log(1), log(max(dmax, 2)), length.out = npts))))))
  fvals <- vapply(grid, function(d) {
    if (d <= 0) return(0)
    zone <- GenomicRanges::reduce(dilate_anchors(anchors, d, orientation),
                                  ignore.strand = TRUE)
    zone <- GenomicRanges::intersect(zone, U, ignore.strand = TRUE)
    sum(as.numeric(GenomicRanges::width(zone))) / Utot
  }, numeric(1))
  list(d = grid, f = pmin(pmax(fvals, 0), 1))
}

## Per-hit signed distance to the nearest anchor (sign in feature frame).
proximity_signed_distances <- function(hits.m, anchors, ignore.strand) {
  out <- rep(NA_real_, length(hits.m))
  if (!length(hits.m) || !length(anchors)) return(out)
  d2n <- GenomicRanges::distanceToNearest(hits.m, anchors,
                                          ignore.strand = TRUE)
  if (length(d2n) == 0L) return(out)
  qi <- S4Vectors::queryHits(d2n)
  si <- S4Vectors::subjectHits(d2n)
  dd <- as.numeric(S4Vectors::mcols(d2n)$distance)
  hc <- (GenomicRanges::start(hits.m) + GenomicRanges::end(hits.m)) / 2
  ac <- (GenomicRanges::start(anchors) + GenomicRanges::end(anchors)) / 2
  raw <- hc[qi] - ac[si]
  ss  <- as.character(GenomicRanges::strand(anchors))[si]
  strand.sign <- ifelse(ss == "-", -1, 1)
  off.sign <- sign(raw)
  off.sign[off.sign == 0] <- 1
  out[qi] <- dd * off.sign * strand.sign
  out
}

## Sample a random anchor set inside U, preserving count and widths.
## Mirrors sample_genomic_universe() in match_bkg.R.
sample_anchors_in_universe <- function(anchors, U, per.chrom) {
  seg.chr   <- as.character(GenomicRanges::seqnames(U))
  seg.start <- as.numeric(GenomicRanges::start(U))
  seg.width <- as.numeric(GenomicRanges::width(U))
  a.w   <- as.numeric(GenomicRanges::width(anchors))
  a.chr <- as.character(GenomicRanges::seqnames(anchors))
  n     <- length(anchors)
  new.start <- numeric(n)
  new.chr   <- character(n)

  ## For each anchor of width w, pick a segment with probability
  ## proportional to its number of legal start positions (seg.width - w + 1),
  ## then a uniform start within. Group by width to vectorise.
  place <- function(ii, seg.idx) {
    sw <- seg.width[seg.idx]; ss <- seg.start[seg.idx]; sc <- seg.chr[seg.idx]
    out.chr <- character(length(ii)); out.start <- numeric(length(ii))
    wii <- a.w[ii]
    for (w in unique(wii)) {
      pos   <- which(wii == w)
      legal <- pmax(0, sw - w + 1)
      if (!any(legal > 0))
        stop("could not place anchors in the universe (a segment is ",
             "narrower than an anchor of width ", w, "); U too small.",
             call. = FALSE)
      sel <- sample.int(length(seg.idx), length(pos), replace = TRUE,
                        prob = legal / sum(legal))
      off <- floor(stats::runif(length(pos)) * legal[sel])
      out.chr[pos]   <- sc[sel]
      out.start[pos] <- ss[sel] + off
    }
    list(chr = out.chr, start = out.start)
  }

  if (per.chrom) {
    for (ch in unique(a.chr)) {
      ii  <- which(a.chr == ch)
      sij <- which(seg.chr == ch)
      if (!length(sij))
        stop("no universe segment on chromosome ", ch, call. = FALSE)
      pl <- place(ii, sij)
      new.chr[ii]   <- pl$chr
      new.start[ii] <- pl$start
    }
  } else {
    pl <- place(seq_len(n), seq_along(seg.chr))
    new.chr   <- pl$chr
    new.start <- pl$start
  }

  GenomicRanges::GRanges(
    seqnames = new.chr,
    ranges   = IRanges::IRanges(start = as.integer(new.start),
                                width = as.integer(a.w)),
    strand   = GenomicRanges::strand(anchors),
    seqinfo  = GenomeInfoDb::seqinfo(U)
  )
}

## One motif: dispatch by method, return a uniform per-motif list.
per_motif_proximity <- function(hits.m, anchors, U, Utot, dthr, method, count,
                                orientation, n.perm, perm.per.chrom,
                                ignore.strand, Ffun) {

  res <- list(
    nhits = if (count == "hits") length(hits.m) else length(anchors),
    best.dist = NA_integer_, hits.near = NA_integer_, hits.far = NA_integer_,
    frac.near = NA_real_, expected.near = NA_real_, enrichment = NA_real_,
    log2.enrichment = NA_real_, pvalue = 1, ks.D = NA_real_,
    dist.to.nearest = proximity_signed_distances(hits.m, anchors, ignore.strand)
  )
  if (length(hits.m) == 0L) return(res)

  ## Copy the per-window descriptive fields from a proximity_scan() result
  ## into `res`, alongside the supplied p-value.
  fill_from_scan <- function(res, sc, pval) {
    res$best.dist       <- sc$d
    res$hits.near       <- sc$k
    res$hits.far        <- as.integer(sc$N - sc$k)
    res$frac.near       <- sc$frac
    res$expected.near   <- sc$expct
    res$enrichment      <- sc$enr
    res$log2.enrichment <- if (!is.na(sc$enr) && sc$enr > 0) log2(sc$enr)
                           else NA_real_
    res$pvalue          <- pval
    res
  }

  if (method == "binomial") {
    sc  <- proximity_scan(hits.m, anchors, U, Utot, dthr, count, orientation)
    res <- fill_from_scan(res, sc, min(1, sc$p * max(sc$n.tested, 1L)))

  } else if (method == "ks") {
    dists <- abs(res$dist.to.nearest)
    dists <- dists[is.finite(dists)]
    if (length(dists)) {
      ks <- suppressWarnings(
        stats::ks.test(dists, Ffun, alternative = "greater"))
      res$ks.D   <- unname(ks$statistic)
      res$pvalue <- ks$p.value
    }

  } else { # permutation -- uses the global RNG (set.seed() upstream)
    sc.obs <- proximity_scan(hits.m, anchors, U, Utot, dthr, count, orientation)
    Tperm  <- vapply(seq_len(n.perm), function(b) {
      ra <- sample_anchors_in_universe(anchors, U, perm.per.chrom)
      proximity_scan(hits.m, ra, U, Utot, dthr, count, orientation)$z
    }, numeric(1))
    emp <- (1 + sum(Tperm >= sc.obs$z, na.rm = TRUE)) / (1 + n.perm)
    res <- fill_from_scan(res, sc.obs, emp)
  }

  res
}

## Bind per-motif results, BH-adjust, filter by qvalue, sort.
assemble_proximity_result <- function(motif_ids, per, method, count,
                                      orientation, qvalue) {
  if (!length(motif_ids)) return(empty_proximity_result())
  num <- function(f) vapply(per, function(x) as.numeric(x[[f]]), numeric(1))
  int <- function(f) vapply(per, function(x) as.integer(x[[f]]), integer(1))

  out <- data.frame(
    motif           = as.character(motif_ids),
    motif.i         = seq_along(motif_ids),
    method          = method,
    count           = count,
    orientation     = orientation,
    nhits           = int("nhits"),
    best.dist       = int("best.dist"),
    hits.near       = int("hits.near"),
    hits.far        = int("hits.far"),
    frac.near       = num("frac.near"),
    expected.near   = num("expected.near"),
    enrichment      = num("enrichment"),
    log2.enrichment = num("log2.enrichment"),
    pvalue          = num("pvalue"),
    stringsAsFactors = FALSE
  )
  out$qvalue          <- stats::p.adjust(out$pvalue, method = "BH")
  out$ks.D            <- num("ks.D")
  out$dist.to.nearest <- lapply(per, function(x) x$dist.to.nearest)

  out <- out[, c("motif", "motif.i", "method", "count", "orientation",
                 "nhits", "best.dist", "hits.near", "hits.far", "frac.near",
                 "expected.near", "enrichment", "log2.enrichment", "pvalue",
                 "qvalue", "ks.D", "dist.to.nearest")]
  out <- out[out$qvalue <= qvalue, , drop = FALSE]
  out <- out[order(out$pvalue, out$motif.i), , drop = FALSE]
  rownames(out) <- NULL
  out
}
