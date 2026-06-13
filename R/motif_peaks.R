#' CentriMo-style positional enrichment of motif hits.
#'
#' `motif_peaks()` tests whether motif hits cluster non-uniformly along
#' input sequences. It implements an analytical CentriMo-style test
#' (Bailey & Machanick 2012): under the null that hit positions are
#' uniformly distributed over `[1, seq.length]`, the count of hits
#' falling inside a candidate window follows a binomial distribution.
#' For each motif and each candidate window, a one-sided binomial
#' p-value is computed; the best-scoring window per motif is reported,
#' Bonferroni-corrected over the number of windows tested, then BH-
#' adjusted across motifs.
#'
#' Two modes are supported. `mode = "central"` (default) tests only
#' windows centred on `seq.length / 2`, which is the appropriate
#' hypothesis for centred input sequences (e.g. ChIP-seq peak summits).
#' `mode = "local"` additionally varies the window centre, scanning the
#' whole sequence for any positionally-enriched region; this is more
#' permissive and uses a larger Bonferroni correction.
#'
#' All sequences are assumed to be of equal length; the function does
#' not enforce centring on a biological reference, but its results are
#' only meaningful if such centring is in place upstream.
#'
#' @param hits Motif hit table from [scan_sequences()] or
#'   [scan_sequences2()]. Accepted as either a `data.frame` or a
#'   `GRanges`. Required columns / metadata: `motif`, `score`,
#'   `start`, and either `end` (scan_sequences2) or `stop`
#'   (scan_sequences). For `data.frame` input, a `sequence` /
#'   `sequence.i` column identifies which sequence each hit belongs
#'   to. For `GRanges` input, the sequence is taken from
#'   `seqnames()`.
#' @param seq.length `integer(1)`. Common length of the input
#'   sequences (in bases). Required for `data.frame` input. For
#'   `GRanges` input, taken from `seqlengths(hits)` if not provided.
#' @param mode `character(1)`. `"central"` (default) tests only
#'   centre-of-sequence windows; `"local"` also varies the window
#'   centre across the sequence.
#' @param seq.strand `NULL` or a named `character` vector mapping each
#'   sequence id (the `sequence` / `seqnames()` values in `hits`) to
#'   `"+"` or `"-"`. When supplied, hit positions on `"-"` sequences are
#'   reflected into the feature 5'->3' frame (`seq.length + 1 - center`),
#'   so that off-centre signals at the same feature-relative position
#'   reinforce across strands instead of splitting into mirror-image
#'   peaks. This only changes `mode = "local"` results; `mode =
#'   "central"` is symmetric about the centre and so is unaffected (the
#'   reflection is still applied, so the reported positions read in
#'   feature frame). Unnecessary when sequences were extracted with
#'   [Biostrings::getSeq()] from stranded ranges, since those are already
#'   oriented. Default `NULL` (strand-agnostic, the original behaviour).
#' @param qvalue `numeric(1)`. BH-adjusted q-value cutoff for the
#'   final reported motifs. Default `0.1`. Set to `1` to return every
#'   motif.
#' @param min.window,max.window `integer(1)`. Bounds of the window
#'   widths to test. `max.window = NULL` (default) uses
#'   `floor(seq.length / 2)`. `min.window` defaults to `10L`.
#' @param window.step `integer(1)`. Step between consecutive window
#'   widths. Default `10L`.
#' @param position.step `integer(1)`. Step between consecutive window
#'   centres in `mode = "local"`. Ignored in `mode = "central"`.
#'   Default `5L`.
#' @param nthreads `numeric(1)`. Number of threads. `nthreads = 0`
#'   uses all available cores. Default `1`.
#'
#' @return A `data.frame` with one row per motif passing the q-value
#'   cutoff, sorted by `pvalue` ascending. Columns: `motif`,
#'   `motif.i`, `mode`, `seq.length`, `nhits`, `best.window`,
#'   `best.center`, `hits.in`, `hits.out`, `expected.in`,
#'   `enrichment`, `log2.enrichment`, `pvalue`, `qvalue`, `centers`
#'   (list-column of integer vectors of hit centres, one entry per
#'   sequence after best-hit-per-sequence dedup; consumed by
#'   [plot_motif_peaks()]).
#'
#' @references
#'
#' Bailey TL, Machanick P (2012). "Inferring direct DNA binding from
#' ChIP-seq." *Nucleic Acids Research*, **40**(17):e128.
#' \doi{10.1093/nar/gks433}.
#'
#' @examples
#' ## Pre-centred sequences (e.g. ChIP-seq peaks centred on their summit):
#' library(Biostrings)
#' set.seed(1)
#' seqs <- create_sequences(seqnum = 200, seqlen = 500, rng.seed = 1)
#' planted <- DNAString("TTGACATA")
#' for (i in seq_len(150)) {
#'   pos <- 246L + sample.int(8, 1) - 1L
#'   subseq(seqs[[i]], start = pos, width = 8) <- planted
#' }
#' m <- create_motif("TTGACATA", name = "test")
#' ## scan_sequences() returns a table that motif_peaks() reads directly;
#' ## pass seq.length since the data.frame does not carry it.
#' hits <- scan_sequences(m, seqs, threshold = 1e-3, threshold.type = "pvalue")
#' motif_peaks(as.data.frame(hits), seq.length = 500)
#'
#' ## Enrichment near a set of genomic coordinates. Given point
#' ## coordinates in a genome (peak summits, TSSs, SNPs, and so on),
#' ## expand each to a fixed-width window centred on the coordinate and
#' ## scan those windows. "Enriched near the anchor" is then exactly the
#' ## central hypothesis motif_peaks() already tests, so nothing beyond
#' ## extracting the windows is required.
#' \dontrun{
#' library(universalmotif)
#' library(GenomicRanges)
#' library(BSgenome.Athaliana.TAIR.TAIR9)
#'
#' ## A set of anchor coordinates (here, arbitrary points on Chr1).
#' anchors <- GRanges("Chr1",
#'                    IRanges(start = seq(1e5, 5e5, by = 1e4), width = 1))
#'
#' ## Expand to fixed 500 bp windows centred on each anchor. The
#' ## equal-width filter matters because motif_peaks() compares positions
#' ## across sequences, so a window clipped at a chromosome end would
#' ## misplace its centre. For stranded anchors (e.g. TSSs), give
#' ## `anchors` a strand and getSeq() will orient each window
#' ## consistently up- or downstream of the anchor.
#' windows <- resize(anchors, width = 500, fix = "center")
#' windows <- trim(windows)
#' windows <- windows[width(windows) == 500]
#'
#' seqs <- getSeq(Athaliana, windows)
#' names(seqs) <- paste0("anchor_", seq_along(seqs))
#'
#' m    <- create_motif("TTGACATA", name = "example")
#' hits <- scan_sequences2(m, seqs, pvalue = 1e-3, return.granges = TRUE)
#' ## seq.length is inferred from seqlengths() here, but pass it
#' ## explicitly when the hit table doesn't carry equal seqlengths.
#' motif_peaks(hits, seq.length = 500)
#' }
#'
#' @seealso [scan_sequences2()], [plot_motif_peaks()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
motif_peaks <- function(hits, seq.length = NULL,
                        mode = c("central", "local"),
                        seq.strand = NULL,
                        qvalue = 0.1,
                        min.window = 10L, max.window = NULL,
                        window.step = 10L, position.step = 5L,
                        nthreads = 1) {

  ## --- arg validation --------------------------------------------------
  if (missing(hits))
    stop("`hits` is required", call. = FALSE)
  mode <- match.arg(mode)
  if (!is.numeric(qvalue) || length(qvalue) != 1L || is.na(qvalue) ||
      qvalue <= 0 || qvalue > 1)
    stop("`qvalue` must be a single numeric in (0, 1]", call. = FALSE)
  if (!is.numeric(min.window) || length(min.window) != 1L || min.window < 1L)
    stop("`min.window` must be a positive integer", call. = FALSE)
  if (!is.numeric(window.step) || length(window.step) != 1L ||
      window.step < 1L)
    stop("`window.step` must be a positive integer", call. = FALSE)
  if (!is.numeric(position.step) || length(position.step) != 1L ||
      position.step < 1L)
    stop("`position.step` must be a positive integer", call. = FALSE)

  nthreads <- resolve_nthreads(nthreads)

  ## --- normalise hits to a flat data.frame ----------------------------
  norm <- normalize_motif_peaks_input(hits, seq.length)
  hits.df    <- norm$hits         # cols: motif, sequence, center, score
  seq.length <- norm$seq.length

  ## Orient '-'-strand sequences into the feature 5'->3' frame. Only
  ## affects mode = "local" (central is symmetric), but applied either way
  ## so reported positions read in feature frame.
  if (!is.null(seq.strand))
    hits.df$center <- reflect_neg_strand_centers(hits.df$sequence,
                                                 hits.df$center,
                                                 seq.length, seq.strand)

  if (is.null(max.window))
    max.window <- floor(seq.length / 2)
  if (max.window > seq.length)
    max.window <- seq.length
  if (max.window < min.window)
    stop("`max.window` (", max.window, ") < `min.window` (", min.window, ")",
         call. = FALSE)

  window_widths <- seq.int(min.window, max.window, by = window.step)
  if (length(window_widths) == 0L)
    stop("no window widths fit between min.window and max.window with step ",
         window.step, call. = FALSE)

  ## --- best-hit per (motif, sequence) ----------------------------------
  ## Keep only the highest-score hit per sequence per motif (CentriMo
  ## convention) so that motifs whose hits cluster within a single
  ## sequence are not double-counted.
  hits.df <- best_hit_per_seq(hits.df)
  if (nrow(hits.df) == 0L) return(empty_peaks_result())

  ## --- per-motif positional test --------------------------------------
  motif_ids <- sort(unique(hits.df$motif))
  per_motif <- if (nthreads > 1L && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(motif_ids,
                       function(m) test_one_motif(hits.df[hits.df$motif == m, ],
                                                  seq.length, window_widths,
                                                  mode, position.step),
                       mc.cores = nthreads)
  } else {
    lapply(motif_ids, function(m) {
      test_one_motif(hits.df[hits.df$motif == m, ], seq.length,
                     window_widths, mode, position.step)
    })
  }

  ## --- assemble output -------------------------------------------------
  out <- data.frame(
    motif         = motif_ids,
    motif.i       = seq_along(motif_ids),
    mode          = mode,
    seq.length    = as.integer(seq.length),
    nhits         = vapply(per_motif, `[[`, integer(1),  "nhits"),
    best.window   = vapply(per_motif, `[[`, integer(1),  "best.window"),
    best.center   = vapply(per_motif, `[[`, numeric(1),  "best.center"),
    hits.in       = vapply(per_motif, `[[`, integer(1),  "hits.in"),
    hits.out      = vapply(per_motif, `[[`, integer(1),  "hits.out"),
    expected.in   = vapply(per_motif, `[[`, numeric(1),  "expected.in"),
    enrichment    = vapply(per_motif, `[[`, numeric(1),  "enrichment"),
    log2.enrichment = vapply(per_motif, `[[`, numeric(1), "log2.enrichment"),
    pvalue        = vapply(per_motif, `[[`, numeric(1),  "pvalue"),
    stringsAsFactors = FALSE
  )
  out$qvalue  <- stats::p.adjust(out$pvalue, method = "BH")
  out$centers <- lapply(per_motif, `[[`, "centers")

  ## filter + sort
  out <- out[out$qvalue <= qvalue, , drop = FALSE]
  out <- out[order(out$pvalue, out$motif.i), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Plot the position distribution of motif hits.
#'
#' Companion to [motif_peaks()]. Returns a `ggplot` faceted by motif,
#' showing a histogram of best-hit centre positions and a shaded
#' rectangle marking the most-enriched window.
#'
#' @param peaks `data.frame` returned by [motif_peaks()]. Must carry
#'   the `centers` list-column.
#' @param motifs `character` or `NULL`. Subset to these motifs (by
#'   name). `NULL` (default) plots every row of `peaks`.
#' @param ncol `integer(1)`. Number of facet columns. Default `2`.
#' @param bins `integer(1)`. Histogram bin count. Default `50`.
#' @param fill.window `character(1)`. Colour of the best-window shade.
#'   Default `"#ffe4a3"` (a pale yellow).
#'
#' @return A `ggplot` object.
#'
#' @examples
#' library(Biostrings)
#' set.seed(1)
#' seqs <- create_sequences(seqnum = 200, seqlen = 500, rng.seed = 1)
#' planted <- DNAString("TTGACATA")
#' for (i in seq_len(150)) {
#'   pos <- 246L + sample.int(8, 1) - 1L
#'   subseq(seqs[[i]], start = pos, width = 8) <- planted
#' }
#' m <- create_motif("TTGACATA", name = "test")
#' hits <- scan_sequences(m, seqs, threshold = 1e-3, threshold.type = "pvalue")
#' peaks <- motif_peaks(as.data.frame(hits), seq.length = 500, qvalue = 1)
#' plot_motif_peaks(peaks)
#'
#' @seealso [motif_peaks()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_motif_peaks <- function(peaks, motifs = NULL, ncol = 2L,
                             bins = 50L, fill.window = "#ffe4a3") {

  if (!is.data.frame(peaks))
    stop("`peaks` must be a data.frame from motif_peaks()", call. = FALSE)
  if (!"centers" %in% names(peaks))
    stop("`peaks` is missing the `centers` list-column ",
         "(was it produced by motif_peaks()?)", call. = FALSE)
  if (nrow(peaks) == 0L)
    stop("`peaks` has 0 rows -- no enriched motifs to plot", call. = FALSE)

  if (!is.null(motifs)) {
    keep <- peaks$motif %in% motifs
    if (!any(keep))
      stop("none of the requested motifs are in `peaks`", call. = FALSE)
    peaks <- peaks[keep, , drop = FALSE]
  }

  ## Flatten centers into a long data.frame for ggplot.
  long <- do.call(rbind, lapply(seq_len(nrow(peaks)), function(i) {
    data.frame(motif = peaks$motif[i], center = peaks$centers[[i]],
               stringsAsFactors = FALSE)
  }))
  long$motif <- factor(long$motif, levels = peaks$motif)

  ## Best-window rectangles, one per motif, for geom_rect().
  rects <- data.frame(
    motif = factor(peaks$motif, levels = peaks$motif),
    xmin  = peaks$best.center - peaks$best.window / 2,
    xmax  = peaks$best.center + peaks$best.window / 2,
    stringsAsFactors = FALSE
  )

  ## Annotate the strip label with the per-motif stats (avoids any
  ## overlap with bars that might land at the corners).
  strip_labels <- setNames(
    sprintf("%s  (p = %.2g | %.1fx)",
            peaks$motif, peaks$pvalue, peaks$enrichment),
    peaks$motif
  )

  ggplot2::ggplot(long, ggplot2::aes(x = .data$center)) +
    ggplot2::geom_rect(data = rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = -Inf, ymax = Inf),
                       fill = fill.window, alpha = 0.5,
                       inherit.aes = FALSE) +
    ggplot2::geom_histogram(bins = bins, colour = NA, fill = "black") +
    ggplot2::facet_wrap(~ motif, ncol = ncol, scales = "free_y",
                        axes = "all_x", axis.labels = "all_x",
                        labeller = ggplot2::as_labeller(strip_labels)) +
    ggplot2::labs(x = "Hit centre position", y = "Count") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background    = ggplot2::element_rect(fill = NA, colour = NA),
      ## no interior grid lines at all
      panel.grid          = ggplot2::element_blank(),
      ## L-shaped axes: kill the full panel border, draw only bottom + left
      panel.border        = ggplot2::element_blank(),
      axis.line.x.bottom  = ggplot2::element_line(colour = "black"),
      axis.line.y.left    = ggplot2::element_line(colour = "black")
    )
}

## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

empty_peaks_result <- function() {
  data.frame(
    motif           = character(0),
    motif.i         = integer(0),
    mode            = character(0),
    seq.length      = integer(0),
    nhits           = integer(0),
    best.window     = integer(0),
    best.center     = numeric(0),
    hits.in         = integer(0),
    hits.out        = integer(0),
    expected.in     = numeric(0),
    enrichment      = numeric(0),
    log2.enrichment = numeric(0),
    pvalue          = numeric(0),
    qvalue          = numeric(0),
    centers         = I(list()),
    stringsAsFactors = FALSE
  )
}

## Normalise either a data.frame (scan_sequences / scan_sequences2) or a
## GRanges into a flat data.frame with cols (motif, sequence, center,
## score). Also resolves seq.length from input when not supplied.
normalize_motif_peaks_input <- function(hits, seq.length) {

  is_gr <- inherits(hits, "GRanges")

  if (is_gr) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE))
      stop("`hits` is a GRanges; the GenomicRanges package must be installed.",
           call. = FALSE)
    if (!requireNamespace("GenomeInfoDb", quietly = TRUE))
      stop("`hits` is a GRanges; the GenomeInfoDb package must be installed.",
           call. = FALSE)
    sl <- GenomeInfoDb::seqlengths(hits)
    if (is.null(seq.length)) {
      if (length(sl) == 0L || any(is.na(sl)))
        stop("could not resolve `seq.length` from seqlengths(hits); ",
             "pass `seq.length` explicitly.", call. = FALSE)
      if (length(unique(sl)) != 1L)
        stop("seqlengths(hits) are not all equal; ",
             "motif_peaks() requires equal-length input sequences. ",
             "Pass `seq.length` explicitly if you've trimmed.",
             call. = FALSE)
      seq.length <- as.integer(sl[1])
    }
    mc <- S4Vectors::mcols(hits)
    if (!"motif" %in% names(mc))
      stop("GRanges input must have a `motif` metadata column", call. = FALSE)
    if (!"score" %in% names(mc))
      stop("GRanges input must have a `score` metadata column", call. = FALSE)
    df <- data.frame(
      motif    = as.character(mc$motif),
      sequence = as.character(GenomicRanges::seqnames(hits)),
      center   = (GenomicRanges::start(hits) + GenomicRanges::end(hits)) / 2,
      score    = as.numeric(mc$score),
      stringsAsFactors = FALSE
    )
  } else if (is.data.frame(hits)) {
    if (is.null(seq.length))
      stop("`seq.length` is required when `hits` is a data.frame.",
           call. = FALSE)
    if (!"motif" %in% names(hits))
      stop("data.frame `hits` must have a `motif` column", call. = FALSE)
    if (!"score" %in% names(hits))
      stop("data.frame `hits` must have a `score` column", call. = FALSE)
    if (!"start" %in% names(hits))
      stop("data.frame `hits` must have a `start` column", call. = FALSE)
    end_col <- if ("end" %in% names(hits)) hits$end          # scan_sequences2
               else if ("stop" %in% names(hits)) hits$stop    # scan_sequences
               else stop("data.frame `hits` must have an `end` or `stop` column",
                         call. = FALSE)
    seq_col <- if ("sequence" %in% names(hits)) hits$sequence
               else if ("sequence.i" %in% names(hits)) as.character(hits$sequence.i)
               else stop("data.frame `hits` must have a `sequence` ",
                         "or `sequence.i` column", call. = FALSE)
    df <- data.frame(
      motif    = as.character(hits$motif),
      sequence = as.character(seq_col),
      center   = (hits$start + end_col) / 2,
      score    = as.numeric(hits$score),
      stringsAsFactors = FALSE
    )
  } else {
    stop("`hits` must be a data.frame or a GRanges", call. = FALSE)
  }

  if (!is.numeric(seq.length) || length(seq.length) != 1L ||
      seq.length < 2L)
    stop("`seq.length` must be a single integer >= 2", call. = FALSE)

  list(hits = df, seq.length = as.integer(seq.length))
}

## Reflect hit centres on '-'-strand sequences into the feature 5'->3'
## frame so up/downstream is consistent across strands. `seq.strand` is a
## named vector (sequence id -> "+"/"-") covering every sequence in `df`.
reflect_neg_strand_centers <- function(sequence, center, seq.length,
                                       seq.strand) {
  if (is.null(names(seq.strand)))
    stop("`seq.strand` must be a named vector mapping sequence ids to ",
         "\"+\"/\"-\"", call. = FALSE)
  sv <- as.character(seq.strand)
  names(sv) <- names(seq.strand)
  miss <- setdiff(unique(sequence), names(sv))
  if (length(miss))
    stop("`seq.strand` is missing entries for sequence(s): ",
         paste(utils::head(miss, 5L), collapse = ", "),
         if (length(miss) > 5L) ", ..." else "", call. = FALSE)
  s <- sv[sequence]
  if (!all(s %in% c("+", "-")))
    stop("`seq.strand` values must be \"+\" or \"-\"", call. = FALSE)
  neg <- s == "-"
  center[neg] <- seq.length + 1 - center[neg]
  center
}

## Keep only the highest-score hit per (motif, sequence) pair.
best_hit_per_seq <- function(df) {
  if (nrow(df) == 0L) return(df)
  key <- paste0(df$motif, "\t", df$sequence)
  keep <- as.logical(
    stats::ave(df$score, key,
               FUN = function(x) seq_along(x) == which.max(x))
  )
  df[keep, , drop = FALSE]
}

## Per-motif positional test. `m_hits` is a slice of the normalised
## data.frame containing only hits for one motif.
test_one_motif <- function(m_hits, seq.length, window_widths, mode,
                           position.step) {
  centers <- m_hits$center
  N       <- length(centers)
  seq_mid <- seq.length / 2

  if (N == 0L) {
    return(list(nhits = 0L, best.window = NA_integer_,
                best.center = NA_real_,
                hits.in = 0L, hits.out = 0L,
                expected.in = NA_real_, enrichment = NA_real_,
                log2.enrichment = NA_real_, pvalue = 1, centers = integer(0)))
  }

  best_p <- 1; best_w <- NA_integer_; best_c <- NA_real_
  best_k <- 0L
  n_tests <- 0L

  if (mode == "central") {
    for (w in window_widths) {
      half <- w / 2
      k    <- sum(abs(centers - seq_mid) <= half)
      p_w  <- w / seq.length
      pv   <- stats::pbinom(k - 1L, N, p_w, lower.tail = FALSE)
      n_tests <- n_tests + 1L
      if (pv < best_p) {
        best_p <- pv; best_w <- w; best_c <- seq_mid; best_k <- k
      }
    }
  } else {
    ## local: vary window centre too
    for (w in window_widths) {
      half <- w / 2
      lo   <- ceiling(half)
      hi   <- floor(seq.length - half)
      if (hi < lo) next
      cs   <- seq.int(lo, hi, by = position.step)
      p_w  <- w / seq.length
      for (cc in cs) {
        k  <- sum(abs(centers - cc) <= half)
        pv <- stats::pbinom(k - 1L, N, p_w, lower.tail = FALSE)
        n_tests <- n_tests + 1L
        if (pv < best_p) {
          best_p <- pv; best_w <- w; best_c <- cc; best_k <- k
        }
      }
    }
  }

  ## Bonferroni-correct over tests performed.
  bonf_p <- min(1, best_p * max(n_tests, 1L))
  exp_in <- N * best_w / seq.length
  enr    <- if (is.finite(exp_in) && exp_in > 0) best_k / exp_in else NA_real_

  list(
    nhits           = as.integer(N),
    best.window     = if (is.na(best_w)) NA_integer_ else as.integer(best_w),
    best.center     = as.numeric(best_c),
    hits.in         = as.integer(best_k),
    hits.out        = as.integer(N - best_k),
    expected.in     = as.numeric(exp_in),
    enrichment      = as.numeric(enr),
    log2.enrichment = if (!is.na(enr) && enr > 0) log2(enr) else NA_real_,
    pvalue          = as.numeric(bonf_p),
    centers         = as.integer(round(centers))
  )
}
