#' Sample composition-matched background sequences from a universe.
#'
#' For each sequence in `sequences`, `match_bkg()` draws one or more
#' background sequences from `universe` that share its GC fraction and
#' length. This is the HOMER-style binned matching used by motif
#' enrichment tools to control for sequence-composition biases (GC,
#' length) that would otherwise inflate enrichment estimates relative
#' to shuffle-only nulls.
#'
#' Algorithm (HOMER-style):
#'   1. Compute per-sequence GC fraction and length for `sequences`
#'      and `universe`.
#'   2. Bin the universe into a 2-D grid: GC bins equal-width over
#'      [0, 1], length bins quantile-based on the universe widths.
#'   3. For each target, locate its (GC, length) bin and sample
#'      `n.per.target` universe sequences from that bin (without
#'      replacement when `unique = TRUE`). If the bin is empty or
#'      undersized, expand outward in Manhattan-distance rings.
#'
#' Use the result as a `bkg.sequences` argument to [enrich_motifs2()]
#' or [motif_finder()] for a composition-controlled null.
#' [shuffle_sequences()] is a complementary alternative that
#' randomises each input in place (preserves k-let composition);
#' `match_bkg()` instead samples real sequences from a larger pool.
#'
#' @param sequences `XStringSet`. Target (positive) sequences. DNA or
#'   RNA only -- GC matching is undefined for other alphabets.
#' @param universe `XStringSet`. Pool of candidate background
#'   sequences, typically extracted from a `BSgenome` over a larger
#'   set of genomic regions (e.g. peak-flanking regions, random
#'   genomic tiles, all promoters). Must share `sequences`'s alphabet
#'   and be at least `n.per.target * length(sequences)` long when
#'   `unique = TRUE`.
#' @param n.per.target `integer(1)`. Number of background sequences
#'   to draw per target. Default `1L`.
#' @param n.bins.gc `integer(1)`. Number of equal-width GC bins.
#'   Default `20L`.
#' @param n.bins.length `integer(1)`. Number of quantile-based length
#'   bins. Default `10L`.
#' @param unique `logical(1)`. If `TRUE` (default), each universe
#'   sequence is drawn at most once. Falls back to with-replacement
#'   sampling (with a warning) if the universe is too small to honour
#'   uniqueness.
#' @param rng.seed `integer(1)`. RNG seed for reproducible sampling.
#' @param return.indices `logical(1)`. If `TRUE`, return a
#'   `data.frame` of per-pair match indices and composition values
#'   instead of the matched `XStringSet`. Default `FALSE`.
#'
#' @return If `return.indices = FALSE` (default): an `XStringSet`
#'   (same subclass as `sequences`) of length
#'   `length(sequences) * n.per.target`. The first `n.per.target`
#'   entries correspond to `sequences[[1]]`, the next to
#'   `sequences[[2]]`, etc. If `return.indices = TRUE`: a
#'   `data.frame` with columns `target.i`, `target.gc`,
#'   `target.length`, `universe.i`, `universe.gc`, `universe.length`.
#'
#' @references
#'
#' Heinz S, Benner C, Spann N, et al. (2010). "Simple combinations of
#' lineage-determining transcription factors prime cis-regulatory
#' elements required for macrophage and B cell identities."
#' *Molecular Cell*, **38**(4):576-589.
#' \doi{10.1016/j.molcel.2010.05.004}.
#'
#' Schep AN, Wu B, Buenrostro JD, Greenleaf WJ (2017). "chromVAR:
#' inferring transcription-factor-associated accessibility from
#' single-cell epigenomic data." *Nature Methods*, **14**:975-978.
#' (Algorithmic neighbour using kNN matching; not ported.)
#'
#' @examples
#' \dontrun{
#' library(universalmotif)
#' data(ArabidopsisPromoters)
#' target   <- ArabidopsisPromoters[1:10]
#' universe <- ArabidopsisPromoters[11:50]
#' bkg <- match_bkg(target, universe, n.per.target = 2)
#' ## Use as a null background for enrichment / discovery
#' # enrich_motifs2(motifs, target, bkg.sequences = bkg)
#' # motif_finder(target, bkg.sequences = bkg)
#' plot_match_bkg(target, bkg)
#' }
#'
#' @seealso [shuffle_sequences()], [enrich_motifs2()], [motif_finder()],
#'   [plot_match_bkg()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
match_bkg <- function(sequences, universe,
                     n.per.target   = 1L,
                     n.bins.gc      = 20L,
                     n.bins.length  = 10L,
                     unique         = TRUE,
                     rng.seed       = sample.int(1e6, 1),
                     return.indices = FALSE) {

  ## --- arg validation --------------------------------------------------
  if (missing(sequences) || missing(universe))
    stop("`sequences` and `universe` are both required", call. = FALSE)
  if (!is(sequences, "XStringSet") || !is(universe, "XStringSet"))
    stop("`sequences` and `universe` must both be XStringSet objects",
         call. = FALSE)
  if (!is.numeric(n.per.target) || length(n.per.target) != 1L ||
      n.per.target < 1L)
    stop("`n.per.target` must be a positive integer", call. = FALSE)
  if (!is.numeric(n.bins.gc) || length(n.bins.gc) != 1L || n.bins.gc < 2L)
    stop("`n.bins.gc` must be >= 2", call. = FALSE)
  if (!is.numeric(n.bins.length) || length(n.bins.length) != 1L ||
      n.bins.length < 2L)
    stop("`n.bins.length` must be >= 2", call. = FALSE)
  if (!isTRUEorFALSE(unique))
    stop("`unique` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(return.indices))
    stop("`return.indices` must be a single logical", call. = FALSE)
  if (!is.numeric(rng.seed) || length(rng.seed) != 1L || is.na(rng.seed))
    stop("`rng.seed` must be a single numeric", call. = FALSE)

  ## --- alphabet checks --------------------------------------------------
  alph.t <- seqtype(sequences)
  alph.u <- seqtype(universe)
  if (alph.t != alph.u)
    stop("`sequences` alphabet (", alph.t, ") and `universe` alphabet (",
         alph.u, ") do not match", call. = FALSE)
  if (!alph.t %in% c("DNA", "RNA"))
    stop("`match_bkg()` only supports DNA/RNA sequences; got `",
         alph.t, "`. GC matching is undefined for other alphabets.",
         call. = FALSE)
  if (any(width(universe) == 0L))
    stop("`universe` must not contain zero-length sequences", call. = FALSE)
  if (length(sequences) == 0L)
    stop("`sequences` is empty", call. = FALSE)

  needed <- as.integer(n.per.target) * length(sequences)
  if (unique && length(universe) < needed)
    warning("`universe` (n=", length(universe), ") is smaller than ",
            "n.per.target * length(sequences) (", needed,
            "); falling back to with-replacement sampling.",
            call. = FALSE)

  ## --- composition vectors ---------------------------------------------
  gc.t  <- compute_gc(sequences)
  gc.u  <- compute_gc(universe)
  len.t <- as.integer(width(sequences))
  len.u <- as.integer(width(universe))

  ## --- bin breaks -------------------------------------------------------
  gc.breaks <- seq(0, 1, length.out = as.integer(n.bins.gc) + 1L)
  ## Quantile-based length breaks; clamp uniques to avoid degenerate
  ## breaks when many lengths tie. If the universe is constant-length,
  ## collapse to a single length bin (no length matching to do).
  len.breaks <- unique(stats::quantile(len.u,
                                       probs = seq(0, 1,
                                                   length.out =
                                                     as.integer(n.bins.length) + 1L),
                                       names = FALSE))
  lo <- min(len.t, len.u)
  hi <- max(len.t, len.u)
  if (length(len.breaks) < 2L) {
    len.breaks <- c(lo - 1L, hi + 1L)
  } else {
    len.breaks[1]                  <- min(len.breaks[1], lo) - 1L
    len.breaks[length(len.breaks)] <- max(len.breaks[length(len.breaks)],
                                          hi) + 1L
  }

  ## --- bin assignments --------------------------------------------------
  gc.bin.u  <- pmin(pmax(findInterval(gc.u,  gc.breaks,  rightmost.closed = TRUE),
                         1L), length(gc.breaks)  - 1L)
  gc.bin.t  <- pmin(pmax(findInterval(gc.t,  gc.breaks,  rightmost.closed = TRUE),
                         1L), length(gc.breaks)  - 1L)
  len.bin.u <- pmin(pmax(findInterval(len.u, len.breaks, rightmost.closed = TRUE),
                         1L), length(len.breaks) - 1L)
  len.bin.t <- pmin(pmax(findInterval(len.t, len.breaks, rightmost.closed = TRUE),
                         1L), length(len.breaks) - 1L)

  ## Build the (GC_bin, length_bin) -> vector-of-universe-indices map.
  bin_idx_u <- (gc.bin.u - 1L) * (length(len.breaks) - 1L) + len.bin.u
  cells <- split(seq_along(gc.u), bin_idx_u)

  n.gc.bins  <- length(gc.breaks)  - 1L
  n.len.bins <- length(len.breaks) - 1L

  ## --- per-target matching ---------------------------------------------
  set.seed(as.integer(rng.seed))
  n.targets <- length(sequences)
  matched_u <- integer(n.targets * n.per.target)
  used_global <- if (unique) integer(0) else NULL
  fallback_used <- 0L

  for (ti in seq_len(n.targets)) {
    cand <- find_candidates(gc.bin.t[ti], len.bin.t[ti],
                            cells, n.gc.bins, n.len.bins,
                            used_global, n.per.target, unique)
    if (cand$rings > 1L) fallback_used <- fallback_used + 1L
    drawn <- cand$picked
    matched_u[((ti - 1L) * n.per.target + 1L):(ti * n.per.target)] <- drawn
    if (unique) used_global <- c(used_global, drawn)
  }

  if (fallback_used > 0L)
    warning(fallback_used, " target(s) required ring-expansion fallback ",
            "(empty or undersized (GC, length) bin). Consider widening ",
            "n.bins.gc / n.bins.length or providing a bigger universe.",
            call. = FALSE)

  ## --- assemble output --------------------------------------------------
  if (return.indices) {
    target.i <- rep(seq_len(n.targets), each = n.per.target)
    return(data.frame(
      target.i        = target.i,
      target.gc       = gc.t[target.i],
      target.length   = len.t[target.i],
      universe.i      = matched_u,
      universe.gc     = gc.u[matched_u],
      universe.length = len.u[matched_u],
      stringsAsFactors = FALSE
    ))
  }

  out <- universe[matched_u]
  if (is.null(names(out)) || any(!nzchar(names(out)))) {
    target.i <- rep(seq_len(n.targets), each = n.per.target)
    k        <- rep(seq_len(n.per.target), times = n.targets)
    names(out) <- sprintf("bkg_%d_%d", target.i, k)
  }
  out
}

#' Visualise composition match between target and background sequences.
#'
#' Companion to [match_bkg()]. Overlays GC-fraction and sequence-length
#' density curves for the target and matched-background sets in a
#' two-panel `ggplot`, useful for visually verifying that the binned
#' matching landed where it was supposed to.
#'
#' @param sequences `XStringSet`. The target sequences passed to
#'   [match_bkg()].
#' @param bkg `XStringSet`. The matched background sequences returned
#'   by [match_bkg()].
#' @param by `character`. Composition axes to plot. Subset of
#'   `c("gc", "length")`. Default both.
#' @param bins `integer(1)`. Ignored when geom = "density"; kept for
#'   API parity with other v2 plotting helpers. Default `30L`.
#'
#' @return A `ggplot` object.
#'
#' @seealso [match_bkg()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_match_bkg <- function(sequences, bkg,
                           by = c("gc", "length"),
                           bins = 30L) {
  if (!is(sequences, "XStringSet") || !is(bkg, "XStringSet"))
    stop("both `sequences` and `bkg` must be XStringSet objects",
         call. = FALSE)
  by <- match.arg(by, c("gc", "length"), several.ok = TRUE)

  long <- rbind(
    data.frame(set = "target", gc = compute_gc(sequences),
               length = as.numeric(width(sequences)),
               stringsAsFactors = FALSE),
    data.frame(set = "matched", gc = compute_gc(bkg),
               length = as.numeric(width(bkg)),
               stringsAsFactors = FALSE)
  )

  ## Pivot to a long form keyed by panel ("gc" or "length").
  pieces <- list()
  if ("gc" %in% by) {
    pieces[["gc"]] <- data.frame(set = long$set, axis = "GC fraction",
                                 value = long$gc,
                                 stringsAsFactors = FALSE)
  }
  if ("length" %in% by) {
    pieces[["length"]] <- data.frame(set = long$set, axis = "Sequence length",
                                     value = long$length,
                                     stringsAsFactors = FALSE)
  }
  long2 <- do.call(rbind, pieces)
  long2$set <- factor(long2$set, levels = c("target", "matched"))

  ggplot2::ggplot(long2, ggplot2::aes(x = .data$value, colour = .data$set,
                                       fill = .data$set)) +
    ggplot2::geom_density(alpha = 0.25, linewidth = 0.6) +
    ggplot2::facet_wrap(~ axis, scales = "free", ncol = length(by)) +
    ggplot2::scale_colour_manual(values = c(target = "#1f77b4",
                                            matched = "#ff7f0e")) +
    ggplot2::scale_fill_manual(values   = c(target = "#1f77b4",
                                            matched = "#ff7f0e")) +
    ggplot2::labs(x = NULL, y = "Density", colour = NULL, fill = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background   = ggplot2::element_rect(fill = NA, colour = NA),
      panel.grid         = ggplot2::element_blank(),
      panel.border       = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(colour = "black"),
      axis.line.y.left   = ggplot2::element_line(colour = "black"),
      legend.position    = "top"
    )
}

## ---------------------------------------------------------------------------
## Internal helpers
## ---------------------------------------------------------------------------

## GC fraction per sequence; works for DNA and RNA (G + C count / length).
compute_gc <- function(seqs) {
  ## letterFrequency with "GC" returns counts of G and C (combined).
  gc <- as.numeric(Biostrings::letterFrequency(seqs, "GC", as.prob = FALSE))
  w  <- as.numeric(width(seqs))
  out <- gc / pmax(w, 1)
  out[!is.finite(out)] <- 0
  out
}

## Ring-expanding candidate finder. Returns picked indices + the ring
## level used (0 = home cell, 1 = first ring, ...).
find_candidates <- function(gc.b, len.b, cells, n.gc.bins, n.len.bins,
                            used_global, n.needed, unique) {
  ring <- 0L
  picked <- integer(0)
  while (length(picked) < n.needed) {
    cand <- ring_candidates(gc.b, len.b, ring, cells, n.gc.bins, n.len.bins)
    if (unique && length(used_global) > 0L)
      cand <- setdiff(cand, used_global)
    cand <- setdiff(cand, picked)
    if (length(cand) > 0L) {
      take <- min(length(cand), n.needed - length(picked))
      picked <- c(picked, sample(cand, take, replace = FALSE))
    }
    ring <- ring + 1L
    max_ring <- max(n.gc.bins, n.len.bins)
    if (ring > max_ring) {
      ## Hit the global ceiling: fall back to with-replacement sampling
      ## from the entire universe to fill the remainder.
      need <- n.needed - length(picked)
      if (need > 0L) {
        all_u <- unlist(cells, use.names = FALSE)
        picked <- c(picked, sample(all_u, need, replace = TRUE))
      }
      break
    }
  }
  list(picked = picked, rings = ring)
}

## Universe-index candidates at Manhattan distance `ring` from (gc.b, len.b).
## ring == 0 returns the home cell.
ring_candidates <- function(gc.b, len.b, ring, cells, n.gc.bins, n.len.bins) {
  if (ring == 0L) {
    key <- (gc.b - 1L) * n.len.bins + len.b
    return(cells[[as.character(key)]] %||% integer(0))
  }
  out <- integer(0)
  for (dg in -ring:ring) {
    dl <- ring - abs(dg)
    for (dlen in c(dl, -dl)) {
      g2 <- gc.b  + dg
      l2 <- len.b + dlen
      if (g2 < 1L || g2 > n.gc.bins) next
      if (l2 < 1L || l2 > n.len.bins) next
      key <- (g2 - 1L) * n.len.bins + l2
      cell <- cells[[as.character(key)]]
      if (!is.null(cell)) out <- c(out, cell)
      if (dlen == 0L) break  # avoid double-counting the (g2, len.b) cell
    }
  }
  unique(out)
}
