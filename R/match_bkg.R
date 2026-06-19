#' Sample composition-matched background sequences from a universe.
#'
#' For each sequence in `sequences`, [match_bkg()] draws one or more
#' background sequences from `universe` that share its GC fraction and
#' length. This is the HOMER-style binned matching used by motif
#' enrichment tools to control for sequence-composition biases (GC,
#' length) that would otherwise inflate enrichment estimates relative
#' to shuffle-only nulls. Matching can be extended to additional
#' external per-sequence covariates (e.g. ChIP signal strength,
#' conservation, mappability) via `covariates`, which simply add more
#' binned axes to the grid.
#'
#' Algorithm (HOMER-style):
#'   1. Compute per-sequence GC fraction and length for `sequences`
#'      and `universe` (plus any supplied `covariates`).
#'   2. Bin the universe into a grid: GC bins equal-width over
#'      \[0, 1\], length bins quantile-based on the universe widths,
#'      and one extra axis per covariate (quantile- or equal-width
#'      binned, see `bin.type.covariates`).
#'   3. For each target, locate its joint bin and sample
#'      `n.per.target` universe sequences from that bin (without
#'      replacement when `unique = TRUE`). If the bin is empty or
#'      undersized, expand outward in Manhattan-distance rings.
#'
#' Use the result as a `bkg.sequences` argument to [enrich_motifs_lite()]
#' or [motif_finder()] for a composition-controlled null.
#' [shuffle_sequences()] is a complementary alternative that
#' randomises each input in place (preserves k-let composition);
#' [match_bkg()] instead samples real sequences from a larger pool.
#'
#' @param sequences `XStringSet`. Target (positive) sequences. DNA or
#'   RNA only -- GC matching is undefined for other alphabets.
#' @param universe `XStringSet` or `NULL`. Pool of candidate
#'   background sequences. Must share `sequences`'s alphabet and be
#'   at least `n.per.target * length(sequences)` long when
#'   `unique = TRUE`. Mutually exclusive with `genome`: supply
#'   exactly one of the two.
#' @param genome A `BSgenome` object, or `NULL`. When supplied, the
#'   universe is built internally by sampling `n.candidates` random
#'   genomic windows whose widths match the target width distribution,
#'   optionally excluding any window that overlaps `exclude`. This is
#'   the typical shortcut for ChIP-seq style backgrounds, where the
#'   universe would otherwise be a manual `sample()` +
#'   `subsetByOverlaps()` + `getSeq()` boilerplate. Requires the
#'   `BSgenome`, `GenomicRanges`, and `IRanges` packages.
#' @param exclude `GRanges` or `NULL`. Only honoured when `genome` is
#'   non-`NULL`. Random windows that overlap any range in `exclude`
#'   are dropped before sampling. Typical use: pass the target peak
#'   set so the background doesn't sample from peak regions.
#' @param n.candidates `integer(1)` or `NULL`. Only honoured when
#'   `genome` is non-`NULL`. Number of random windows to sample
#'   before any `exclude` filtering. Default `NULL` picks
#'   `max(10000, length(sequences) * 10)`.
#' @param chromosomes `character` or `NULL`. Only honoured when
#'   `genome` is non-`NULL`. Names of seqlevels to sample from.
#'   Default `NULL` uses `GenomeInfoDb::standardChromosomes(genome)`
#'   when available (typically the autosomes plus the sex
#'   chromosomes), falling back to all seqlevels of the genome.
#' @param n.per.target `integer(1)`. Number of background sequences
#'   to draw per target. Default `1L`.
#' @param n.bins.gc `integer(1)`. Number of equal-width GC bins.
#'   Default `20L`.
#' @param n.bins.length `integer(1)`. Number of quantile-based length
#'   bins. Default `10L`.
#' @param covariates `NULL`, or a `data.frame` / `matrix` / numeric
#'   vector of external per-sequence covariates to additionally match
#'   on, with one row per sequence in `sequences` and one numeric
#'   column per covariate. A bare numeric vector is treated as a single
#'   covariate named `"covariate1"`. Each covariate becomes an extra
#'   binned matching axis. Column names must not be `"i"`, `"gc"` or
#'   `"length"` (these collide with the `return.indices` columns).
#'   Covariate values must be finite (no `NA`/`NaN`/`Inf`). Mutually
#'   exclusive with `genome`: internally-sampled genomic windows have
#'   no covariate values, so covariate matching requires an explicit
#'   `universe` paired with `universe.covariates`. Default `NULL`.
#' @param universe.covariates `NULL`, or the covariate table for the
#'   `universe` sequences. Required (and only used) when `covariates`
#'   is supplied. Must have one row per `universe` sequence and the
#'   same column names, in the same order, as `covariates`. Default
#'   `NULL`.
#' @param n.bins.covariates `integer(1)` applied to every covariate, or
#'   an integer vector named by covariate column (naming every
#'   covariate). Number of bins per covariate axis; each must be `>= 2`.
#'   Default `10L`.
#' @param bin.type.covariates `character`. How to bin each covariate
#'   axis: `"quantile"` (default, breaks from the universe-covariate
#'   quantiles, like length) or `"equalwidth"` (equal-width breaks over
#'   the combined target/universe range, like GC). A single value
#'   applies to every covariate, or supply a character vector named by
#'   covariate column. Quantile binning is rank-based, hence scale- and
#'   shift-invariant; no covariate normalisation is needed.
#' @param unique `logical(1)`. If `TRUE` (default), each universe
#'   sequence is drawn at most once. Falls back to with-replacement
#'   sampling (with a warning) if the universe is too small to honour
#'   uniqueness.
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
#'   `target.length`, `universe.i`, `universe.gc`, `universe.length`,
#'   followed (when `covariates` is supplied) by a
#'   `target.<name>` / `universe.<name>` pair for each covariate.
#'
#' @references
#'
#' Heinz S, Benner C, Spann N, et al. (2010). "Simple combinations of
#' lineage-determining transcription factors prime cis-regulatory
#' elements required for macrophage and B cell identities."
#' *Molecular Cell*, **38**(4):576-589.
#' \doi{10.1016/j.molcel.2010.05.004}.
#'
#' @examples
#' data(ArabidopsisPromoters)
#' target   <- ArabidopsisPromoters[1:10]
#' universe <- ArabidopsisPromoters[11:50]
#' bkg <- match_bkg(target, universe, n.per.target = 2)
#' ## Use as a null background for enrichment / discovery
#' # enrich_motifs_lite(motifs, target, bkg.sequences = bkg)
#' # motif_finder(target, bkg.sequences = bkg)
#' plot_match_bkg(target, bkg)
#'
#' ## Additionally match on an external covariate (e.g. a ChIP signal
#' ## score aligned to each sequence). Covariate matching needs an
#' ## explicit universe + universe.covariates.
#' sig.t <- runif(length(target))
#' sig.u <- runif(length(universe))
#' idx <- match_bkg(target, universe,
#'                  covariates          = data.frame(signal = sig.t),
#'                  universe.covariates = data.frame(signal = sig.u),
#'                  return.indices      = TRUE)
#' plot_match_bkg(indices = idx)   # GC, length and signal panels
#'
#' \dontrun{
#' ## With a BSgenome, sample the universe from random genomic windows
#' ## whose widths match the targets:
#' library(BSgenome.Athaliana.TAIR.TAIR9)
#' bkg <- match_bkg(target, genome = Athaliana, n.per.target = 2)
#' }
#'
#' @seealso [shuffle_sequences()], [enrich_motifs_lite()], [motif_finder()],
#'   [plot_match_bkg()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
match_bkg <- function(sequences, universe = NULL,
                     genome              = NULL,
                     exclude             = NULL,
                     n.candidates        = NULL,
                     chromosomes         = NULL,
                     n.per.target        = 1L,
                     n.bins.gc           = 20L,
                     n.bins.length       = 10L,
                     covariates          = NULL,
                     universe.covariates = NULL,
                     n.bins.covariates   = 10L,
                     bin.type.covariates = "quantile",
                     unique              = TRUE,
                     return.indices      = FALSE) {

  ## --- arg validation --------------------------------------------------
  if (missing(sequences))
    stop("`sequences` is required", call. = FALSE)
  if (is.null(universe) && is.null(genome))
    stop("supply exactly one of `universe` or `genome`", call. = FALSE)
  if (!is.null(universe) && !is.null(genome))
    stop("`universe` and `genome` are mutually exclusive; supply only one",
         call. = FALSE)
  if (!is(sequences, "XStringSet"))
    stop("`sequences` must be an XStringSet object", call. = FALSE)
  if (!is.null(covariates) && !is.null(genome))
    stop("`covariates` cannot be combined with `genome`: internally-",
         "sampled genomic windows have no covariate values. Build the ",
         "universe yourself (sample windows, compute covariates per ",
         "window) and pass an explicit `universe` + `universe.covariates`.",
         call. = FALSE)

  ## --- BSgenome shortcut: build the universe internally ---------------
  if (!is.null(genome)) {
    universe <- sample_genomic_universe(genome, sequences,
                                         n.candidates = n.candidates,
                                         chromosomes  = chromosomes,
                                         exclude      = exclude)
  }

  if (!is(universe, "XStringSet"))
    stop("`universe` must be an XStringSet object", call. = FALSE)
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

  ## --- covariate validation / coercion ---------------------------------
  cov <- validate_covariates(covariates, universe.covariates,
                             length(sequences), length(universe),
                             n.bins.covariates, bin.type.covariates)

  ## --- matching axes ----------------------------------------------------
  ## Axis 1 = GC (equal-width over [0, 1]), axis 2 = length (quantile),
  ## built the same way whether or not covariates are added.
  ax.gc  <- make_match_axis(gc.t,  gc.u,  n.bins.gc,     "equalwidth", "gc",
                            range = c(0, 1))
  ax.len <- make_match_axis(len.t, len.u, n.bins.length, "quantile",   "length")

  n.targets <- length(sequences)

  ## Build the cell map plus a per-target ring-candidate gatherer
  ## (`ring_fun`) and matrix of target home coordinates (`homes`). With
  ## no covariates this is the fast two-axis (GC x length) path: the
  ## joint bin index is computed inline and the rings are the plain 2-D
  ## Manhattan rings, identical to the original algorithm. Covariates
  ## switch to the general N-axis path (mixed-radix index, N-dimensional
  ## rings), which carries a small constant-factor overhead but is only
  ## reached when covariates are actually supplied. Both feed the same
  ## per-target loop below.
  if (is.null(cov)) {
    n.gc.bins  <- ax.gc$n.bins
    n.len.bins <- ax.len$n.bins
    bin_idx_u  <- (ax.gc$bin.u - 1L) * n.len.bins + ax.len$bin.u
    cells      <- split(seq_along(ax.gc$bin.u), bin_idx_u)
    homes      <- cbind(ax.gc$bin.t, ax.len$bin.t)
    max.ring   <- MAX_LOCAL_RING
    ring_fun   <- function(home, ring)
      ring_candidates(home[1L], home[2L], ring, cells, n.gc.bins, n.len.bins)
  } else {
    axes <- list(ax.gc, ax.len)
    for (j in seq_along(cov$names)) {
      nm <- cov$names[j]
      ax <- make_match_axis(cov$cov.t[[nm]], cov$cov.u[[nm]],
                            cov$n.bins[j], cov$bin.type[j], nm)
      if (ax$n.bins < 2L)
        message("covariate `", nm, "` is constant; not used for matching")
      axes[[length(axes) + 1L]] <- ax
    }
    n.axes     <- length(axes)
    n.bins.vec <- vapply(axes, function(a) a$n.bins, integer(1))
    homes      <- do.call(cbind, lapply(axes, function(a) a$bin.t))
    bins.u     <- do.call(cbind, lapply(axes, function(a) a$bin.u))
    bin_idx_u  <- encode_bins(bins.u, n.bins.vec)
    cells      <- split(seq_len(length(universe)), bin_idx_u)
    ## A fixed radius covers far less ground per axis as dimensions grow,
    ## so raise the ring cap with the number of axes (floor of two keeps
    ## the two-axis default unchanged).
    max.ring   <- max(MAX_LOCAL_RING, n.axes)
    ring_fun   <- function(home, ring) ring_cells(home, ring, cells, n.bins.vec)
  }

  ## --- per-target matching (shared) ------------------------------------
  ## RNG comes from the caller's global state (set.seed() before the call
  ## for reproducibility).
  matched_u     <- integer(n.targets * n.per.target)
  used_global   <- if (unique) integer(0) else NULL
  fallback_used <- 0L

  for (ti in seq_len(n.targets)) {
    cand <- find_candidates(homes[ti, ], ring_fun, cells,
                            used_global, n.per.target, unique, max.ring)
    if (cand$fellback) fallback_used <- fallback_used + 1L
    drawn <- cand$picked
    matched_u[((ti - 1L) * n.per.target + 1L):(ti * n.per.target)] <- drawn
    if (unique) used_global <- c(used_global, drawn)
  }

  if (fallback_used > 0L)
    warning(fallback_used, " target(s) required ring-expansion fallback ",
            "(empty or undersized bin). Consider widening the n.bins.* ",
            "arguments or providing a bigger universe.", call. = FALSE)

  ## --- assemble output --------------------------------------------------
  if (return.indices) {
    target.i <- rep(seq_len(n.targets), each = n.per.target)
    df <- data.frame(
      target.i        = target.i,
      target.gc       = gc.t[target.i],
      target.length   = len.t[target.i],
      universe.i      = matched_u,
      universe.gc     = gc.u[matched_u],
      universe.length = len.u[matched_u],
      stringsAsFactors = FALSE
    )
    if (!is.null(cov)) {
      cov.cols <- list()
      for (nm in cov$names) {
        cov.cols[[paste0("target.", nm)]]   <- cov$cov.t[[nm]][target.i]
        cov.cols[[paste0("universe.", nm)]] <- cov$cov.u[[nm]][matched_u]
      }
      df <- cbind(df, data.frame(cov.cols, stringsAsFactors = FALSE,
                                 check.names = FALSE))
    }
    return(df)
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
#' Companion to [match_bkg()]. Overlays density curves for the target
#' and matched-background sets in a faceted `ggplot`, useful for
#' visually verifying that the binned matching landed where it was
#' supposed to. Called with two `XStringSet`s it plots GC fraction and
#' sequence length; called with the `indices` data.frame from
#' `match_bkg(..., return.indices = TRUE)` it plots every matched axis,
#' including any covariates.
#'
#' @param sequences `XStringSet` or `NULL`. The target sequences passed
#'   to [match_bkg()]. Ignored when `indices` is supplied.
#' @param bkg `XStringSet` or `NULL`. The matched background sequences
#'   returned by [match_bkg()]. Ignored when `indices` is supplied.
#' @param by `character` or `NULL`. Axes to plot. For the
#'   `sequences`/`bkg` form, a subset of `c("gc", "length")` (default
#'   both). For the `indices` form, a subset of the available axis names
#'   (`"gc"`, `"length"`, and each covariate name); `NULL` (default)
#'   plots them all.
#' @param bins `integer(1)`. Ignored when geom = "density"; kept for
#'   API parity with other v2 plotting helpers. Default `30L`.
#' @param indices `data.frame` or `NULL`. The output of
#'   `match_bkg(..., return.indices = TRUE)`. When supplied, every
#'   `target.<axis>` / `universe.<axis>` column pair (apart from the
#'   index column) becomes a density panel, so covariate matches can be
#'   inspected alongside GC and length. Default `NULL`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' data(ArabidopsisPromoters)
#' target   <- ArabidopsisPromoters[1:10]
#' universe <- ArabidopsisPromoters[11:50]
#' bkg <- match_bkg(target, universe, n.per.target = 2)
#'
#' ## GC and length density of the target vs. its matched background:
#' plot_match_bkg(target, bkg)
#'
#' @seealso [match_bkg()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_match_bkg <- function(sequences = NULL, bkg = NULL,
                           by = NULL, bins = 30L, indices = NULL) {
  label_for <- function(a)
    switch(a, gc = "GC fraction", length = "Sequence length", a)

  if (!is.null(indices)) {
    if (!is.data.frame(indices))
      stop("`indices` must be the data.frame returned by ",
           "match_bkg(..., return.indices = TRUE)", call. = FALSE)
    axes <- setdiff(sub("^target\\.", "",
                        grep("^target\\.", names(indices), value = TRUE)),
                    "i")
    if (!is.null(by)) axes <- intersect(axes, by)
    axes <- axes[vapply(axes, function(a)
      all(c(paste0("target.", a), paste0("universe.", a)) %in% names(indices)),
      logical(1))]
    if (!length(axes))
      stop("no plottable target.*/universe.* column pairs found in `indices`",
           call. = FALSE)
    pieces <- lapply(axes, function(a) {
      lab <- label_for(a)
      rbind(
        data.frame(set = "target",  axis = lab,
                   value = indices[[paste0("target.", a)]],
                   stringsAsFactors = FALSE),
        data.frame(set = "matched", axis = lab,
                   value = indices[[paste0("universe.", a)]],
                   stringsAsFactors = FALSE)
      )
    })
    long2 <- do.call(rbind, pieces)
    long2$axis <- factor(long2$axis,
                         levels = vapply(axes, label_for, character(1)))
    n.panels <- length(axes)
  } else {
    if (!is(sequences, "XStringSet") || !is(bkg, "XStringSet"))
      stop("both `sequences` and `bkg` must be XStringSet objects (or ",
           "supply `indices`)", call. = FALSE)
    if (is.null(by)) by <- c("gc", "length")
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
    n.panels <- length(by)
  }

  long2$set <- factor(long2$set, levels = c("target", "matched"))

  ggplot2::ggplot(long2, ggplot2::aes(x = .data$value, colour = .data$set,
                                       fill = .data$set)) +
    ggplot2::geom_density(alpha = 0.25, linewidth = 0.6) +
    ggplot2::facet_wrap(~ axis, scales = "free", ncol = n.panels) +
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

## Ring-expanding candidate finder. Returns picked indices + a logical
## `fellback` flag (TRUE only when the with-replacement fallback below
## was actually used). `ring_fun(home, ring)` returns the universe
## indices at Manhattan distance `ring` from the target's home cell;
## the caller supplies the two-axis or N-axis gatherer, so the
## fallback logic here is shared between both paths.
##
## Try the home cell first, then expand outward by at most `max.ring`
## Manhattan-distance rings. If we still don't have enough candidates
## after that (typical when `unique = TRUE` and the target distribution
## is narrow), fall back to with-replacement sampling from the
## cumulative *local* pool. This keeps the composition match tight at
## the cost of letting a few universe sequences repeat, rather than
## expanding to far-distant (compositionally unmatched) cells.
##
## MAX_LOCAL_RING is the two-axis default and the floor used by the
## caller (which raises the cap to the number of axes for higher
## dimensions).
MAX_LOCAL_RING <- 2L
find_candidates <- function(home, ring_fun, cells, used_global,
                            n.needed, unique, max.ring) {
  ring <- 0L
  picked     <- integer(0)
  local_pool <- integer(0)
  fellback   <- FALSE
  repeat {
    cand <- ring_fun(home, ring)
    local_pool <- c(local_pool, cand)
    if (unique && length(used_global) > 0L)
      cand <- setdiff(cand, used_global)
    cand <- setdiff(cand, picked)
    if (length(cand) > 0L) {
      take <- min(length(cand), n.needed - length(picked))
      ## Use cand[sample.int(...)] instead of sample(cand, ...) to
      ## sidestep R's `sample()` length-1 gotcha (sample(N, 1) is
      ## sample.int(N, 1) when N is a single number, not "return N").
      picked <- c(picked, cand[sample.int(length(cand), take,
                                          replace = FALSE)])
    }
    if (length(picked) >= n.needed) break
    ring <- ring + 1L
    if (ring > max.ring) {
      ## Cap reached. Sample with replacement from the accumulated
      ## local pool (preserves composition match). If that's somehow
      ## empty too, fall back to global random as a last resort.
      need <- n.needed - length(picked)
      if (need > 0L) {
        fellback <- TRUE
        pool <- unique(local_pool)
        if (length(pool) == 0L)
          pool <- unlist(cells, use.names = FALSE)
        picked <- c(picked, pool[sample.int(length(pool), need,
                                            replace = TRUE)])
      }
      break
    }
  }
  list(picked = picked, fellback = fellback)
}

## Sample a candidate-universe XStringSet from a BSgenome by drawing
## `n.candidates` random windows whose widths match the target sequence
## width distribution, optionally excluding any window that overlaps
## `exclude`. Returns an XStringSet ready to be passed to the standard
## bin-and-match path below.
sample_genomic_universe <- function(genome, sequences,
                                     n.candidates = NULL,
                                     chromosomes  = NULL,
                                     exclude      = NULL) {
  for (pkg in c("BSgenome", "GenomicRanges", "IRanges", "GenomeInfoDb"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("package `", pkg, "` is required when `genome` is supplied",
           call. = FALSE)
  if (!is(genome, "BSgenome"))
    stop("`genome` must be a BSgenome object", call. = FALSE)
  if (!is.null(exclude) && !is(exclude, "GRanges"))
    stop("`exclude` must be a GRanges object or NULL", call. = FALSE)

  if (is.null(chromosomes)) {
    chromosomes <- tryCatch(
      GenomeInfoDb::standardChromosomes(genome),
      error = function(e) GenomeInfoDb::seqnames(genome)
    )
  }
  if (!length(chromosomes))
    stop("no chromosomes available in `genome`", call. = FALSE)

  chr.lens <- GenomeInfoDb::seqlengths(genome)[chromosomes]
  chr.lens <- chr.lens[!is.na(chr.lens) & chr.lens > 0L]
  if (!length(chr.lens))
    stop("no usable chromosomes (all NA / zero length)", call. = FALSE)

  target.widths <- as.integer(width(sequences))
  if (any(target.widths < 1L))
    stop("`sequences` contains zero-width entries", call. = FALSE)

  if (is.null(n.candidates))
    n.candidates <- as.integer(max(10000L, length(sequences) * 10L))
  n.candidates <- as.integer(n.candidates)
  if (n.candidates < length(sequences))
    stop("`n.candidates` (", n.candidates, ") must be at least ",
         "length(sequences) (", length(sequences), ")", call. = FALSE)

  ## Oversample to absorb exclude losses.
  oversample <- if (is.null(exclude)) 1.0 else 2.0
  n.try <- as.integer(ceiling(n.candidates * oversample))

  chr.sample    <- sample(names(chr.lens), n.try, replace = TRUE,
                          prob = as.numeric(chr.lens) / sum(chr.lens))
  width.sample  <- sample(target.widths, n.try, replace = TRUE)
  max.start     <- as.numeric(chr.lens[chr.sample]) - width.sample + 1
  fits          <- max.start >= 1
  chr.sample    <- chr.sample[fits]
  width.sample  <- width.sample[fits]
  max.start     <- max.start[fits]
  start.sample  <- as.integer(floor(stats::runif(length(max.start)) *
                                      max.start)) + 1L

  gr <- GenomicRanges::GRanges(
    seqnames = chr.sample,
    ranges   = IRanges::IRanges(start = start.sample, width = width.sample)
  )

  if (!is.null(exclude))
    gr <- IRanges::subsetByOverlaps(gr, exclude, invert = TRUE)

  if (length(gr) < n.candidates) {
    if (length(gr) < length(sequences))
      stop("could not sample enough non-excluded windows ",
           "(got ", length(gr), ", need at least ", length(sequences),
           "); try a larger `n.candidates` or a smaller `exclude`.",
           call. = FALSE)
  } else {
    gr <- gr[seq_len(n.candidates)]
  }

  Biostrings::getSeq(genome, gr)
}

## Build one matching axis. `type` is "equalwidth" (equal-width breaks;
## pass an explicit `range` to pin them, as GC does over [0, 1]) or
## "quantile" (breaks from the universe values' quantiles, like
## length). Mirrors the original per-axis GC and length binning,
## including the constant-axis collapse to a single bin. Returns a list
## with per-target/per-universe 1-based bin indices and the bin count.
make_match_axis <- function(values.t, values.u, n.bins, type, name,
                            range = NULL) {
  n.bins <- as.integer(n.bins)
  if (type == "equalwidth") {
    if (is.null(range)) range <- range(c(values.t, values.u))
    breaks <- seq(range[1L], range[2L], length.out = n.bins + 1L)
    ## Degenerate (constant) range: collapse to a single bin.
    if (length(unique(breaks)) < 2L) {
      lo <- min(c(values.t, values.u))
      hi <- max(c(values.t, values.u))
      breaks <- c(lo - 1, hi + 1)
    }
  } else {
    ## Quantile-based breaks; clamp uniques to avoid degenerate breaks
    ## when many values tie. If the universe is constant on this axis,
    ## collapse to a single bin (nothing to match on).
    breaks <- unique(stats::quantile(values.u,
                                     probs = seq(0, 1,
                                                 length.out = n.bins + 1L),
                                     names = FALSE))
    lo <- min(values.t, values.u)
    hi <- max(values.t, values.u)
    if (length(breaks) < 2L) {
      breaks <- c(lo - 1, hi + 1)
    } else {
      breaks[1]              <- min(breaks[1], lo) - 1
      breaks[length(breaks)] <- max(breaks[length(breaks)], hi) + 1
    }
  }
  bin.t <- pmin(pmax(findInterval(values.t, breaks, rightmost.closed = TRUE),
                     1L), length(breaks) - 1L)
  bin.u <- pmin(pmax(findInterval(values.u, breaks, rightmost.closed = TRUE),
                     1L), length(breaks) - 1L)
  list(name = name, bin.t = as.integer(bin.t), bin.u = as.integer(bin.u),
       n.bins = length(breaks) - 1L)
}

## Mixed-radix linear bin index for the N-axis path, axis 1
## slowest-varying. `bins` is an n x n.axes matrix of 1-based per-axis
## bin numbers (a length-n.axes vector is accepted as a single row).
## For two axes this is exactly (bins[, 1] - 1) * n.bins.vec[2] +
## bins[, 2], the same index the two-axis path computes inline.
encode_bins <- function(bins, n.bins.vec) {
  n.ax <- length(n.bins.vec)
  if (!is.matrix(bins)) bins <- matrix(as.integer(bins), ncol = n.ax)
  key <- bins[, 1L] - 1L
  if (n.ax >= 2L)
    for (a in 2:n.ax) key <- key * n.bins.vec[a] + (bins[, a] - 1L)
  as.integer(key + 1L)
}

## All integer offset vectors of length `n.axes` whose absolute values
## sum to `ring`: the n-dimensional generalisation of a Manhattan-
## distance ring. `ring == 0` is the single zero vector (home cell).
## A nonzero allocation to an axis is emitted as both +/- (a single 0
## otherwise), so no offset vector is generated twice.
ring_offsets <- function(ring, n.axes) {
  if (ring == 0L) return(list(rep(0L, n.axes)))
  out <- list()
  recurse <- function(axis, remaining, acc) {
    if (axis == n.axes) {
      if (remaining == 0L) {
        out[[length(out) + 1L]] <<- c(acc, 0L)
      } else {
        out[[length(out) + 1L]] <<- c(acc,  remaining)
        out[[length(out) + 1L]] <<- c(acc, -remaining)
      }
      return(invisible(NULL))
    }
    for (m in 0:remaining) {
      if (m == 0L) {
        recurse(axis + 1L, remaining, c(acc, 0L))
      } else {
        recurse(axis + 1L, remaining - m, c(acc,  m))
        recurse(axis + 1L, remaining - m, c(acc, -m))
      }
    }
  }
  recurse(1L, ring, integer(0))
  out
}

## N-axis ring gatherer: universe-index candidates in all cells at
## Manhattan distance `ring` from the target's home bin coordinate
## `home` (a length-n.axes integer vector). `ring == 0` returns the
## home cell. The cell key is computed inline (same mixed-radix formula
## as encode_bins()) to avoid allocating a matrix per offset.
ring_cells <- function(home, ring, cells, n.bins.vec) {
  offs <- ring_offsets(ring, length(home))
  n.ax <- length(n.bins.vec)
  out <- integer(0)
  for (d in offs) {
    coord <- home + d
    if (any(coord < 1L) || any(coord > n.bins.vec)) next
    key <- coord[1L] - 1L
    for (a in 2:n.ax) key <- key * n.bins.vec[a] + (coord[a] - 1L)
    cell <- cells[[as.character(key + 1L)]]
    if (!is.null(cell)) out <- c(out, cell)
  }
  unique(out)
}

## Two-axis ring gatherer: universe-index candidates at Manhattan
## distance `ring` from (gc.b, len.b). `ring == 0` returns the home
## cell. Keys are computed inline; this is the original fast path used
## when no covariates are supplied.
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

## Coerce a covariate argument (data.frame / matrix / numeric vector)
## to a data.frame with named numeric columns.
coerce_covariates <- function(x, argname) {
  if (is.data.frame(x)) {
    out <- x
  } else if (is.matrix(x)) {
    out <- as.data.frame(x, stringsAsFactors = FALSE)
    if (is.null(colnames(x)))
      names(out) <- paste0("covariate", seq_len(ncol(x)))
  } else if (is.numeric(x) && is.null(dim(x))) {
    out <- data.frame(covariate1 = x, stringsAsFactors = FALSE)
  } else {
    stop("`", argname, "` must be a data.frame, matrix, or numeric vector",
         call. = FALSE)
  }
  out
}

## Resolve a per-covariate argument: a single value applied to every
## covariate, or a vector named by covariate column (naming all of
## them).
resolve_per_cov <- function(x, nm, argname, char = FALSE) {
  conv <- function(v) if (char) as.character(v) else as.integer(v)
  if (!is.null(names(x)) && any(nzchar(names(x)))) {
    if (!setequal(names(x), nm))
      stop("named `", argname, "` must name every covariate exactly; got {",
           paste(names(x), collapse = ", "), "} for covariates {",
           paste(nm, collapse = ", "), "}", call. = FALSE)
    out <- conv(x[nm])
  } else if (length(x) == 1L) {
    out <- rep(conv(x), length(nm))
  } else {
    stop("`", argname, "` must be a single value or a vector named by ",
         "covariate", call. = FALSE)
  }
  names(out) <- nm
  out
}

## Validate and coerce the covariate tables. Returns NULL when no
## covariates are supplied, else a list with the coerced target/universe
## tables, per-covariate bin counts and bin types, and column names.
validate_covariates <- function(covariates, universe.covariates,
                                n.targets, n.universe,
                                n.bins.covariates, bin.type.covariates) {
  if (is.null(covariates) && is.null(universe.covariates))
    return(NULL)
  if (is.null(covariates))
    stop("`universe.covariates` supplied without `covariates`", call. = FALSE)
  if (is.null(universe.covariates))
    stop("`covariates` requires a matching `universe.covariates`",
         call. = FALSE)

  cov.t <- coerce_covariates(covariates, "covariates")
  cov.u <- coerce_covariates(universe.covariates, "universe.covariates")

  if (nrow(cov.t) != n.targets)
    stop("`covariates` has ", nrow(cov.t), " row(s) but `sequences` has ",
         n.targets, call. = FALSE)
  if (nrow(cov.u) != n.universe)
    stop("`universe.covariates` has ", nrow(cov.u), " row(s) but `universe` ",
         "has ", n.universe, call. = FALSE)
  if (ncol(cov.t) < 1L)
    stop("`covariates` must have at least one column", call. = FALSE)
  if (ncol(cov.t) != ncol(cov.u))
    stop("`covariates` and `universe.covariates` must have the same number ",
         "of columns", call. = FALSE)
  if (!identical(names(cov.t), names(cov.u)))
    stop("`covariates` and `universe.covariates` must have identical column ",
         "names in the same order; got (",
         paste(names(cov.t), collapse = ", "), ") vs (",
         paste(names(cov.u), collapse = ", "), ")", call. = FALSE)

  reserved <- c("i", "gc", "length")
  bad <- intersect(names(cov.t), reserved)
  if (length(bad))
    stop("covariate column name(s) {", paste(bad, collapse = ", "),
         "} are reserved (collide with the return.indices columns); ",
         "rename them", call. = FALSE)

  for (nm in names(cov.t)) {
    if (!is.numeric(cov.t[[nm]]) || !is.numeric(cov.u[[nm]]))
      stop("covariate column `", nm, "` must be numeric", call. = FALSE)
    if (any(!is.finite(cov.t[[nm]])) || any(!is.finite(cov.u[[nm]])))
      stop("covariate column `", nm, "` contains NA/NaN/Inf; impute or drop ",
           "before calling match_bkg()", call. = FALSE)
  }

  nm <- names(cov.t)
  nb <- resolve_per_cov(n.bins.covariates, nm, "n.bins.covariates")
  if (any(nb < 2L))
    stop("`n.bins.covariates` must be >= 2 for every covariate", call. = FALSE)
  bt <- resolve_per_cov(bin.type.covariates, nm, "bin.type.covariates",
                        char = TRUE)
  bt <- vapply(bt, function(x) match.arg(x, c("quantile", "equalwidth")),
               character(1))

  list(cov.t = cov.t, cov.u = cov.u, n.bins = nb, bin.type = bt, names = nm)
}
