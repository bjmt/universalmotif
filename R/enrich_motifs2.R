#' Fast minimalist motif enrichment.
#'
#' `enrich_motifs2()` is a deliberately stripped-down counterpart to
#' [enrich_motifs()] whose default surface mirrors the command-line tool
#' [yamtk](https://github.com/bjmt/yamtk). It exposes a
#' single p-value cutoff for hits, a single q-value cutoff for results, two
#' test modes (`"seqs"` and `"sites"`, both Fisher's exact), and uses
#' [scan_sequences2()] under the hood to scan target and background
#' sequences. P-value adjustment is hard-coded to Benjamini-Hochberg.
#'
#' Use [enrich_motifs()] when you need any of: q-value adjustment methods
#' other than BH; multifreq / gapped motifs; multiple-testing E-values;
#' the `respect.strand`, `allow.nonfinite`, or non-pvalue `threshold.type`
#' machinery; or amino-acid motifs.
#'
#' @param motifs See [convert_motifs()] for accepted motif formats. DNA or
#'   RNA only; amino-acid and custom alphabet motifs are rejected.
#' @param sequences `XStringSet`. Target sequences.
#' @param bkg.sequences `XStringSet` or `NULL`. Background sequences. If
#'   `NULL` (default), target sequences are shuffled k-let-conserving via
#'   [shuffle_sequences()] with `k = shuffle.k`, matching both
#'   `enrich_motifs()` and `yamtk enr` default behaviour.
#' @param pvalue `numeric(1)`. P-value cutoff passed to [scan_sequences2()]
#'   for reporting individual hits. Default `1e-4`.
#' @param qvalue `numeric(1)`. Result-level q-value cutoff. Motifs whose
#'   BH-adjusted enrichment p-value exceeds this are filtered out. Default
#'   `0.1`.
#' @param test `character(1)`. One of `"seqs"` (default) or `"sites"`.
#'   `"seqs"` runs Fisher's exact on per-sequence hit presence (does this
#'   sequence contain >=1 hit?), equivalent to [enrich_motifs()]'s
#'   `mode = "seq.hits"`. `"sites"` runs Fisher's exact on per-position
#'   rates (how many hits per scannable position?), equivalent to
#'   `mode = "total.hits"`.
#' @param RC `logical(1)`. If `TRUE` (default), scan both strands.
#' @param shuffle.k `integer(1)`. K-let size for the background shuffle
#'   (only used when `bkg.sequences = NULL`). Default `2L`.
#' @param rng.seed `integer(1)`. RNG seed for background shuffling.
#'   Default `sample.int(1e6, 1)`; set explicitly for reproducible runs.
#' @param pseudocount `integer(1)`. Pseudocount added uniformly to each
#'   cell of the Fisher 2x2 table; useful when one side has zero counts.
#'   Default `1L` (matches yamtk's `-p 1` default).
#' @param nthreads `numeric(1)`. Number of threads passed to
#'   [scan_sequences2()] and [shuffle_sequences()]. `nthreads = 0` uses
#'   all available threads.
#'
#' @return A `data.frame` with one row per significant motif, sorted by
#'   `pvalue` ascending. Columns:
#'
#' \itemize{
#'   \item `motif`: motif name
#'   \item `motif.i`: 1-based index into the input `motifs`
#'   \item `consensus`: IUPAC consensus of the motif
#'   \item `target.seq.n`: number of target sequences scanned
#'   \item `target.seq.hits`: number of target sequences containing >=1 hit
#'   \item `target.site.hits`: total hit count across all target sequences
#'   \item `bkg.seq.n`: number of background sequences scanned
#'   \item `bkg.seq.hits`: number of background sequences containing >=1 hit
#'   \item `bkg.site.hits`: total hit count across all background sequences
#'   \item `enrichment`: fold-enrichment (rate ratio; floor of 0.5/n in the
#'     denominator to avoid division by zero), per yamtk's formula
#'   \item `log2.enrichment`: `log2(enrichment)`
#'   \item `pvalue`: Fisher's exact one-sided ("greater") p-value
#'   \item `qvalue`: Benjamini-Hochberg adjusted p-value across motifs
#' }
#'
#' @details
#' Internally, [scan_sequences2()] is run once on the target set and once on
#' the background set; the resulting hit tables are aggregated per motif into
#' a 2x2 contingency table whose entries depend on `test`:
#'
#' \describe{
#'   \item{`"seqs"`}{
#'     `a` = target sequences with >=1 hit;
#'     `b` = background sequences with >=1 hit;
#'     `c` = target sequences with no hits;
#'     `d` = background sequences with no hits.
#'   }
#'   \item{`"sites"`}{
#'     `a` = total hits in target;
#'     `b` = total hits in background;
#'     `c` = scannable target positions minus `a`;
#'     `d` = scannable background positions minus `b`.
#'     Scannable positions are computed per motif as
#'     `sum(pmax(width(seqs) - motif_width + 1, 0))`, doubled if `RC = TRUE`.
#'   }
#' }
#'
#' Each cell has `pseudocount` added before [stats::fisher.test()] is run
#' with `alternative = "greater"`. Q-values are computed via
#' `p.adjust(method = "BH")` across all input motifs (not just the ones
#' that pass the filter), then rows with `qvalue > qvalue` are dropped.
#'
#' @references
#'
#' Tremblay BJM (2026). yamtk: Yet Another Motif ToolKit.
#' \url{https://github.com/bjmt/yamtk}.
#'
#' McLeay R, Bailey TL (2010). "Motif Enrichment Analysis: A unified
#' framework and method evaluation." *BMC Bioinformatics*, **11**.
#'
#' @examples
#' library(universalmotif)
#' data(ArabidopsisPromoters)
#' data(ArabidopsisMotif)
#' if (R.Version()$arch != "i386") {
#'   enrich_motifs2(ArabidopsisMotif, ArabidopsisPromoters,
#'                  pvalue = 1e-3, qvalue = 0.5, rng.seed = 1)
#' }
#'
#' @seealso [enrich_motifs()], [scan_sequences2()], [shuffle_sequences()],
#'     [match_bkg()], [motif_pvalue()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
enrich_motifs2 <- function(motifs, sequences, bkg.sequences = NULL,
                           pvalue = 1e-4, qvalue = 0.1,
                           test = c("seqs", "sites"),
                           RC = TRUE, shuffle.k = 2L,
                           rng.seed = sample.int(1e6, 1),
                           pseudocount = 1L, nthreads = 1) {

  ## --- argument validation ---------------------------------------------
  if (missing(motifs) || missing(sequences))
    stop("need both `motifs` and `sequences`", call. = FALSE)
  test <- match.arg(test)
  if (!is.numeric(pvalue) || length(pvalue) != 1L || is.na(pvalue) ||
      pvalue <= 0 || pvalue >= 1)
    stop("`pvalue` must be a single numeric in (0, 1)", call. = FALSE)
  if (!is.numeric(qvalue) || length(qvalue) != 1L || is.na(qvalue) ||
      qvalue <= 0 || qvalue > 1)
    stop("`qvalue` must be a single numeric in (0, 1]", call. = FALSE)
  if (!isTRUEorFALSE(RC))
    stop("`RC` must be a single logical", call. = FALSE)
  if (!is.numeric(shuffle.k) || length(shuffle.k) != 1L || shuffle.k < 1L)
    stop("`shuffle.k` must be a positive integer", call. = FALSE)
  shuffle.k <- as.integer(shuffle.k)
  if (!is.numeric(rng.seed) || length(rng.seed) != 1L || is.na(rng.seed))
    stop("`rng.seed` must be a single numeric", call. = FALSE)
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      is.na(pseudocount) || pseudocount < 0)
    stop("`pseudocount` must be a non-negative numeric", call. = FALSE)

  nthreads <- resolve_nthreads(nthreads)

  ## --- normalise motifs ------------------------------------------------
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) != 1L)
    stop("all motifs must share the same alphabet", call. = FALSE)
  mot.alph <- unique(mot.alphs)
  if (!mot.alph %in% c("DNA", "RNA"))
    stop("`enrich_motifs2()` only supports DNA/RNA motifs; got `",
         mot.alph, "`. Use `enrich_motifs()` for other alphabets.",
         call. = FALSE)

  mot.names  <- vapply(motifs, function(x) x@name, character(1))
  if (any(duplicated(mot.names)))
    mot.names <- make.unique(mot.names)
  mot.widths <- vapply(motifs, function(x) ncol(x@motif), integer(1))

  ## consensus comes straight from the motif's @consensus slot (always
  ## populated by universalmotif's validity machinery); fall back to
  ## building it column-by-column from the score matrix if the slot is
  ## empty for any reason.
  mot.cons <- vapply(motifs, function(x) {
    if (length(x@consensus) == 1L && nzchar(x@consensus))
      return(x@consensus)
    paste0(vapply(seq_len(ncol(x@motif)),
                  function(j) get_consensus(x@motif[, j], mot.alph,
                                            x@type, x@pseudocount),
                  character(1)),
           collapse = "")
  }, character(1))

  ## --- background sequences --------------------------------------------
  if (is.null(bkg.sequences)) {
    bkg.sequences <- shuffle_sequences(sequences, k = shuffle.k,
                                       method = "euler",
                                       nthreads = nthreads,
                                       rng.seed = rng.seed)
  }

  if (seqtype(sequences) != mot.alph)
    stop("motif alphabet (", mot.alph, ") and `sequences` alphabet (",
         seqtype(sequences), ") do not match", call. = FALSE)
  if (seqtype(bkg.sequences) != mot.alph)
    stop("motif alphabet (", mot.alph, ") and `bkg.sequences` alphabet (",
         seqtype(bkg.sequences), ") do not match", call. = FALSE)

  ## --- scan -------------------------------------------------------------
  hits.tgt <- scan_sequences2(motifs, sequences,  pvalue = pvalue, RC = RC,
                              nthreads = nthreads, return.granges = FALSE)
  hits.bkg <- scan_sequences2(motifs, bkg.sequences, pvalue = pvalue, RC = RC,
                              nthreads = nthreads, return.granges = FALSE)

  n.tgt.seq <- length(sequences)
  n.bkg.seq <- length(bkg.sequences)
  tgt.widths <- as.integer(width(sequences))
  bkg.widths <- as.integer(width(bkg.sequences))
  rc.mult <- if (RC) 2L else 1L

  ## --- per-motif aggregation -------------------------------------------
  n.motifs <- length(motifs)
  target.site.hits <- integer(n.motifs)
  target.seq.hits  <- integer(n.motifs)
  bkg.site.hits    <- integer(n.motifs)
  bkg.seq.hits     <- integer(n.motifs)

  if (nrow(hits.tgt) > 0L) {
    tg <- table(factor(hits.tgt$motif.i, levels = seq_len(n.motifs)))
    target.site.hits <- as.integer(tg)
    sg <- vapply(seq_len(n.motifs), function(i) {
      r <- hits.tgt$motif.i == i
      if (!any(r)) 0L
      else length(unique(hits.tgt$sequence.i[r]))
    }, integer(1))
    target.seq.hits <- sg
  }
  if (nrow(hits.bkg) > 0L) {
    bg <- table(factor(hits.bkg$motif.i, levels = seq_len(n.motifs)))
    bkg.site.hits <- as.integer(bg)
    sb <- vapply(seq_len(n.motifs), function(i) {
      r <- hits.bkg$motif.i == i
      if (!any(r)) 0L
      else length(unique(hits.bkg$sequence.i[r]))
    }, integer(1))
    bkg.seq.hits <- sb
  }

  ## --- per-motif Fisher's exact + enrichment -------------------------------
  pvalues <- numeric(n.motifs)
  enrichments <- numeric(n.motifs)
  l2effs  <- numeric(n.motifs)

  eps <- 1e-9
  for (i in seq_len(n.motifs)) {
    if (test == "seqs") {
      a <- target.seq.hits[i]
      b <- bkg.seq.hits[i]
      c <- n.tgt.seq - a
      d <- n.bkg.seq - b
      pos.rate <- a / n.tgt.seq
      neg.rate <- b / n.bkg.seq
      neg.den  <- if (neg.rate > eps) neg.rate else 0.5 / n.bkg.seq
      eff <- pos.rate / neg.den
      l2e <- if (eff > 0) log2(eff) else -Inf
    } else {
      ## sites mode: per-motif scannable position counts
      mw  <- mot.widths[i]
      ppos <- sum(pmax(tgt.widths - mw + 1L, 0L)) * rc.mult
      npos <- sum(pmax(bkg.widths - mw + 1L, 0L)) * rc.mult
      a <- target.site.hits[i]
      b <- bkg.site.hits[i]
      c <- if (a <= ppos) ppos - a else 0L
      d <- if (b <= npos) npos - b else 0L
      pos.rate <- if (ppos > 0) a / ppos else 0
      neg.den  <- if (npos > 0 && b > 0) b / npos else
                  0.5 / (if (npos > 0) npos else 1)
      eff <- if (pos.rate > 0) pos.rate / neg.den else 0
      l2e <- if (eff > 0) log2(eff) else -Inf
    }

    tbl <- matrix(c(a, c, b, d), nrow = 2L, byrow = FALSE) +
           as.integer(pseudocount)
    pvalues[i] <- fisher.test(tbl, alternative = "greater")$p.value
    enrichments[i] <- eff
    l2effs[i]  <- l2e
  }

  qvalues <- p.adjust(pvalues, method = "BH")

  ## --- assemble + filter ------------------------------------------------
  out <- data.frame(
    motif            = mot.names,
    motif.i          = seq_len(n.motifs),
    consensus        = mot.cons,
    target.seq.n     = n.tgt.seq,
    target.seq.hits  = target.seq.hits,
    target.site.hits = target.site.hits,
    bkg.seq.n        = n.bkg.seq,
    bkg.seq.hits     = bkg.seq.hits,
    bkg.site.hits    = bkg.site.hits,
    enrichment           = enrichments,
    log2.enrichment      = l2effs,
    pvalue           = pvalues,
    qvalue           = qvalues,
    stringsAsFactors = FALSE
  )

  out <- out[out$qvalue <= qvalue, , drop = FALSE]
  out <- out[order(out$pvalue, out$motif.i), , drop = FALSE]
  rownames(out) <- NULL
  out
}
