#' Discover _de novo_ motifs in a set of sequences.
#'
#' `motif_finder()` is a minimalist _de novo_ motif discovery
#' function. Its defaults mirror the command-line tool `yamtk me` (see
#' [yamtk](https://github.com/bjmt/yamtk)). The pipeline iterates over a
#' user-controlled range of motif widths; at each width it enumerates
#' over-represented k-mer seeds (Fisher's exact test on per-sequence
#' presence vs. a background set), aligns Hamming-1 neighbours of the best
#' seed to form a PPM, refines the PPM via two re-scan passes, accepts the
#' motif if its Fisher's exact p-value passes `stop.pvalue`, masks the
#' covered positions and repeats up to `nmotifs` times per width. After
#' all widths complete, motifs are cross-width-deduplicated by coverage
#' overlap, BH-adjusted, IC-trimmed at the flanks, and returned in
#' p-value order.
#'
#' @param sequences `XStringSet`. DNA or RNA target sequences.
#' @param bkg.sequences `XStringSet` or `NULL`. Background sequences.
#'   If `NULL` (default), target sequences are shuffled k-let-conserving
#'   via [shuffle_sequences()] with `k = shuffle.k`, matching `yamtk me`
#'   default behaviour.
#' @param min.width `integer(1)`. Minimum motif width. Default `6L`.
#' @param max.width `integer(1)`. Maximum motif width. Default `15L`.
#'   Hard ceiling 30.
#' @param nmotifs `integer(1)`. Maximum number of motifs to discover per
#'   width before moving on. Default `10L`.
#' @param hit.pvalue `numeric(1)`. P-value threshold for scoring a single
#'   hit during scanning. Default `1e-4`.
#' @param stop.pvalue `numeric(1)`. Per-motif Fisher's exact p-value below
#'   which a candidate motif is accepted. Default `1e-3`.
#' @param qvalue `numeric(1)`. Final BH-adjusted q-value cutoff. Motifs
#'   with `qvalue` above this are dropped. Default `1e-3`.
#' @param dedup.overlap `numeric(1)`. Cross-width dedup overlap fraction.
#'   When the intersection of two motifs' positional coverage divided by
#'   the smaller motif's coverage exceeds this, the higher-p-value motif
#'   is dropped. Default `0.5`.
#' @param RC `logical(1)`. If `TRUE` (default), scan both strands.
#' @param shuffle.k `integer(1)`. K-let size for the background shuffle
#'   (only used when `bkg.sequences = NULL`). Default `2L`.
#' @param rng.seed `integer(1)`. RNG seed for background shuffling.
#'   Default `sample.int(1e6, 1)`; set explicitly for reproducible runs.
#' @param pseudocount `integer(1)`. Pseudocount used when converting
#'   per-position counts to log-odds scores. Default `1L`.
#' @param nthreads `numeric(1)`. Number of threads. Parallelism is at
#'   the width level. `nthreads = 0` uses all available threads.
#'
#' @return A `universalmotif_df` (see [to_df()]) with one row per
#'   discovered motif, sorted by `pvalue` ascending. The standard
#'   `to_df()` columns (`name`, `altname`, `family`, `organism`,
#'   `consensus`, `alphabet`, `type`, `nsites`, `pval`, `eval`, `motif`,
#'   ...) are populated from each discovered motif. The yamtk-specific
#'   stats are carried as additional columns:
#'
#' \itemize{
#'   \item `rank`: 1-based discovery rank
#'   \item `width`: motif width (post IC-trim)
#'   \item `seqs_pos`: positive sequences with >=1 hit
#'   \item `seqs_neg`: background sequences with >=1 hit
#'   \item `sites_pos`: total target hits
#'   \item `sites_neg`: total background hits
#'   \item `n_pos`: number of positive sequences
#'   \item `n_neg`: number of background sequences
#'   \item `pvalue`: Fisher's exact one-sided p-value
#'   \item `qvalue`: Benjamini-Hochberg q-value
#' }
#'
#' Use [to_list()] on the returned object to recover a plain list of
#' `universalmotif` S4 objects.
#'
#' @details
#' Algorithm and defaults are a faithful port of `yamtk me`
#' (https://github.com/bjmt/yamtk), whose own design is in turn based
#' on the STREME algorithm (Bailey 2021): seed enumeration via word
#' counting with per-sequence Fisher's exact ranking, iterative PPM
#' refinement on positive sequences, per-motif Fisher's exact
#' significance against a shuffled or user-supplied background, and
#' position masking between discoveries.
#'
#' The C++ inner loops (seed enumeration, PPM refinement, per-motif
#' Fisher's evaluation, cross-width dedup) live in
#' `src/motif_finder.cpp` and are parallelised across motif widths via
#' `RcppThread::parallelFor`.
#'
#' Motifs with `qvalue > qvalue` are dropped from the result. To see
#' all discovered motifs regardless of significance, set
#' `qvalue = 1`.
#'
#' @references
#'
#' Bailey TL (2021). "STREME: accurate and versatile sequence motif
#' discovery." *Bioinformatics*, **37**(18), 2834-2840.
#' \doi{10.1093/bioinformatics/btab203}.
#'
#' Tremblay BJM (2026). yamtk: Yet Another Motif ToolKit.
#' \url{https://github.com/bjmt/yamtk}.
#'
#' @examples
#' \dontrun{
#' library(universalmotif)
#' library(Biostrings)
#' set.seed(1)
#' planted <- DNAString("TTGACATA")
#' seqs <- create_sequences(seqnum = 100, seqlen = 200, rng.seed = 1)
#' ## Plant the motif at random positions in 80% of sequences:
#' for (i in seq_len(80)) {
#'   pos <- sample.int(200 - 8, 1)
#'   seqs[[i]][pos:(pos + 7)] <- planted
#' }
#' motifs <- motif_finder(seqs, qvalue = 1, rng.seed = 1)
#' to_list(motifs)[[1]]
#' }
#'
#' @seealso [scan_sequences2()], [enrich_motifs2()], [shuffle_sequences()],
#'   [create_motif()], [to_df()], [to_list()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
motif_finder <- function(sequences, bkg.sequences = NULL,
                         min.width   = 6L,
                         max.width   = 15L,
                         nmotifs     = 10L,
                         hit.pvalue  = 1e-4,
                         stop.pvalue = 1e-3,
                         qvalue      = 1e-3,
                         dedup.overlap = 0.5,
                         RC          = TRUE,
                         shuffle.k   = 2L,
                         rng.seed    = sample.int(1e6, 1),
                         pseudocount = 1L,
                         nthreads    = 1) {

  ## --- argument validation ---------------------------------------------
  if (missing(sequences))
    stop("`sequences` is required", call. = FALSE)
  if (!is.numeric(min.width) || length(min.width) != 1L || min.width < 3L)
    stop("`min.width` must be a single integer >= 3", call. = FALSE)
  if (!is.numeric(max.width) || length(max.width) != 1L || max.width > 30L ||
      max.width < min.width)
    stop("`max.width` must be a single integer in [min.width, 30]",
         call. = FALSE)
  if (!is.numeric(nmotifs) || length(nmotifs) != 1L || nmotifs < 1L)
    stop("`nmotifs` must be a positive integer", call. = FALSE)
  if (!is.numeric(hit.pvalue) || length(hit.pvalue) != 1L || is.na(hit.pvalue) ||
      hit.pvalue <= 0 || hit.pvalue >= 1)
    stop("`hit.pvalue` must be a single numeric in (0, 1)", call. = FALSE)
  if (!is.numeric(stop.pvalue) || length(stop.pvalue) != 1L ||
      is.na(stop.pvalue) || stop.pvalue <= 0 || stop.pvalue > 1)
    stop("`stop.pvalue` must be a single numeric in (0, 1]", call. = FALSE)
  if (!is.numeric(qvalue) || length(qvalue) != 1L || is.na(qvalue) ||
      qvalue <= 0 || qvalue > 1)
    stop("`qvalue` must be a single numeric in (0, 1]", call. = FALSE)
  if (!is.numeric(dedup.overlap) || length(dedup.overlap) != 1L ||
      dedup.overlap < 0 || dedup.overlap > 1)
    stop("`dedup.overlap` must be a single numeric in [0, 1]", call. = FALSE)
  if (!isTRUEorFALSE(RC))
    stop("`RC` must be a single logical", call. = FALSE)
  if (!is.numeric(shuffle.k) || length(shuffle.k) != 1L || shuffle.k < 1L)
    stop("`shuffle.k` must be a positive integer", call. = FALSE)
  if (!is.numeric(rng.seed) || length(rng.seed) != 1L || is.na(rng.seed))
    stop("`rng.seed` must be a single numeric", call. = FALSE)
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L ||
      pseudocount < 0)
    stop("`pseudocount` must be a non-negative integer", call. = FALSE)

  nthreads <- resolve_nthreads(nthreads)

  ## --- alphabet ---------------------------------------------------------
  seq.alph <- seqtype(sequences)
  if (!seq.alph %in% c("DNA", "RNA"))
    stop("`motif_finder()` only supports DNA/RNA sequences; got `",
         seq.alph, "`.", call. = FALSE)

  ## --- background sequences --------------------------------------------
  if (is.null(bkg.sequences)) {
    bkg.sequences <- shuffle_sequences(sequences, k = as.integer(shuffle.k),
                                       method = "euler",
                                       nthreads = nthreads,
                                       rng.seed = rng.seed)
  }
  if (seqtype(bkg.sequences) != seq.alph)
    stop("`sequences` alphabet (", seq.alph,
         ") and `bkg.sequences` alphabet (", seqtype(bkg.sequences),
         ") do not match", call. = FALSE)

  ## --- background frequencies ------------------------------------------
  ## Compute background from the actual base composition of sequences +
  ## bkg.sequences combined, matching yamtk me's default (yamme.c
  ## compute_bkg_from_counts). The order is A, C, G, T/U.
  combined <- c(sequences, bkg.sequences)
  letter_set <- if (seq.alph == "DNA") c("A", "C", "G", "T")
                else                   c("A", "C", "G", "U")
  freqs <- Biostrings::letterFrequency(combined, letters = letter_set,
                                       OR = 0, as.prob = FALSE)
  totals <- colSums(freqs)
  bkg <- as.numeric(totals / sum(totals))
  ## Guard against zero / very low values (yamtk applies a 0.001 floor).
  bkg[bkg < 1e-3] <- 1e-3
  bkg <- bkg / sum(bkg)
  names(bkg) <- letter_set

  ## --- call C++ pipeline -----------------------------------------------
  raw <- motif_finder_cpp(as.character(sequences),
                          as.character(bkg.sequences),
                          as.integer(min.width),
                          as.integer(max.width),
                          as.integer(nmotifs),
                          as.numeric(hit.pvalue),
                          as.numeric(stop.pvalue),
                          as.numeric(dedup.overlap),
                          as.logical(RC),
                          as.numeric(bkg),
                          as.integer(pseudocount),
                          as.integer(nthreads))

  if (length(raw) == 0L) {
    return(empty_finder_df())
  }

  ## --- build universalmotif S4 objects ---------------------------------
  alph.letters <- if (seq.alph == "DNA") c("A","C","G","T") else c("A","C","G","U")
  mlist <- vector("list", length(raw))
  for (i in seq_along(raw)) {
    r <- raw[[i]]
    ## C++ returns ppm as a (width x 4) matrix; universalmotif expects
    ## (alphabet x width), so transpose.
    ppm <- t(r$ppm)
    rownames(ppm) <- alph.letters
    mlist[[i]] <- create_motif(ppm,
                               alphabet = seq.alph,
                               type     = "PPM",
                               name     = paste0("motif_", i),
                               nsites   = as.numeric(r$nsites),
                               pseudocount = as.numeric(pseudocount),
                               bkg      = bkg,
                               pval     = as.numeric(r$pvalue),
                               eval     = as.numeric(r$qvalue))
  }

  ## --- filter by qvalue ------------------------------------------------
  qvals <- vapply(raw, function(x) as.numeric(x$qvalue), numeric(1))
  keep  <- qvals <= qvalue
  if (!any(keep)) return(empty_finder_df())
  mlist <- mlist[keep]
  raw   <- raw[keep]

  ## --- build universalmotif_df + attach yamtk stats -------------------
  ## `to_df()` already pulls the motif's @consensus slot into the
  ## `consensus` column; we don't add a yamtk-style consensus to avoid
  ## clobbering universalmotif's slot-driven value.
  df <- to_df(mlist)
  df$rank      <- seq_along(mlist)
  df$width     <- vapply(raw, function(x) as.integer(x$width),     integer(1))
  df$seqs_pos  <- vapply(raw, function(x) as.integer(x$seqs_pos),  integer(1))
  df$seqs_neg  <- vapply(raw, function(x) as.integer(x$seqs_neg),  integer(1))
  df$sites_pos <- vapply(raw, function(x) as.integer(x$sites_pos), integer(1))
  df$sites_neg <- vapply(raw, function(x) as.integer(x$sites_neg), integer(1))
  df$n_pos     <- vapply(raw, function(x) as.integer(x$n_pos),     integer(1))
  df$n_neg     <- vapply(raw, function(x) as.integer(x$n_neg),     integer(1))
  df$pvalue    <- vapply(raw, function(x) as.numeric(x$pvalue),    numeric(1))
  df$qvalue    <- vapply(raw, function(x) as.numeric(x$qvalue),    numeric(1))
  df
}

# Empty-result return: a 0-row universalmotif_df-shaped data.frame with the
# documented yamtk-specific columns appended. We produce it by running to_df()
# on a single placeholder motif and then trimming to zero rows.
empty_finder_df <- function() {
  placeholder <- list(create_motif("A", name = "placeholder"))
  df <- to_df(placeholder)[0, , drop = FALSE]
  df$rank      <- integer(0)
  df$width     <- integer(0)
  df$seqs_pos  <- integer(0)
  df$seqs_neg  <- integer(0)
  df$sites_pos <- integer(0)
  df$sites_neg <- integer(0)
  df$n_pos     <- integer(0)
  df$n_neg     <- integer(0)
  df$pvalue    <- numeric(0)
  df$qvalue    <- numeric(0)
  df
}
