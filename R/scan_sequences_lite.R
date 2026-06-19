#' Minimalist motif scanner aligned with `yamtk scan` defaults.
#'
#' `scan_sequences_lite()` is a deliberately pared-down counterpart to
#' [scan_sequences()], with a default surface that mirrors the command-line
#' tool [yamtk](https://github.com/bjmt/yamtk). It exposes a single threshold
#' (a P-value), always scans both strands by default, always computes a
#' per-hit P-value, and returns either a `GRanges` (preferred, when
#' GenomicRanges is installed) or a `data.frame`. Use [scan_sequences()] when
#' you need any of multifreq scoring, gapped motifs, q-values, exhaustive
#' P-values, `respect.strand`, `allow.nonfinite`, or `threshold.type`s other
#' than a P-value.
#'
#' The scanner re-uses the very same C++ core as [scan_sequences()]; the
#' speed advantage simply comes from skipping the R-side bookkeeping for the
#' features this function does not support.
#'
#' P-values are computed via the dynamic-programming algorithm of
#' [motif_pvalue()] (FIMO-style; Grant et al. 2011), which is also what
#' `yamtk scan` uses.
#'
#' @param motifs See [convert_motifs()] for accepted motif formats. DNA or RNA
#'   only -- amino-acid and custom alphabet motifs are rejected.
#' @param sequences `XStringSet`. DNA or RNA sequences to scan; the alphabet
#'   must match the motif alphabet.
#' @param pvalue `numeric(1)`. P-value cutoff for reporting hits. Default
#'   `1e-4`.
#' @param RC `logical(1)`. If `TRUE` (default), scan both strands.
#' @param nthreads `numeric(1)`. Number of threads. `nthreads = 0` uses all
#'   available threads. Parallelisation happens inside `scan_sequences_cpp`
#'   over the motif x sequence cross-product.
#' @param return.granges `logical(1)` or `NULL`. When `NULL` (default), returns
#'   a `GRanges` if the `GenomicRanges` package is installed, otherwise a
#'   `data.frame`. Set `TRUE` / `FALSE` to force one or the other.
#' @param no.overlaps `logical(1)`. If `TRUE`, drop overlapping hits within
#'   each `(sequence, motif, strand)` group using a greedy
#'   algorithm (see [dedup_hits()]). Default `FALSE`.
#' @param no.overlaps.by `character(1)`. Which column to use as the priority
#'   for breaking ties between overlapping hits. `"pvalue"` (default) keeps
#'   the lowest p-value; `"score"` keeps the highest log-odds score. Only
#'   used when `no.overlaps = TRUE`.
#' @param no.overlaps.by.strand `logical(1)`. If `TRUE`, overlapping hits on
#'   opposite strands compete with each other. If `FALSE` (default), each
#'   strand is deduplicated independently.
#' @param no.overlaps.by.motif `logical(1)`. If `TRUE`, overlapping hits from
#'   different motifs compete with each other. If `FALSE` (default), each
#'   motif is deduplicated independently.
#'
#' @return Either a `GRanges` or a `data.frame`.
#'
#' Common fields (both shapes):
#' `motif`, `motif.i`, `sequence.i`, `score`, `score.pct`, `match`, `pvalue`.
#'
#' Coordinate fields:
#' - `GRanges`: `seqnames` (sequence name), `start`, `end`, `strand`.
#' - `data.frame`: `sequence`, `start`, `end`, `strand`.
#'
#' Coordinates are always 1-based, inclusive, with `start <= end` regardless of
#' strand. For hits on the `-` strand the `match` column is the reverse
#' complement of the sequence substring; i.e. the sequence as matched against
#' the motif. 
#'
#' @references
#'
#' Grant CE, Bailey TL, Noble WS (2011). "FIMO: scanning for occurrences of a
#' given motif." *Bioinformatics*, **27**(7), 1017-1018.
#'
#' Tremblay BJM (2026). yamtk: Yet Another Motif ToolKit.
#' \url{https://github.com/bjmt/yamtk}.
#'
#' @examples
#' library(universalmotif)
#' motifs <- create_motif(c("TATAAA", "CACGTG"))
#' seqs <- create_sequences(seqnum = 5, seqlen = 200, rng.seed = 1)
#' hits <- scan_sequences_lite(motifs, seqs, pvalue = 1e-3, return.granges = FALSE)
#' head(hits)
#'
#' @seealso [scan_sequences()], [motif_pvalue()], [convert_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @family lite motif functions
#' @export
scan_sequences_lite <- function(motifs, sequences, pvalue = 1e-4, RC = TRUE,
                            nthreads = 1, return.granges = NULL,
                            no.overlaps = FALSE,
                            no.overlaps.by = c("pvalue", "score"),
                            no.overlaps.by.strand = FALSE,
                            no.overlaps.by.motif  = FALSE) {

  if (missing(motifs) || missing(sequences))
    stop("need both `motifs` and `sequences`", call. = FALSE)
  if (!is.numeric(pvalue) || length(pvalue) != 1L || is.na(pvalue) ||
      pvalue <= 0 || pvalue >= 1)
    stop("`pvalue` must be a single numeric in (0, 1)", call. = FALSE)
  if (!isTRUEorFALSE(RC))
    stop("`RC` must be a single logical", call. = FALSE)
  if (!is.null(return.granges) && !isTRUEorFALSE(return.granges))
    stop("`return.granges` must be TRUE, FALSE, or NULL", call. = FALSE)
  if (!isTRUEorFALSE(no.overlaps))
    stop("`no.overlaps` must be a single logical", call. = FALSE)
  no.overlaps.by <- match.arg(no.overlaps.by)
  if (!isTRUEorFALSE(no.overlaps.by.strand))
    stop("`no.overlaps.by.strand` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(no.overlaps.by.motif))
    stop("`no.overlaps.by.motif` must be a single logical", call. = FALSE)

  nthreads <- resolve_nthreads(nthreads)

  ## --- normalise motifs ------------------------------------------------
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) != 1L)
    stop("all motifs must share the same alphabet", call. = FALSE)
  mot.alph <- unique(mot.alphs)
  if (!mot.alph %in% c("DNA", "RNA"))
    stop("`scan_sequences_lite()` only supports DNA/RNA motifs; got `",
         mot.alph, "`. Use `scan_sequences()` for other alphabets.",
         call. = FALSE)

  seq.alph <- seqtype(sequences)
  if (seq.alph != mot.alph)
    stop("motif alphabet (", mot.alph, ") and sequence alphabet (",
         seq.alph, ") do not match", call. = FALSE)

  ## --- PWM construction -------------------------------------------------
  motifs.pwm <- convert_type_internal(motifs, "PWM")
  ## Add a pseudocount to any motif that ended up with -Inf PWM entries,
  ## otherwise the downstream P-value CDF (which integerises scores) is
  ## undefined. Same approach as scan_sequences() (R/scan_sequences.R:253-260).
  needsfix <- vapply(motifs.pwm, function(x) any(is.infinite(x@motif)),
                     logical(1))
  if (any(needsfix)) {
    for (i in which(needsfix)) {
      motifs[[i]]     <- suppressMessages(normalize(motifs[[i]]))
      motifs.pwm[[i]] <- convert_type_internal(motifs[[i]], "PWM")
    }
  }
  score.mats <- lapply(motifs.pwm, function(x) x@motif)
  mot.names  <- vapply(motifs, function(x) x@name, character(1))
  if (any(duplicated(mot.names)))
    mot.names <- make.unique(mot.names)
  max.scores <- vapply(score.mats, function(m) sum(apply(m, 2, max)),
                       numeric(1))

  ## --- P-value -> per-motif score threshold -----------------------------
  thresholds <- motif_pvalue(motifs.pwm, pvalue = pvalue,
                             method = "dynamic", nthreads = nthreads)

  ## --- RC duplication ---------------------------------------------------
  strands     <- rep("+", length(score.mats))
  mot.indices <- seq_along(motifs)
  if (RC) {
    score.mats.rc <- lapply(score.mats,
      function(x) matrix(rev(as.numeric(x)), nrow = nrow(x)))
    score.mats  <- c(score.mats, score.mats.rc)
    thresholds  <- c(thresholds, thresholds)
    max.scores  <- c(max.scores, max.scores)
    strands     <- c(strands, rep("-", length(motifs)))
    mot.indices <- c(mot.indices, seq_along(motifs))
  }

  ## --- alphabet string for C++ ------------------------------------------
  alph <- switch(mot.alph, "DNA" = "ACGT", "RNA" = "ACGU")
  seq.names <- names(sequences)
  if (is.null(seq.names))
    seq.names <- sprintf("seq%d", seq_len(length(sequences)))

  ## --- scan -------------------------------------------------------------
  res <- scan_sequences_cpp(score.mats, as.character(sequences),
                            k = 1L, alph, thresholds,
                            nthreads = nthreads,
                            allow_nonfinite = FALSE, warnNA = FALSE)

  if (nrow(res) == 0L)
    return(empty_scan2_result(seq.names, sequences, return.granges))

  ## --- attach names, strand, original motif index -----------------------
  scan.idx        <- res$motif                  # 1..length(score.mats)
  res$motif       <- mot.names[mot.indices[scan.idx]]
  res$motif.i     <- mot.indices[scan.idx]
  res$sequence    <- seq.names[res$sequence]
  res$strand      <- strands[scan.idx]

  ## stop -> end (no swap on '-' strand)
  res$end  <- res$stop
  res$stop <- NULL

  ## --- per-hit pvalue (one batched motif_pvalue() call) -----------------
  ## Group scores by motif index and dispatch once; motif_pvalue() then makes
  ## a single C++ call parallelised across motifs. Only pass the motifs that
  ## actually got hits -- the batched C++ would otherwise build a CDF for
  ## every input motif, including ones with empty score vectors.
  hit_motif_i <- sort(unique(res$motif.i))
  hit_groups  <- split(res$score, factor(res$motif.i, levels = hit_motif_i))
  pvals_by_motif <- motif_pvalue(motifs.pwm[hit_motif_i],
                                 score    = hit_groups,
                                 method   = "dynamic",
                                 nthreads = nthreads)
  res$pvalue <- unsplit(pvals_by_motif,
                        factor(res$motif.i, levels = hit_motif_i))

  ## --- score.pct --------------------------------------------------------
  res$score.pct <- 100 * res$score / max.scores[scan.idx]

  ## --- reverse-complement match on '-' strand ---------------------------
  rev.rows <- res$strand == "-"
  if (any(rev.rows)) {
    rc <- if (mot.alph == "DNA")
            as.character(reverseComplement(DNAStringSet(res$match[rev.rows])))
          else
            as.character(reverseComplement(RNAStringSet(res$match[rev.rows])))
    res$match[rev.rows] <- rc
  }

  ## --- order & shape ----------------------------------------------------
  res <- res[order(res$motif.i, res$sequence.i, res$start), , drop = FALSE]
  rownames(res) <- NULL
  res <- res[, c("motif", "motif.i", "sequence", "sequence.i",
                 "start", "end", "strand",
                 "score", "score.pct", "match", "pvalue")]

  ## --- optional yamtk-style dedup --------------------------------------
  if (no.overlaps && nrow(res) > 0L) {
    res <- dedup_hits(res,
                      by            = no.overlaps.by,
                      ignore.strand = no.overlaps.by.strand,
                      ignore.motif  = no.overlaps.by.motif)
    rownames(res) <- NULL
  }

  build_scan2_result(res, sequences, return.granges)
}

# Decide whether to return a GRanges or a data.frame.
use_granges <- function(return.granges) {
  if (isTRUE(return.granges)) return(TRUE)
  if (isFALSE(return.granges)) return(FALSE)
  requireNamespace("GenomicRanges", quietly = TRUE)
}

build_scan2_result <- function(res, sequences, return.granges) {
  if (!use_granges(return.granges)) return(res)
  gr <- GenomicRanges::GRanges(
    seqnames = res$sequence,
    ranges   = IRanges::IRanges(start = res$start, end = res$end),
    strand   = res$strand,
    seqlengths = structure(width(sequences),
                           names = names(sequences) %||%
                             sprintf("seq%d", seq_len(length(sequences))))
  )
  S4Vectors::mcols(gr) <- res[, c("motif", "motif.i", "sequence.i",
                                  "score", "score.pct", "match", "pvalue"),
                              drop = FALSE]
  gr
}

empty_scan2_result <- function(seq.names, sequences, return.granges) {
  empty_df <- data.frame(
    motif       = character(0),
    motif.i     = integer(0),
    sequence    = character(0),
    sequence.i  = integer(0),
    start       = integer(0),
    end         = integer(0),
    strand      = character(0),
    score       = numeric(0),
    score.pct   = numeric(0),
    match       = character(0),
    pvalue      = numeric(0),
    stringsAsFactors = FALSE
  )
  if (!use_granges(return.granges)) return(empty_df)
  GenomicRanges::GRanges(
    seqnames = character(0),
    ranges   = IRanges::IRanges(),
    strand   = character(0),
    seqlengths = structure(width(sequences), names = seq.names)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x
