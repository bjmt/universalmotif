#' Deduplicate overlapping motif hits.
#'
#' Given a set of motif hits, typically the output of
#' [scan_sequences()] or [scan_sequences_lite()], or a `GRanges` /
#' `data.frame` from any other source, drop overlapping hits within
#' each `(sequence, motif, strand)` group using a greedy
#' algorithm: cluster overlapping hits, then within each cluster keep
#' the hit with the best priority and, after that, any other hit that
#' does not overlap an already-kept hit.
#'
#' Note that this differs from [scan_sequences()]'s built-in
#' `no.overlaps` argument, which collapses every connected-overlap
#' cluster to a single survivor (hierarchical clustering with
#' `cutree`). The greedy approach used here keeps legitimately distinct
#' hits even when an intermediate noisy hit bridges them, which is
#' usually what users want.
#'
#' Coordinate convention: rows must have `start <= end`. Overlaps are
#' inclusive: two hits overlap if `a$start <= b$end & b$start <= a$end`.
#'
#' @param hits `data.frame` or `GRanges`. For a `data.frame`, must have
#'   columns `start`, `end`, and the columns named by `by`, `seq`,
#'   `motif`, `strand`. For a `GRanges`, the `seqnames`/`start`/`end`/
#'   `strand` slots are used; the `motif` and priority columns are read
#'   from `mcols()`.
#' @param by `character(1)`. Name of the priority column. The semantic
#'   is auto-detected: a column named `pvalue` (or anything starting
#'   with `pval`) is treated as lower-is-better; anything else (e.g.
#'   `score`) is treated as higher-is-better. Use `reverse = TRUE` to
#'   flip.
#' @param reverse `logical(1)`. Flip the priority direction.
#' @param ignore.strand `logical(1)`. If `TRUE`, `+` and `-` hits
#'   compete with each other at overlapping coordinates.
#' @param ignore.motif `logical(1)`. If `TRUE`, hits from different
#'   motifs compete with each other at overlapping coordinates.
#' @param seq `character(1)`. Name of the sequence-id column in a
#'   `data.frame` input. Default `"sequence"`. Ignored for `GRanges`.
#' @param motif `character(1)`. Name of the motif-id column in a
#'   `data.frame` input. Default `"motif"`. Ignored for `GRanges`.
#' @param strand `character(1)`. Name of the strand column in a
#'   `data.frame` input. Default `"strand"`. Ignored for `GRanges`.
#'
#' @return The input `hits`, subsetted to the kept rows in the original
#'   order.
#'
#' @examples
#' library(universalmotif)
#' motifs <- list(create_motif("TATAAA", name = "M1"),
#'                create_motif("CACGTG", name = "M2"))
#' seqs <- create_sequences(seqnum = 5, seqlen = 400, rng.seed = 1)
#' hits <- scan_sequences_lite(motifs, seqs, pvalue = 5e-2,
#'                         return.granges = FALSE)
#' deduped <- dedup_hits(hits)
#' nrow(hits)
#' nrow(deduped)
#'
#' @seealso [scan_sequences_lite()], [scan_sequences()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
dedup_hits <- function(hits,
                       by = "pvalue",
                       reverse = FALSE,
                       ignore.strand = FALSE,
                       ignore.motif  = FALSE,
                       seq    = "sequence",
                       motif  = "motif",
                       strand = "strand") {

  if (missing(hits)) stop("`hits` is required", call. = FALSE)
  if (!is.character(by) || length(by) != 1L)
    stop("`by` must be a single column name", call. = FALSE)
  for (arg_name in c("reverse", "ignore.strand", "ignore.motif")) {
    if (!isTRUEorFALSE(get(arg_name)))
      stop(sprintf("`%s` must be a single logical", arg_name),
           call. = FALSE)
  }

  is_gr <- inherits(hits, "GRanges")

  if (is_gr) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE))
      stop("GRanges input requires the GenomicRanges package", call. = FALSE)
    mc       <- GenomicRanges::mcols(hits)
    seq_v    <- as.character(GenomicRanges::seqnames(hits))
    start_v  <- as.integer(GenomicRanges::start(hits))
    end_v    <- as.integer(GenomicRanges::end(hits))
    strand_v <- as.character(GenomicRanges::strand(hits))
    if (!by %in% colnames(mc))
      stop(sprintf("priority column `%s` not found in mcols(hits)", by),
           call. = FALSE)
    priority_v <- as.numeric(mc[[by]])
    motif_v <- if ("motif" %in% colnames(mc))
                 as.character(mc[["motif"]])
               else
                 rep("motif", length(hits))
  } else {
    if (!is.data.frame(hits))
      stop("`hits` must be a data.frame or GRanges", call. = FALSE)
    needed <- c(seq, motif, strand, "start", "end", by)
    miss <- setdiff(needed, colnames(hits))
    if (length(miss))
      stop("`hits` is missing column(s): ",
           paste(miss, collapse = ", "), call. = FALSE)
    seq_v      <- as.character(hits[[seq]])
    motif_v    <- as.character(hits[[motif]])
    strand_v   <- as.character(hits[[strand]])
    start_v    <- as.integer(hits$start)
    end_v      <- as.integer(hits$end)
    priority_v <- as.numeric(hits[[by]])
  }

  if (anyNA(start_v) || anyNA(end_v))
    stop("`start` and `end` cannot contain NA", call. = FALSE)
  if (any(start_v > end_v))
    stop("all rows must satisfy start <= end", call. = FALSE)
  if (anyNA(priority_v))
    stop(sprintf("priority column `%s` cannot contain NA", by),
         call. = FALSE)

  # Direction: p-value-like columns are lower-is-better; everything else
  # is higher-is-better. Caller can flip via `reverse`.
  lower_is_better <- grepl("^pval", by, ignore.case = TRUE)
  if (reverse) lower_is_better <- !lower_is_better
  if (lower_is_better) priority_v <- -priority_v

  # Encode strings as 0-based integer codes for the C++ side.
  seq_codes    <- match(seq_v,    unique(seq_v))    - 1L
  motif_codes  <- match(motif_v,  unique(motif_v))  - 1L
  strand_codes <- match(strand_v, c("+", "-", "*")) - 1L
  strand_codes[is.na(strand_codes)] <- -1L  # treat other values as a separate group

  keep <- dedup_hits_cpp(seq_codes, motif_codes, start_v, end_v,
                         strand_codes, priority_v,
                         ignore_strand = ignore.strand,
                         ignore_motif  = ignore.motif)

  if (is_gr) hits[keep] else hits[keep, , drop = FALSE]
}
