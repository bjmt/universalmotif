#' Implant sampled motif instances into sequences at known positions.
#'
#' [implant_motifs()] overwrites bases in an existing `XStringSet` with
#' samples drawn column-by-column from one or more motif PPMs at chosen
#' positions, returning the modified sequences together (optionally)
#' with a ground-truth data.frame of where each implant was placed.
#' This produces benchmark-grade positive sequences for testing
#' [scan_sequences2()], [motif_finder()], [enrich_motifs2()],
#' [motif_peaks()], or any other discovery / scanning pipeline against
#' a known answer.
#'
#' Three insertion modes (mutually exclusive; later non-`NULL` argument
#' wins):
#'
#'   1. Fixed count (`n.per.seq`): plant exactly `n.per.seq`
#'      instances in every sequence. Matches HOMER-style ground-truth
#'      simulations.
#'   2. Poisson rate (`rate`): per-sequence count drawn from
#'      `rpois(1, rate * width(seq))`. Density scales naturally with
#'      sequence length.
#'   3. Explicit positions (`positions`): a `list` of length
#'      `length(sequences)`; element `i` is a 1-based integer vector of
#'      start positions for that sequence. Bypasses `centre.bias`,
#'      `min.spacing`, and `max.retries`.
#'
#' For each implant the function (a) picks a motif uniformly at random
#' from `motifs`, (b) chooses a strand (uniform `+`/`-` if
#' `strand = "both"`), (c) picks a start position (uniform or centre-
#' biased via Irwin-Hall, retried up to `max.retries` times if it
#' collides with an already-placed implant within `min.spacing`),
#' (d) samples one instance from the motif's PPM column-by-column,
#' reverse-complements if needed, and (e) overwrites that span in the
#' target sequence.
#'
#' @param motifs See [convert_motifs()] for accepted motif formats.
#'   DNA or RNA only, amino-acid and custom alphabet motifs are
#'   rejected. A single motif or a list of motifs.
#' @param sequences `XStringSet`. Host sequences to plant into. Must
#'   share the motif alphabet.
#' @param n.per.seq `integer(1)`. Mode A. Number of instances to plant
#'   in each sequence. Default `1L`.
#' @param rate `numeric(1)` or `NULL`. Mode B. Per-bp Poisson rate.
#'   If non-`NULL`, overrides `n.per.seq`. Default `NULL`.
#' @param positions `list` or `NULL`. Mode C. List of length
#'   `length(sequences)`; element `i` is an integer vector of 1-based
#'   start positions. Overrides `n.per.seq` and `rate`. Default `NULL`.
#' @param strand `"both"` or `"+"`. Whether to also draw from the
#'   reverse complement. Default `"both"`.
#' @param centre.bias `integer(1)`. Irwin-Hall N. `1L` (default) is
#'   uniform; larger values cluster implants near the centre of each
#'   sequence. Ignored in `positions` mode.
#' @param min.spacing `integer(1)`. Minimum gap (in bp) between the
#'   end of one implant and the start of the next within the same
#'   sequence. `0L` (default) allows abutting but not overlapping
#'   implants. Ignored in `positions` mode.
#' @param max.retries `integer(1)`. Per-implant cap on placement
#'   attempts before giving up (with a warning at the end if any
#'   targets fell short). Default `100L`.
#' @param return.indices `logical(1)`. If `FALSE` (default), return
#'   the modified `XStringSet`. If `TRUE`, return a `data.frame` of
#'   per-implant ground-truth metadata instead.
#'
#' @return If `return.indices = FALSE`: an `XStringSet` of the same
#'   length, width, and class as `sequences`, with bases overwritten
#'   at the implant positions. Names are preserved. If
#'   `return.indices = TRUE`: a `data.frame` with columns
#'   `sequence.i` (1-based), `motif.i` (1-based index into `motifs`),
#'   `start` (1-based), `width`, `strand` (`"+"` / `"-"`), and
#'   `planted` (the actual implanted string, already reverse-
#'   complemented when `strand == "-"`).
#'
#' @examples
#' \dontrun{
#' library(universalmotif)
#' data(ArabidopsisMotif)
#' seqs  <- create_sequences("DNA", seqnum = 100, seqlen = 500)
#' set.seed(1)
#' out   <- implant_motifs(ArabidopsisMotif, seqs, n.per.seq = 1)
#' set.seed(1)
#' truth <- implant_motifs(ArabidopsisMotif, seqs, n.per.seq = 1,
#'                         return.indices = TRUE)
#' ## Recover-rate check
#' hits <- scan_sequences(ArabidopsisMotif, out, threshold = 0.85,
#'                        threshold.type = "logodds")
#' ## Density check
#' peaks <- motif_peaks(out, motif = ArabidopsisMotif)
#' plot_motif_peaks(peaks)
#' }
#'
#' @seealso [sample_sites()], [create_sequences()], [scan_sequences2()],
#'   [motif_finder()], [motif_peaks()], [enrich_motifs2()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
implant_motifs <- function(motifs, sequences,
                           n.per.seq      = 1L,
                           rate           = NULL,
                           positions      = NULL,
                           strand         = c("both", "+"),
                           centre.bias    = 1L,
                           min.spacing    = 0L,
                           max.retries    = 100L,
                           return.indices = FALSE) {

  ## --- arg validation --------------------------------------------------
  if (missing(motifs) || missing(sequences))
    stop("`motifs` and `sequences` are both required", call. = FALSE)
  if (!is(sequences, "XStringSet"))
    stop("`sequences` must be an XStringSet object", call. = FALSE)
  if (!isTRUEorFALSE(return.indices))
    stop("`return.indices` must be a single logical", call. = FALSE)
  strand <- match.arg(strand)
  if (!is.numeric(centre.bias) || length(centre.bias) != 1L ||
      centre.bias < 1L)
    stop("`centre.bias` must be a positive integer (1 = uniform)",
         call. = FALSE)
  if (!is.numeric(min.spacing) || length(min.spacing) != 1L ||
      min.spacing < 0L)
    stop("`min.spacing` must be a non-negative integer", call. = FALSE)
  if (!is.numeric(max.retries) || length(max.retries) != 1L ||
      max.retries < 1L)
    stop("`max.retries` must be a positive integer", call. = FALSE)

  ## --- motif normalisation + alphabet checks ---------------------------
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  if (length(motifs) == 0L)
    stop("`motifs` must contain at least one motif", call. = FALSE)
  motifs <- lapply(motifs, convert_type_internal, "PPM")
  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) != 1L)
    stop("all motifs must share the same alphabet", call. = FALSE)
  mot.alph <- unique(mot.alphs)
  if (!mot.alph %in% c("DNA", "RNA"))
    stop("`implant_motifs()` only supports DNA/RNA motifs; got `",
         mot.alph, "`.", call. = FALSE)
  seq.alph <- seqtype(sequences)
  if (seq.alph != mot.alph)
    stop("motif alphabet (", mot.alph, ") and sequence alphabet (",
         seq.alph, ") do not match", call. = FALSE)
  if (strand == "both" && mot.alph != "DNA" && mot.alph != "RNA")
    stop("strand = \"both\" requires DNA/RNA motifs", call. = FALSE)

  ## --- pick mode -------------------------------------------------------
  mode <- "count"
  if (!is.null(rate))      mode <- "rate"
  if (!is.null(positions)) mode <- "positions"

  if (mode == "count") {
    if (!is.numeric(n.per.seq) || length(n.per.seq) != 1L || n.per.seq < 0L)
      stop("`n.per.seq` must be a non-negative integer", call. = FALSE)
  } else if (mode == "rate") {
    if (!is.numeric(rate) || length(rate) != 1L || rate < 0)
      stop("`rate` must be a non-negative numeric", call. = FALSE)
  } else {
    if (!is.list(positions) || length(positions) != length(sequences))
      stop("`positions` must be a list of length equal to `sequences`",
           call. = FALSE)
  }

  ## --- precompute per-motif PPMs and widths ----------------------------
  mot.ppms   <- lapply(motifs, function(x) x@motif)
  mot.widths <- vapply(motifs, function(x) ncol(x@motif), integer(1))
  mot.alph.l <- rownames(mot.ppms[[1]])   # alphabet letter order

  seq.widths <- as.integer(width(sequences))
  n.seq      <- length(sequences)

  if (any(seq.widths < min(mot.widths)))
    warning("some sequences are shorter than the narrowest motif; ",
            "no implants will be placed in them", call. = FALSE)

  ## --- main loop -------------------------------------------------------
  ## RNG comes from the caller's global state (set.seed() before the
  ## call for reproducibility).
  seq.chars  <- as.character(sequences)   # mutable buffer (R strings copy-on-mod)
  meta.list  <- vector("list", n.seq)
  short_capacity <- FALSE

  for (i in seq_len(n.seq)) {
    L <- seq.widths[i]
    ## decide how many to place
    n.i <- switch(mode,
      count     = as.integer(n.per.seq),
      rate      = stats::rpois(1L, rate * L),
      positions = length(positions[[i]])
    )
    if (n.i == 0L) next

    placed.starts <- integer(0)
    placed.ends   <- integer(0)
    placed.mot    <- integer(0)
    placed.strand <- character(0)
    placed.str    <- character(0)

    for (j in seq_len(n.i)) {
      ## sample motif identity once per implant
      m.idx <- if (length(motifs) == 1L) 1L
               else                       sample.int(length(motifs), 1L)
      mw    <- mot.widths[m.idx]
      if (mw > L) next   # narrowest motif doesn't fit; skip silently

      ## strand
      s <- if (strand == "+") "+"
           else if (sample.int(2L, 1L) == 1L) "+" else "-"

      ## start position
      if (mode == "positions") {
        start.j <- as.integer(positions[[i]][j])
        if (is.na(start.j) || start.j < 1L || start.j + mw - 1L > L) {
          warning("dropping out-of-bounds position for sequence ", i,
                  " (start=", start.j, ", width=", mw, ", seq.len=", L,
                  ")", call. = FALSE)
          next
        }
        end.j <- start.j + mw - 1L
        if (length(placed.starts) &&
            collides(start.j, end.j, placed.starts, placed.ends, 0L)) {
          warning("position ", start.j, " in sequence ", i,
                  " overlaps a previous implant; planting anyway",
                  call. = FALSE)
        }
      } else {
        ## random / centre-biased with collision retry
        retries <- 0L
        repeat {
          retries <- retries + 1L
          u <- if (centre.bias == 1L) stats::runif(1L)
               else                    mean(stats::runif(as.integer(centre.bias)))
          start.j <- 1L + as.integer(floor(u * (L - mw + 1L)))
          if (start.j < 1L) start.j <- 1L
          if (start.j > L - mw + 1L) start.j <- L - mw + 1L
          end.j   <- start.j + mw - 1L
          if (!collides(start.j, end.j, placed.starts, placed.ends,
                        as.integer(min.spacing))) break
          if (retries >= as.integer(max.retries)) {
            short_capacity <- TRUE
            start.j <- NA_integer_
            break
          }
        }
        if (is.na(start.j)) next
      }

      ## sample one instance from the PPM
      planted <- sample_one_instance(mot.ppms[[m.idx]], mot.alph.l)
      if (s == "-") planted <- revcomp_string(planted, mot.alph)

      ## overwrite in the sequence buffer
      substr(seq.chars[i], start.j, end.j) <- planted

      placed.starts <- c(placed.starts, start.j)
      placed.ends   <- c(placed.ends,   end.j)
      placed.mot    <- c(placed.mot,    m.idx)
      placed.strand <- c(placed.strand, s)
      placed.str    <- c(placed.str,    planted)
    }

    if (length(placed.starts) > 0L)
      meta.list[[i]] <- data.frame(
        sequence.i = i,
        motif.i    = placed.mot,
        start      = placed.starts,
        width      = placed.ends - placed.starts + 1L,
        strand     = placed.strand,
        planted    = placed.str,
        stringsAsFactors = FALSE
      )
  }

  if (short_capacity)
    warning("one or more implants could not be placed within `max.retries` ",
            "attempts; consider increasing sequence length, lowering ",
            "`n.per.seq`/`rate`, or reducing `min.spacing`", call. = FALSE)

  if (return.indices) {
    meta <- do.call(rbind, meta.list)
    if (is.null(meta))
      meta <- data.frame(sequence.i = integer(0), motif.i = integer(0),
                         start = integer(0), width = integer(0),
                         strand = character(0), planted = character(0),
                         stringsAsFactors = FALSE)
    rownames(meta) <- NULL
    return(meta)
  }

  ## rebuild XStringSet, preserving names and subclass
  out <- switch(mot.alph,
                DNA = Biostrings::DNAStringSet(seq.chars),
                RNA = Biostrings::RNAStringSet(seq.chars))
  names(out) <- names(sequences)
  out
}

## ---- internal helpers ---------------------------------------------------

## Sample one motif instance: per-column draw from the PPM probabilities.
## Returns a single character string of width = ncol(mat).
sample_one_instance <- function(mat, alph) {
  collapse_cpp(vapply(seq_len(ncol(mat)),
                      function(j) sample(alph, 1L, prob = mat[, j]),
                      character(1)))
}

## Reverse-complement a single string for DNA or RNA.
revcomp_string <- function(s, alph) {
  x <- if (alph == "DNA") Biostrings::DNAString(s)
       else                Biostrings::RNAString(s)
  as.character(Biostrings::reverseComplement(x))
}

## Linear scan for interval collision with min-spacing `gap`.
## (Linear is fine for the typical few-dozen implants per sequence.)
collides <- function(s, e, ps, pe, gap) {
  if (length(ps) == 0L) return(FALSE)
  any((s <= pe + gap) & (e + gap >= ps))
}
