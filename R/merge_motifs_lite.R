#' Merge a list of motifs into a single consensus motif.
#'
#' `merge_motifs_lite()` is a faster minimalist counterpart to [merge_motifs()] that
#' aligns motifs using `compare_motifs_lite()`'s alignment finder and averages
#' their position-probability columns in a single shared coordinate frame.
#'
#' @param motifs `list` of motifs. See [convert_motifs()] for accepted
#'   formats. DNA or RNA only.
#' @param min.overlap `integer(1)`. Minimum overlap (in columns) for an
#'   alignment to be considered. Default `5L`. Forwarded to
#'   [compare_motifs_lite()].
#' @param RC `logical(1)`. If `TRUE` (default), also test reverse-complement
#'   alignments. Motifs whose best alignment hits the `-` strand of the
#'   anchor are reverse-complemented before averaging, so the merged
#'   motif is always returned in the anchor's orientation.
#' @param new.name `character(1)` or `NULL`. Name for the merged motif.
#'   If `NULL` (default), set to `paste0("merged_", paste(names, collapse = "+"))`.
#' @param weighted `logical(1)`. If `TRUE`, weight each input motif's
#'   column contribution by its `@nsites` slot. Default `FALSE` (simple
#'   unweighted mean). Motifs whose `@nsites` is missing or zero fall back
#'   to weight 1 in either mode.
#' @param nthreads `numeric(1)`. Number of threads passed to
#'   [compare_motifs_lite()]. `nthreads = 0` uses all available threads.
#'
#' @return A single `universalmotif` S4 object of type `"PPM"`.
#'
#' @details
#' For each non-anchor motif, [compare_motifs_lite()] gives the best
#' alignment offset, strand, and overlap against the anchor. Each
#' aligned motif (reverse-complemented if its best alignment is `-`
#' strand) is then placed in a unified column coordinate frame whose
#' origin coincides with the anchor's first column. The output's column
#' at frame position `p` is the (optionally `@nsites`-weighted) mean of
#' the PPM columns from all input motifs that cover `p`; positions
#' covered by only some motifs are still emitted (no automatic
#' IC-trimming). Run [trim_motifs()] on the result if you want to
#' shrink low-information edges.
#'
#' The merged motif's `@bkg` is the unweighted mean of the input
#' motifs' `@bkg` vectors (or uniform if missing). `@nsites` is the
#' sum of input `@nsites`.
#'
#' @examples
#' library(universalmotif)
#' m1 <- create_motif("TTGACATA", name = "a")
#' m2 <- create_motif("CTTGACAT", name = "b")
#' m3 <- create_motif("TGACATAT", name = "c")
#' merged <- merge_motifs_lite(list(m1, m2, m3))
#' merged
#'
#' @seealso [merge_motifs()], [merge_similar_lite()], [compare_motifs_lite()],
#'   [trim_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @family lite motif functions
#' @export
merge_motifs_lite <- function(motifs, min.overlap = 5L, RC = TRUE,
                          new.name = NULL, weighted = FALSE,
                          nthreads = 1) {

  ## --- arg validation ---------------------------------------------------
  if (missing(motifs))
    stop("`motifs` is required", call. = FALSE)
  if (!is.numeric(min.overlap) || length(min.overlap) != 1L ||
      min.overlap < 1L)
    stop("`min.overlap` must be a positive integer", call. = FALSE)
  if (!isTRUEorFALSE(RC))
    stop("`RC` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(weighted))
    stop("`weighted` must be a single logical", call. = FALSE)
  if (!is.null(new.name) &&
      (!is.character(new.name) || length(new.name) != 1L))
    stop("`new.name` must be NULL or a single string", call. = FALSE)

  nthreads <- resolve_nthreads(nthreads)

  ## --- normalise input --------------------------------------------------
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  if (length(motifs) == 0L)
    stop("`motifs` must contain at least one motif", call. = FALSE)

  alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(alphs)) != 1L)
    stop("all motifs must share the same alphabet", call. = FALSE)
  mot.alph <- unique(alphs)
  if (!mot.alph %in% c("DNA", "RNA"))
    stop("`merge_motifs_lite()` only supports DNA/RNA motifs; got `",
         mot.alph, "`. Use `merge_motifs()` for other alphabets.",
         call. = FALSE)

  ## --- single-motif fast path ------------------------------------------
  if (length(motifs) == 1L) {
    out <- convert_type_internal(motifs[[1L]], "PPM")
    if (!is.null(new.name)) out@name <- new.name
    return(out)
  }

  motifs.ppm <- convert_type_internal(motifs, "PPM")
  mot.names  <- vapply(motifs.ppm, function(x) x@name, character(1))
  mot.mats   <- lapply(motifs.ppm, function(x) x@motif)
  widths     <- vapply(mot.mats, ncol, integer(1))

  ## per-motif nsites (fallback = 1 when missing / 0 / NA)
  nsites <- vapply(motifs.ppm, function(x) {
    n <- x@nsites
    if (length(n) == 1L && !is.na(n) && n > 0) n else 1
  }, numeric(1))

  ## --- choose anchor (highest IC) --------------------------------------
  ic <- vapply(motifs.ppm, function(x) x@icscore, numeric(1))
  anchor <- which.max(ic)

  ## --- align every motif to the anchor in ONE C++ call -----------------
  qi <- rep.int(anchor, length(motifs.ppm))
  ti <- seq_along(motifs.ppm)
  al <- compare_motifs2_align_cpp(mot.mats,
                                  qi          = as.integer(qi),
                                  ti          = as.integer(ti),
                                  min_overlap = as.integer(min.overlap),
                                  RC          = as.logical(RC),
                                  nthreads    = as.integer(nthreads))

  ## --- compute frame placement per motif --------------------------------
  ## Anchor sits at frame_start = 0, frame_end = w_anchor.
  ## For each non-anchor motif i (potentially RC'd if strand=="-"), the
  ## aligned region maps anchor cols [qstart_orig, qstart_orig+L) to
  ## (oriented) motif cols [t_start_oriented, t_start_oriented+L).
  ## frame_offset = qstart_orig - t_start_oriented places motif col 0 at
  ## frame position (-t_start_oriented + qstart_orig). Motifs that don't
  ## meet min.overlap (overlap == 0) get a degenerate placement and are
  ## dropped from the average.
  w_anchor      <- widths[anchor]
  oriented_mats <- vector("list", length(motifs.ppm))
  frame_offset  <- integer(length(motifs.ppm))
  used          <- logical(length(motifs.ppm))

  for (i in seq_along(motifs.ppm)) {
    if (i == anchor) {
      oriented_mats[[i]] <- mot.mats[[i]]
      frame_offset[i]    <- 0L
      used[i]            <- TRUE
      next
    }
    L_i <- al$overlap[i]
    if (is.na(L_i) || L_i < min.overlap) {
      used[i] <- FALSE
      next
    }
    qstart_oriented <- al$q_start_oriented[i]
    t_start         <- max(0L, -al$offset[i])
    strand_neg      <- isTRUE(al$strand[i] == 1L)
    w_m             <- widths[i]

    if (strand_neg) {
      ## RC the motif. Compare_motifs2 reports q_start_oriented in
      ## RC-anchor coordinates; the corresponding original-anchor start
      ## is (w_anchor - q_start_oriented - L).
      oriented_mats[[i]] <- rev_comp_mat(mot.mats[[i]])
      qstart_orig        <- w_anchor - qstart_oriented - L_i
      ## After RC'ing the target, its t_start in the RC'd orientation
      ## flips too: new t_start = w_m - t_start - L.
      t_start_oriented   <- w_m - t_start - L_i
      frame_offset[i]    <- qstart_orig - t_start_oriented
    } else {
      oriented_mats[[i]] <- mot.mats[[i]]
      frame_offset[i]    <- qstart_oriented - t_start
    }
    used[i] <- TRUE
  }

  ## Drop motifs that failed alignment.
  if (!any(used))
    stop("no motif could be aligned to the anchor at min.overlap = ",
         min.overlap, call. = FALSE)
  use.idx <- which(used)

  ## --- unified frame extent --------------------------------------------
  frame_starts <- frame_offset[use.idx]
  frame_ends   <- frame_offset[use.idx] +
                  vapply(oriented_mats[use.idx], ncol, integer(1))
  frame_lo <- min(frame_starts)
  frame_hi <- max(frame_ends)
  out_w    <- frame_hi - frame_lo

  ## --- per-position average --------------------------------------------
  alph_letters <- if (mot.alph == "DNA") c("A","C","G","T")
                  else                   c("A","C","G","U")
  out_mat <- matrix(0.0, nrow = 4L, ncol = out_w,
                    dimnames = list(alph_letters, NULL))
  weights <- if (weighted) nsites[use.idx] else rep.int(1, length(use.idx))

  for (p in seq_len(out_w)) {
    frame_p <- frame_lo + p - 1L  # 0-based frame position
    contrib_sum <- numeric(4L)
    weight_sum  <- 0
    for (k in seq_along(use.idx)) {
      i <- use.idx[k]
      col_local <- frame_p - frame_offset[i]  # 0-based motif col
      if (col_local < 0L || col_local >= ncol(oriented_mats[[i]])) next
      contrib_sum <- contrib_sum + weights[k] * oriented_mats[[i]][, col_local + 1L]
      weight_sum  <- weight_sum + weights[k]
    }
    if (weight_sum == 0) {
      out_mat[, p] <- 0.25
    } else {
      out_mat[, p] <- contrib_sum / weight_sum
    }
  }

  ## --- aggregated metadata ---------------------------------------------
  total_nsites <- sum(nsites[use.idx])
  ## Average bkg vector (in alphabet letter order)
  bkgs <- lapply(motifs.ppm[use.idx], function(x) {
    b <- x@bkg[alph_letters]
    if (any(is.na(b))) rep.int(0.25, 4L) else as.numeric(b)
  })
  bkg_avg <- Reduce(`+`, bkgs) / length(bkgs)
  names(bkg_avg) <- alph_letters

  if (is.null(new.name)) {
    new.name <- paste0("merged_",
                       paste(mot.names[use.idx], collapse = "+"))
  }

  create_motif(out_mat,
               alphabet    = mot.alph,
               type        = "PPM",
               name        = new.name,
               nsites      = total_nsites,
               bkg         = bkg_avg,
               pseudocount = motifs.ppm[[anchor]]@pseudocount)
}

## ---------------------------------------------------------------------------
## Helper: column reverse-complement of a 4xN DNA/RNA PPM. Rownames must be in
## {A,C,G,T} or {A,C,G,U} order (the alphabetical default produced by
## convert_motifs() / create_motif()).
## ---------------------------------------------------------------------------
rev_comp_mat <- function(mat) {
  ## Row 1 = A, 2 = C, 3 = G, 4 = T/U.
  ## Reverse columns, then swap (A,T) and (C,G) by permuting rows. The
  ## `m[c(4,3,2,1), ]` step also carries the input row *names* over in
  ## the new order (T,G,C,A), so for callers that read the result by
  ## row name (e.g. the view_motifs() logo renderer) we restore the
  ## alphabetical A,C,G,T order. Callers that read positionally (e.g.
  ## merge_motifs_lite()'s per-position accumulator) are unaffected because
  ## the underlying data placement is unchanged.
  rn <- rownames(mat)
  m <- mat[, rev(seq_len(ncol(mat))), drop = FALSE]
  m <- m[c(4L, 3L, 2L, 1L), , drop = FALSE]
  rownames(m) <- rn
  m
}
