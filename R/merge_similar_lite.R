#' Cluster similar motifs by significance and merge each cluster.
#'
#' `merge_similar_lite()` is a faster minimalist counterpart to [merge_similar()]
#' that builds a pairwise similarity-significance graph using
#' [compare_motifs_lite()]'s empirical or parametric null q-values, finds
#' connected components on that graph, and collapses each component
#' (cluster) into a single merged motif via [merge_motifs_lite()].
#'
#' Unlike [merge_similar()] (which builds a distance matrix from raw
#' similarity scores and applies hierarchical clustering at an absolute
#' score cutoff) this function clusters on statistical
#' significance: two motifs are linked if their pairwise q-value
#' meets the user's `qvalue` cutoff. This matches the STAMP /
#' TOMTOM-clustering semantics ("group all motifs that are
#' significantly similar"), is deterministic, has no linkage-method
#' choice to expose, and lets the user reason about the cutoff in
#' q-value units rather than abstract distance.
#'
#' @param motifs `list` of motifs. DNA or RNA only.
#' @param qvalue `numeric(1)`. BH-adjusted q-value cutoff for linking
#'   two motifs in the significance graph. Default `0.01`.
#' @param null `character(1)`. `"empirical"` (default) or `"parametric"`,
#'   forwarded to [compare_motifs_lite()].
#' @param min.overlap `integer(1)`. Forwarded to [compare_motifs_lite()] /
#'   [merge_motifs_lite()]. Default `5L`.
#' @param RC `logical(1)`. If `TRUE` (default), test reverse-complement
#'   alignments.
#' @param bkg `numeric(4)` or `NULL`. Background frequencies, forwarded
#'   to [compare_motifs_lite()]. `NULL` uses uniform.
#' @param weighted `logical(1)`. Forwarded to [merge_motifs_lite()]. Default
#'   `FALSE`.
#' @param return.clusters `logical(1)`. If `TRUE`, return a
#'   `universalmotif_df` of cluster assignments instead of merged
#'   motifs. Default `FALSE`.
#' @param nthreads `numeric(1)`. Threads forwarded to [compare_motifs_lite()]
#'   and [merge_motifs_lite()]. `nthreads = 0` uses all available threads.
#'
#' @return If `return.clusters = FALSE` (default): a `list` of
#'   `universalmotif` S4 objects with similar motifs collapsed (the
#'   list is shorter than the input whenever any merging occurred).
#'   Singletons pass through unchanged.
#'
#'   If `return.clusters = TRUE`: a `universalmotif_df` (as produced
#'   by [to_df()]) describing the input motifs -- one row per input
#'   motif, with the `motif` column holding the motif objects and
#'   `name` their names -- augmented with two extra columns:
#'   `motif.i` (1-based input index) and `cluster` (1-based cluster
#'   id; motifs in the same cluster share an id).
#'
#' @details
#' The significance graph is built from a symmetrised q-value matrix:
#' `Qsym[i,j] = pmin(Q[i,j], Q[j,i])`. [compare_motifs_lite()]'s per-query BH
#' adjustment makes its native q-value matrix non-symmetric; the min-
#' symmetrisation defines a pair as linked if either direction meets
#' the cutoff (more permissive of the two, which is appropriate when
#' the goal is to group similar motifs).
#'
#' Clustering is by connected components on the resulting graph,
#' implemented with a small union-find. Each component (size >= 2)
#' becomes a single merged motif via [merge_motifs_lite()]; singletons
#' (size == 1) pass through unchanged.
#'
#' @examples
#' ## A clear cluster of 3 related motifs + 1 unrelated motif
#' m1 <- create_motif("TTGACATA", name = "a")
#' m2 <- create_motif("CTTGACAT", name = "b")
#' m3 <- create_motif("TGACATAT", name = "c")
#' m4 <- create_motif("GGGCCCCC", name = "unrelated")
#' merged <- merge_similar_lite(list(m1, m2, m3, m4), qvalue = 0.05)
#' length(merged)  # 2: one merged + the unrelated singleton
#'
#' @seealso [merge_similar()], [merge_motifs_lite()], [compare_motifs_lite()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @family lite motif functions
#' @export
merge_similar_lite <- function(motifs, qvalue = 0.01,
                           null = c("empirical", "parametric"),
                           min.overlap = 5L, RC = TRUE,
                           bkg = NULL, weighted = FALSE,
                           return.clusters = FALSE,
                           nthreads = 1) {

  ## --- arg validation ---------------------------------------------------
  if (missing(motifs))
    stop("`motifs` is required", call. = FALSE)
  null <- match.arg(null)
  if (!is.numeric(qvalue) || length(qvalue) != 1L || is.na(qvalue) ||
      qvalue <= 0 || qvalue > 1)
    stop("`qvalue` must be a single numeric in (0, 1]", call. = FALSE)
  if (!isTRUEorFALSE(RC))
    stop("`RC` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(weighted))
    stop("`weighted` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(return.clusters))
    stop("`return.clusters` must be a single logical", call. = FALSE)

  nthreads <- resolve_nthreads(nthreads)

  ## --- early empty-input handling (before convert_motifs, which would error) ---
  if (is.list(motifs) && length(motifs) == 0L) {
    if (return.clusters) {
      empty <- to_df(list(create_motif()))[0, , drop = FALSE]
      empty$motif.i <- integer(0)
      empty$cluster <- integer(0)
      return(empty)
    }
    return(list())
  }

  ## --- normalise input --------------------------------------------------
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(alphs)) != 1L)
    stop("all motifs must share the same alphabet", call. = FALSE)
  mot.alph <- unique(alphs)
  if (!mot.alph %in% c("DNA", "RNA"))
    stop("`merge_similar_lite()` only supports DNA/RNA motifs; got `",
         mot.alph, "`. Use `merge_similar()` for other alphabets.",
         call. = FALSE)

  n_mot <- length(motifs)
  mot.names <- vapply(motifs, function(x) x@name, character(1))

  ## --- single-motif fast path ------------------------------------------
  if (n_mot == 1L) {
    if (return.clusters) {
      return(clusters_to_df(motifs, 1L))
    }
    return(motifs)
  }

  ## --- pairwise q-values via compare_motifs_lite ---------------------------
  Q <- suppressWarnings(compare_motifs_lite(motifs,
                                        qvalue      = 1,
                                        min.overlap = min.overlap,
                                        RC          = RC,
                                        null        = null,
                                        bkg         = bkg,
                                        matrix.out  = "qvalue",
                                        nthreads    = nthreads))
  ## Symmetrise: pair linked if either direction is significant.
  Qsym <- pmin(Q, t(Q), na.rm = FALSE)

  ## --- union-find on significance edges --------------------------------
  parent <- seq_len(n_mot)
  find <- function(x) {
    ## Iterative path compression.
    root <- x
    while (parent[root] != root) root <- parent[root]
    while (parent[x] != root) { nxt <- parent[x]; parent[x] <<- root; x <- nxt }
    root
  }
  unite <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if (ra != rb) parent[rb] <<- ra
  }

  for (i in seq_len(n_mot - 1L)) {
    for (j in seq.int(i + 1L, n_mot)) {
      q_ij <- Qsym[i, j]
      if (!is.na(q_ij) && q_ij <= qvalue) unite(i, j)
    }
  }

  ## --- collect clusters in stable order --------------------------------
  roots    <- vapply(seq_len(n_mot), find, integer(1))
  ## Renumber clusters by first-occurrence order so the output is
  ## independent of root-selection details.
  cluster_id <- integer(n_mot)
  seen <- integer(0)
  for (i in seq_len(n_mot)) {
    r <- roots[i]
    pos <- match(r, seen)
    if (is.na(pos)) {
      seen <- c(seen, r)
      cluster_id[i] <- length(seen)
    } else {
      cluster_id[i] <- pos
    }
  }

  if (return.clusters) {
    return(clusters_to_df(motifs, cluster_id))
  }

  ## --- merge each multi-motif cluster ----------------------------------
  out <- vector("list", max(cluster_id))
  for (cid in seq_len(max(cluster_id))) {
    members <- which(cluster_id == cid)
    if (length(members) == 1L) {
      out[[cid]] <- motifs[[members]]
    } else {
      out[[cid]] <- merge_motifs_lite(motifs[members],
                                  min.overlap = min.overlap,
                                  RC          = RC,
                                  weighted    = weighted,
                                  nthreads    = nthreads)
    }
  }
  out
}

## Build the return.clusters output: a universalmotif_df (via to_df())
## of the input motifs, augmented with `motif.i` (1-based input index)
## and `cluster` (1-based cluster id) columns. The `motif` column holds
## the input motif objects, `name` their names, etc.
clusters_to_df <- function(motifs, cluster_id) {
  df <- to_df(motifs)
  df$motif.i <- seq_len(nrow(df))
  df$cluster <- as.integer(cluster_id)
  df
}
