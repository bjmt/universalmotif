#' Trim CWM motifs by absolute column sum.
#'
#' Trim the low-contribution edge columns of one or more
#' Contribution Weight Matrix (CWM) motifs. CWMs do not have a
#' probabilistic interpretation, so [trim_motifs()]'s
#' information-content threshold is not the natural cutoff;
#' instead, [trim_cwm()] drops edge columns whose absolute
#' column sum (`sum_i |cwm[i,j]|`) falls below a threshold.
#'
#' Two threshold modes are supported. By default, `trim_cwm()`
#' uses the TF-MoDISco-lite *fraction-of-peak* rule: drop edge
#' columns whose absolute column sum is below
#' `trim.threshold * max_j(sum_i |cwm[i,j]|)`, with
#' `trim.threshold = 0.3` matching TF-MoDISco-lite's
#' `--trim_threshold` default. When `abs.threshold` is supplied
#' (non-`NULL`), the function switches to an absolute-value cutoff:
#' edge columns with absolute column sum below `abs.threshold` are
#' dropped, ignoring `trim.threshold`.
#'
#' In both modes the trimmer walks inward only from the chosen
#' edges (`trim.from`), matching [trim_motifs()]'s edge-only
#' semantics; interior low-contribution columns are kept.
#'
#' @param motifs A CWM motif or list of CWM motifs (`type = "CWM"`).
#' @param trim.threshold `numeric(1)`. Fraction-of-peak cutoff,
#'   used when `abs.threshold` is `NULL`. Default `0.3`.
#' @param abs.threshold `numeric(1)` or `NULL`. Absolute-value
#'   cutoff on `sum_i |cwm[i,j]|`. When non-`NULL`, takes
#'   precedence over `trim.threshold`. Default `NULL`.
#' @param trim.from `character(1)`. One of `"both"`, `"left"`, or
#'   `"right"`. Direction of edge trimming.
#'
#' @return A trimmed CWM motif (single input) or list of trimmed
#'   CWM motifs (list input). The `type = "CWM"` tag is preserved.
#'
#' @examples
#' library(universalmotif)
#'
#' ## Synthetic CWM with a strong centre and weak flanks.
#' m <- matrix(c(0.02, -0.01,  0.01,  0.00,
#'               0.90, -0.30, -0.30, -0.30,
#'              -0.20,  0.80, -0.20, -0.20,
#'              -0.30, -0.30, -0.30,  0.90,
#'               0.03,  0.00, -0.02,  0.01),
#'             nrow = 4, byrow = FALSE,
#'             dimnames = list(c("A","C","G","T"), NULL))
#' cwm <- create_motif(m, type = "CWM", name = "demo")
#'
#' ## Default: drop edge columns where |colsum| < 0.3 * peak.
#' trimmed <- trim_cwm(cwm)
#' ncol(trimmed["motif"])
#'
#' ## Explicit absolute-value cutoff instead.
#' trim_cwm(cwm, abs.threshold = 0.5)
#'
#' ## One-sided trim.
#' trim_cwm(cwm, trim.from = "left")
#'
#' @seealso [trim_motifs()], [create_motif()], [convert_type()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
trim_cwm <- function(motifs, trim.threshold = 0.3, abs.threshold = NULL,
                     trim.from = c("both", "left", "right")) {

  ## --- arg validation ---------------------------------------------------
  if (!is.numeric(trim.threshold) || length(trim.threshold) != 1L ||
      is.na(trim.threshold) || trim.threshold < 0)
    stop("`trim.threshold` must be a single non-negative numeric", call. = FALSE)
  if (!is.null(abs.threshold) &&
      (!is.numeric(abs.threshold) || length(abs.threshold) != 1L ||
       is.na(abs.threshold) || abs.threshold < 0))
    stop("`abs.threshold` must be NULL or a single non-negative numeric",
         call. = FALSE)
  trim.from <- match.arg(trim.from)

  ## --- normalise input --------------------------------------------------
  was_list <- is.list(motifs)
  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, "character")
  else CLASS_IN <- .internal_convert(motifs)

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  in_types <- vapply(motifs, function(x) x@type, character(1))
  if (any(in_types != "CWM"))
    stop("`trim_cwm()` requires CWM motifs; got type `",
         paste(unique(in_types[in_types != "CWM"]), collapse = "`, `"),
         "`. Use `trim_motifs()` for non-CWM types.", call. = FALSE)

  ## --- trim each motif --------------------------------------------------
  trimmed <- lapply(motifs, function(m) trim_cwm_single(m,
                                                        trim.threshold,
                                                        abs.threshold,
                                                        trim.from))

  ## Drop fully-trimmed motifs with a message (matching trim_motifs()).
  dont_keep <- vapply(trimmed, is.null, logical(1))
  if (any(dont_keep)) {
    mot_names <- vapply(motifs[dont_keep], function(x) x@name, character(1))
    if (all(dont_keep))
      stop("All motifs were completely trimmed", call. = FALSE)
    message("The following motifs were completely trimmed: ",
            paste(mot_names, collapse = ", "))
  }
  trimmed <- trimmed[!dont_keep]

  if (length(trimmed) == 1 && !was_list) trimmed <- trimmed[[1]]
  trimmed <- .internal_convert(trimmed, unique(CLASS_IN))
  trimmed
}

## Single-motif trimmer. Returns NULL when no column survives.
trim_cwm_single <- function(motif, trim.threshold, abs.threshold, trim.from) {

  mat <- motif@motif
  if (ncol(mat) == 0L) return(NULL)

  col_abs <- colSums(abs(mat))

  ## Resolve the per-motif threshold.
  if (is.null(abs.threshold)) {
    peak <- max(col_abs)
    if (peak == 0) return(NULL)
    cutoff <- trim.threshold * peak
  } else {
    cutoff <- abs.threshold
  }

  keep <- col_abs >= cutoff

  ## Edge-only trim: keep all columns from the first surviving column
  ## to the last surviving column, depending on trim.from.
  surviving <- which(keep)
  if (length(surviving) == 0L) return(NULL)

  left  <- if (trim.from %in% c("both", "left"))  min(surviving) else 1L
  right <- if (trim.from %in% c("both", "right")) max(surviving) else ncol(mat)
  keep_range <- left:right
  if (length(keep_range) == 0L) return(NULL)

  motif@motif <- mat[, keep_range, drop = FALSE]

  ## Refresh derived slots so the trimmed motif validates cleanly.
  alph <- motif@alphabet
  pseudo <- motif@pseudocount
  bkg <- motif@bkg[rownames(motif@motif)]
  if (any(is.na(bkg))) bkg <- rep(1 / nrow(motif@motif), nrow(motif@motif))

  motif@icscore <- sum(apply(motif@motif, 2, position_icscoreC,
                             bkg = bkg, type = "CWM",
                             pseudocount = pseudo,
                             nsites = if (length(motif@nsites) > 0) motif@nsites else 100))

  if (alph %in% c("DNA", "RNA")) {
    consensus <- apply(motif@motif, 2, get_consensusC,
                       alphabet = alph, type = "CWM",
                       pseudocount = pseudo)
    colnames(motif@motif) <- consensus
    motif@consensus <- paste(consensus, collapse = "")
  } else if (alph == "AA") {
    consensus <- apply(motif@motif, 2, get_consensusAAC,
                       type = "CWM", pseudocount = pseudo)
    colnames(motif@motif) <- consensus
    motif@consensus <- paste(consensus, collapse = "")
  }

  validObject_universalmotif(motif)
  motif
}
