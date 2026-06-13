#' Plot motif logos using `compare_motifs2()` alignment (v2).
#'
#' `view_motifs2()` is the leaner counterpart of [view_motifs()]: it aligns
#' the input motifs using the same Pearson-correlation backend used by
#' [compare_motifs2()] and [merge_motifs2()].
#' The rendering (logo polygons, `facet_wrap()` layout, `sort.by` permutation,
#' `names.pos` placement) is identical to [view_motifs()].
#'
#' DNA and RNA only. For amino-acid or custom-alphabet motifs use
#' [view_motifs()].
#'
#' @param motifs See [convert_motifs()] for accepted motif formats.
#'   DNA or RNA only.
#' @param use.type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#'   Controls the matrix type used for display.
#' @param tryRC `logical(1)` Also test reverse-complement alignments.
#'   Default `TRUE`.
#' @param min.overlap `integer(1)` Minimum overlap (in columns) for an
#'   alignment to be accepted. Default `6`.
#' @param relative_entropy `logical(1)` If `TRUE`, compute information
#'   content using the motif's background frequencies (relative-entropy
#'   ICM). Used only for display when `use.type = "ICM"`.
#' @param nthreads `numeric(1)` Number of threads for comparison.
#'   `nthreads = 0` uses all available threads.
#' @param return.raw `logical(1)` Instead of returning a plot, return
#'   the aligned named matrices used to generate the plot.
#' @param dedup.names `logical(1)` Allow duplicate motif names to be
#'   silently disambiguated for plotting.
#' @param show.positions `logical(1)` Show x-axis position tick labels.
#' @param show.positions.once `logical(1)` When plotting multiple
#'   motifs, show x-axis position tick labels only once.
#' @param sort.by `character(1)`. Display order for multi-motif plots.
#'   One of `"none"` (default), `"ic"`, or `"similarity"`. Alignment
#'   is always performed with the highest-IC motif as the anchor;
#'   `sort.by` controls only the order in which the aligned motifs
#'   are shown.
#' @param show.names `logical(1)` Add motif names when plotting
#'   multiple motifs.
#' @param names.pos `character(1)` Motif name locations, either
#'   `"top"` or `"right"`.
#' @param use.freq `numeric(1)` Plot higher-order motifs from the
#'   `multifreq` slot.
#' @param colour.scheme `character` Named character vector of colour
#'   names. Default colours are provided for DNA, RNA, and AA motifs
#'   if `NULL`.
#' @param fontDF `data.frame` or `DataFrame` Polygon data for letters.
#' @param min.height `numeric(1)` Minimum letter height to show, as a
#'   fraction of total plot height.
#' @param x.spacer `numeric(1)` Horizontal spacing between letters.
#' @param y.spacer `numeric(1)` Vertical spacing between letters.
#' @param sort.positions `logical(1)` Sort letters vertically per
#'   position by height.
#' @param sort.positions.decreasing `logical(1)` Sort direction.
#' @param text.size `numeric(1)` Text size for labels.
#' @param fit.to.height `numeric(1)` Normalise per-position height.
#' @param RC.text `character(1)` Text appended to motif names shown
#'   as their reverse complement.
#' @param flip.neg `logical(1)` Flip letters with negative heights.
#'
#' @return A `ggplot` object. If `return.raw = TRUE`, a list of
#'   matrices.
#'
#' @details
#' The alignment step uses the same
#' backend as [merge_motifs2()] and [merge_similar2()]. The
#' highest-information-content input motif is chosen as the anchor,
#' every other motif is aligned against it once, motifs whose best
#' alignment lands on the reverse complement strand are
#' reverse-complemented before display, and all motifs are placed in
#' a shared column frame with zero-padding on either side as needed.
#' The padded matrices then flow through the same logo-rendering
#' pipeline used by [view_motifs()].
#'
#' See [view_motifs()] for the original version, which supports
#' Euclidean / similarity metrics and arbitrary alphabets.
#'
#' @examples
#' m1 <- create_motif("TTGACATA", name = "a")
#' m2 <- create_motif("CTTGACAT", name = "b")
#' m3 <- create_motif("TGACATAT", name = "c")
#' view_motifs2(list(m1, m2, m3))
#' view_motifs2(list(m1, m2, m3), sort.by = "similarity",
#'              names.pos = "right")
#'
#' @seealso [view_motifs()], [view_logo()], [compare_motifs2()],
#'   [merge_motifs2()], [merge_similar2()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
view_motifs2 <- function(motifs, use.type = "ICM", tryRC = TRUE,
  min.overlap = 6, relative_entropy = FALSE, nthreads = 1,
  return.raw = FALSE, dedup.names = TRUE, show.positions = TRUE,
  show.positions.once = TRUE, sort.by = c("none", "ic", "similarity"),
  show.names = TRUE, names.pos = c("top", "right"),
  use.freq = 1, colour.scheme = NULL, fontDF = NULL, min.height = 0.01,
  x.spacer = if (use.freq == 1) 0.04 else 0.1,
  y.spacer = 0.01, sort.positions = !use.type %in% c("PCM", "PPM"),
  sort.positions.decreasing = TRUE, text.size = 16,
  fit.to.height = if (use.type == "PPM") 1 else NULL, RC.text = " [RC]",
  flip.neg = FALSE) {

  ## --- arg validation ---------------------------------------------------
  if (!use.type %in% c("PPM", "ICM", "PWM", "PCM", "CWM"))
    stop("`use.type` must be one of `PCM`, `PPM`, `PWM`, `ICM`, or `CWM`",
         call. = FALSE)
  if (!isTRUEorFALSE(tryRC))
    stop("`tryRC` must be a single logical", call. = FALSE)
  if (!is.numeric(min.overlap) || length(min.overlap) != 1L ||
      min.overlap < 1L)
    stop("`min.overlap` must be a positive integer", call. = FALSE)
  if (!isTRUEorFALSE(relative_entropy))
    stop("`relative_entropy` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(return.raw))
    stop("`return.raw` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(dedup.names))
    stop("`dedup.names` must be a single logical", call. = FALSE)
  if (!isTRUEorFALSE(flip.neg))
    stop("`flip.neg` must be a single logical", call. = FALSE)

  nthreads <- resolve_nthreads(nthreads)

  names.pos <- match.arg(names.pos)
  sort.by   <- match.arg(sort.by)

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  ## When the user asks for CWM display, capture the original CWM
  ## matrices *before* convert_type_internal coerces them to PPM, so
  ## we can splice them back at display time. CWM display requires
  ## every input motif to actually be a CWM.
  cwm.mats <- NULL
  if (use.type == "CWM") {
    in_types <- vapply(motifs, function(x) x@type, character(1))
    if (any(in_types != "CWM"))
      stop("`use.type = \"CWM\"` requires every input motif to have ",
           "type = \"CWM\"", call. = FALSE)
    cwm.mats <- lapply(motifs, function(x) x@motif)
  }

  motifs <- convert_type_internal(motifs, "PPM")

  ## DNA / RNA only (matches merge_motifs2 / compare_motifs2_align_cpp).
  alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  mot.alph <- unique(alphs)
  if (length(mot.alph) > 1L)
    stop("all motifs must share the same alphabet", call. = FALSE)
  if (!mot.alph %in% c("DNA", "RNA"))
    stop("`view_motifs2()` supports DNA/RNA motifs only; got `",
         mot.alph, "`. Use `view_motifs()` for other alphabets.",
         call. = FALSE)

  ## Sort by IC for alignment: highest-IC motif becomes the anchor that
  ## every other motif is aligned to. `sort.by` controls only the
  ## *display* order, applied after alignment further below.
  if (length(motifs) > 1L) {
    ic.scores <- vapply(motifs, function(x) x@icscore, numeric(1))
    ic.order  <- order(ic.scores, decreasing = TRUE)
    motifs    <- motifs[ic.order]
    if (!is.null(cwm.mats)) cwm.mats <- cwm.mats[ic.order]
  } else {
    ic.order  <- 1L
  }

  if (use.type == "PWM" && !is.null(fit.to.height))
    stop("`fit.to.height` must be NULL if `use.type = \"PWM\"`")

  mot.nsites <- lapply(motifs, function(x) x@nsites)
  mot.pseudo <- lapply(motifs, function(x) x@pseudocount)

  ## Keep PPM matrices for the v2 aligner; build display matrices for
  ## the requested use.type.
  if (use.freq == 1) {
    mot.mats.ppm <- lapply(motifs, function(x) x@motif)
    mot.bkgs <- get_bkgs(motifs)
    mot.mats <- vector("list", length(mot.mats.ppm))
    for (i in seq_along(mot.mats)) {
      if (use.type == "CWM") {
        ## Splice the original CWM matrix back in instead of deriving
        ## a display matrix from the PPM-converted version.
        mot.mats[[i]] <- cwm.mats[[i]]
      } else {
        mot.mats[[i]] <- convert_mat_type_from_ppm(mot.mats.ppm[[i]], use.type,
          mot.nsites[[i]], mot.bkgs[[i]], mot.pseudo[[i]], relative_entropy)
      }
    }
  } else if (use.freq > 1) {
    if (any(vapply(motifs, function(x) length(x@multifreq) == 0, logical(1))))
      stop("missing multifreq slots")
    check_multi <- vapply(
      motifs,
      function(x) any(names(x@multifreq) %in% as.character(use.freq)),
      logical(1)
    )
    if (!any(check_multi))
      stop("not all motifs have corresponding multifreq matrix")
    mot.mats <- lapply(motifs,
      function(x) x@multifreq[[as.character(use.freq)]])
    mot.mats.ppm <- mot.mats
    if (use.type != "PPM") {
      for (i in seq_along(mot.mats)) {
        rn <- rownames(mot.mats[[i]])
        if (use.type == "PWM") {
          mot.mats[[i]] <- MATRIX_ppm_to_pwm(
            mot.mats[[i]], nsites = motifs[[i]]@nsites,
            pseudocount = motifs[[i]]@pseudocount,
            bkg = motifs[[i]]@bkg[rownames(mot.mats[[i]])])
        } else if (use.type == "ICM") {
          mot.mats[[i]] <- MATRIX_ppm_to_icm(
            mot.mats[[i]], bkg = motifs[[i]]@bkg[rownames(mot.mats[[i]])],
            relative_entropy = relative_entropy
          )
        } else if (use.type == "PCM") {
          mot.mats[[i]] <- MATRIX_ppm_to_pcm(
            mot.mats[[i]], nsites = motifs[[i]]@nsites
          )
        }
        rownames(mot.mats[[i]]) <- rn
      }
    }
    mot.bkgs <- get_bkgs(motifs)
  } else {
    stop("`use.freq` must be a positive integer")
  }

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  if (length(mot.names) != length(unique(mot.names))) {
    if (!dedup.names) stop(wmsg(
      "All motifs must have unique names when `dedup.names=FALSE`."
    ), call. = FALSE)
    mot.names <- make.unique(mot.names)
  }

  alph <- switch(mot.alph, "DNA" = DNA_BASES, "RNA" = RNA_BASES,
    rownames(mot.mats[[1]]))

  ylim2 <- NULL
  breaks <- NULL
  switch(use.type,
    "PPM" = {
      yname <- "probability"
      ylim2 <- 0:1
      breaks <- c(0, 0.5, 1)
    },
    "ICM" = {
      yname <- "bits"
      ylim2 <- c(0, log2(nrow(mot.mats[[1]])))
      breaks <- unique(round(c(0, ylim2[2] * 0.5, ylim2[2])))
    },
    "PWM" = {
      yname <- "logodds"
    },
    "PCM" = {
      yname <- "counts"
    },
    "CWM" = {
      yname <- "contribution"
    },
    stop("'use.type' must be one of 'PCM', 'PPM', 'PWM', 'ICM', 'CWM'")
  )

  if (is.null(fontDF)) fontDF <- fontDFroboto
  fontDF <- as.data.frame(fontDF)
  if (!any(c("x", "y", "order", "group") %in% colnames(fontDF)))
    stop(wmsg("Expected columns [x, y, order, group] in `fontDF` object"))
  fontDF <- fontDF[, c("x", "y", "order", "group")]

  if (is.null(colour.scheme)) {
    colour.scheme <- switch(mot.alph, "DNA" = DNA_COLOURS, "RNA" = RNA_COLOURS,
      NULL)
  } else {
    if (any(!names(colour.scheme) %in% alph))
      stop("colour.scheme must be a named vector with all possible letters present")
  }

  if (use.type %in% c("PWM", "ICM") && !is.null(fit.to.height)) {
    warning("`fit.to.height` is ignored for `use.type = c(\"ICM\", \"PWM\")`")
    fit.to.height <- NULL
  }

  if (length(motifs) == 1) {

    if (return.raw) {
      colnames(mot.mats[[1]]) <- NULL
      names(mot.mats) <- mot.names
      return(mot.mats)
    }

    if (use.type == "PWM" && any(is.infinite(mot.mats[[1]])))
      stop("Found non-finite values in motif")

    plotobj <- prep_single_motif_plot_data(mot.mats[[1]], use.type, fontDF,
      min.height, x.spacer, y.spacer, sort.positions, sort.positions.decreasing,
      fit.to.height, 1, flip.neg = flip.neg)

    breaks2 <- seq_len(ncol(mot.mats[[1]]))
    limits2 <- c(min(breaks2) - 0.501, max(breaks2) + 0.501)

    if (use.type == "PCM") {
      breaks <- c(0, round(max(colSums(mot.mats[[1]])) / 2), max(colSums(mot.mats[[1]])))
      ylim2 <- breaks[-2]
      breaks <- unique(breaks)
    } else if (use.type %in% c("PWM", "CWM")) {
      mot.mats.tmp <- mot.mats[[1]]
      mot.mats.tmp[mot.mats.tmp > 0] <- 0
      minval <- min(colSums(mot.mats.tmp))
      mot.mats.tmp <- mot.mats[[1]]
      mot.mats.tmp[mot.mats.tmp < 0] <- 0
      maxval <- max(colSums(mot.mats.tmp))
      breaks <- c(minval, minval / 2, 0, maxval / 2, maxval)
      if (use.type == "PWM") breaks <- round(breaks)
      breaks[1] <- breaks[1] - max(c(0.1 * abs(breaks[1]), 0.1))
      ylim2 <- c(breaks[1], maxval)
      breaks <- unique(breaks)
    }

    plotobj$y <- plotobj$y * 0.999

    plotobj <- ggplot(plotobj, aes(.data$x, .data$y, group = .data$letter.id, fill = .data$group)) +
      geom_polygon() +
      ylab(yname) +
      xlab(NULL) +
      theme_minimal() +
      scale_y_continuous(breaks = breaks, limits = ylim2, expand = c(0, 0)) +
      scale_x_continuous(breaks = breaks2, limits = limits2, expand = c(0.02, 0)) +
      theme(axis.line.y = element_line(linewidth = 0.25),
        panel.grid = element_blank(),
        text = element_text(size = text.size),
        axis.ticks.y = element_line(linewidth = 0.25), legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black", margin = margin(r = 1)))

    if (!show.positions) plotobj <- plotobj + theme(axis.text.x = element_blank())

    if (!is.null(colour.scheme)) plotobj <- plotobj +
      scale_fill_manual(values = colour.scheme[alph], limits = alph)

    if (use.type %in% c("PWM", "CWM"))
      plotobj <- plotobj + geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey75")

  } else {

    ## ---- v2 alignment: PCC via compare_motifs2_align_cpp -------------
    res <- view_motifs2_align(mot.mats.ppm, mot.mats, min.overlap, tryRC,
                              nthreads)
    which.rc <- res$motIsRC
    mots <- res$motifs

    for (i in seq_along(which.rc)) {
      if (which.rc[i]) mot.names[i + 1] <- paste0(mot.names[i + 1], RC.text)
    }
    names(mots) <- mot.names

    ## Apply display permutation. Alignment is done; `sort.by` decides
    ## the order in which the aligned motifs are shown.
    disp.perm <- switch(sort.by,
      "ic"         = seq_along(mots),
      "none"       = order(ic.order),
      "similarity" = {
        n <- length(mots)
        sim_mat <- matrix(1, n, n)
        for (i in seq_len(n - 1L)) {
          vi <- as.vector(mots[[i]])
          for (j in seq.int(i + 1L, n)) {
            r <- suppressWarnings(cor(vi, as.vector(mots[[j]])))
            if (!is.finite(r)) r <- 0
            sim_mat[i, j] <- r
            sim_mat[j, i] <- r
          }
        }
        hclust(as.dist(1 - sim_mat))$order
      }
    )
    if (!identical(as.integer(disp.perm), seq_along(mots))) {
      mots        <- mots[disp.perm]
      mot.names   <- mot.names[disp.perm]
      mot.mats    <- mot.mats[disp.perm]
      res$offsets <- res$offsets[disp.perm]
    }

    if (return.raw) return(mots)

    if (use.type == "PWM") {
      for (i in seq_along(mots)) {
        if (any(is.infinite(mots[[i]])))
          stop("Found non-finite values in motif")
      }
    }

    plotobj_multi <- vector("list", length(mots))
    for (i in seq_along(plotobj_multi)) {
      plotobj_multi[[i]] <- prep_single_motif_plot_data(
        mots[[i]], use.type, fontDF, min.height, x.spacer, y.spacer, sort.positions,
        sort.positions.decreasing, fit.to.height, 1, flip.neg = flip.neg
      )
      plotobj_multi[[i]]$motif.id <- names(mots)[i]
    }
    plotobj <- do.call(rbind, plotobj_multi)
    plotobj$y <- plotobj$y * 0.999

    breaks2 <- seq_len(ncol(mots[[1]]))
    limits2 <- c(min(breaks2) - 0.501, max(breaks2) + 0.501)

    mot_colsums <- vapply(mots, function(x) max(colSums(x)), numeric(1))
    if (use.type == "PCM") {
      breaks <- c(0, round(max(mot_colsums) / 2), max(mot_colsums))
      ylim2 <- breaks[-2]
      breaks <- unique(breaks)
    } else if (use.type %in% c("PWM", "CWM")) {
      mots.tmp <- mots
      minvals <- numeric(length(mots))
      for (i in seq_along(mots.tmp)) {
        mots.tmp[[i]][mots.tmp[[i]] > 0] <- 0
        minvals[i] <- min(colSums(mots.tmp[[i]]))
      }
      mots.tmp <- mots
      maxvals <- numeric(length(mots.tmp))
      for (i in seq_along(mots.tmp)) {
        mots.tmp[[i]][mots.tmp[[i]] < 0] <- 0
        maxvals[i] <- max(colSums(mots.tmp[[i]]))
      }
      breaks <- c(min(minvals), min(minvals) / 2, 0, max(maxvals) / 2, max(maxvals))
      if (use.type == "PWM") breaks <- round(breaks)
      breaks[1] <- breaks[1] - max(c(0.1 * abs(breaks[1]), 0.1))
      ylim2 <- c(breaks[1], max(maxvals))
      breaks <- unique(breaks)
    }

    plotobj$motif.id <- factor(plotobj$motif.id, levels = mot.names)

    if (!show.positions.once) {
      plabs <- res$offsets
      plotcounter <- new.env()
      plotcounter$p <- 0
      labfun <- function(x) {
        pcount <- get("p", envir = plotcounter)
        pcount <- pcount + 1
        assign("p", pcount, envir = plotcounter)
        y <- plabs[[pcount]]
        z <- seq(from = x[1], to = x[length(x)], by = 1)
        y <- c(rep(TRUE, y), rep(FALSE, ncol(mot.mats[[pcount]])))
        y <- c(y, rep(TRUE, length(x) - length(y)))
        z[y] <- ""
        z[!y] <- seq_len(sum(!y))
        z
      }
    }
    plotobj <- ggplot(plotobj, aes(.data$x, .data$y, group = .data$letter.id,
        fill = .data$group)) +
      geom_polygon() +
      ylab(yname) +
      xlab(NULL) +
      theme_minimal() +
      scale_y_continuous(breaks = breaks, limits = ylim2, expand = c(0, 0)) +
      scale_x_continuous(breaks = breaks2, limits = limits2, expand = c(0.02, 0),
        labels = if (!show.positions.once) labfun else waiver()) +
      theme(axis.line.y = element_line(linewidth = 0.25),
        panel.grid = element_blank(),
        text = element_text(size = text.size),
        axis.ticks.y = element_line(linewidth = 0.25), legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black", margin = margin(r = 1))) +
      facet_wrap(~motif.id, ncol = 1, strip.position = names.pos,
        scales = if (!show.positions.once) "free_x" else "fixed")

    if (!show.names && (!show.positions || show.positions.once))
      plotobj <- plotobj + theme(panel.spacing = unit(1, "lines"))
    else if ((!show.names || names.pos == "right") &&
             show.positions && show.positions.once)
      plotobj <- plotobj + theme(panel.spacing = unit(1, "lines"))

    if (show.names && names.pos == "right")
      plotobj <- plotobj +
        theme(strip.text.y.right = element_text(angle = 0, hjust = 0))

    if (!show.positions) plotobj <- plotobj + theme(axis.text.x = element_blank())

    if (!is.null(colour.scheme)) plotobj <- plotobj +
      scale_fill_manual(values = colour.scheme[alph], limits = alph)

    if (use.type %in% c("PWM", "CWM"))
      plotobj <- plotobj + geom_hline(yintercept = 0, size = 0.25, colour = "grey75")

    if (!show.names) plotobj <- plotobj + theme(strip.text = element_blank())

  }

  plotobj

}

#-----------------------------------------------------------------------------
# Private helper: align a list of motifs to a common column frame using
# compare_motifs2_align_cpp. Returns a list with the same shape as
# view_motifs_prep() so the rendering pipeline can be reused unchanged.
#
# - mot.mats.ppm: PPM matrices for the v2 aligner.
# - mot.mats.disp: display matrices (whatever use.type the caller asked
#                  for); the function pads / RC's these without further
#                  type conversion.
# - The first motif is taken as the anchor (caller IC-sorts).
#-----------------------------------------------------------------------------

view_motifs2_align <- function(mot.mats.ppm, mot.mats.disp, min.overlap,
                               tryRC, nthreads) {

  n <- length(mot.mats.disp)
  if (n == 1L) {
    return(list(motifs = mot.mats.disp, motIsRC = logical(0), offsets = 0L))
  }

  anchor <- 1L
  qi <- rep.int(anchor, n)
  ti <- seq_len(n)
  al <- compare_motifs2_align_cpp(mot.mats.ppm,
                                  qi          = as.integer(qi),
                                  ti          = as.integer(ti),
                                  min_overlap = as.integer(min.overlap),
                                  RC          = as.logical(tryRC),
                                  nthreads    = as.integer(nthreads))

  widths <- vapply(mot.mats.disp, ncol, integer(1))
  w.anchor <- widths[anchor]

  oriented_mats <- vector("list", n)
  frame_offset  <- integer(n)
  which.rc      <- logical(n)

  oriented_mats[[anchor]] <- mot.mats.disp[[anchor]]
  frame_offset[anchor] <- 0L

  for (i in seq_len(n)) {
    if (i == anchor) next
    L_i <- al$overlap[i]
    if (is.na(L_i) || L_i < min.overlap) {
      ## No alignment at min.overlap: drop the motif at offset 0
      ## (anchor-aligned by left edge). Matches view_motifs() v1
      ## behaviour for motifs that don't meet the threshold.
      oriented_mats[[i]] <- mot.mats.disp[[i]]
      frame_offset[i] <- 0L
      next
    }
    qstart_oriented <- al$q_start_oriented[i]
    t_start         <- max(0L, -al$offset[i])
    strand_neg      <- isTRUE(al$strand[i] == 1L)
    w_m             <- widths[i]

    if (strand_neg) {
      oriented_mats[[i]] <- rev_comp_mat(mot.mats.disp[[i]])
      qstart_orig        <- w.anchor - qstart_oriented - L_i
      t_start_oriented   <- w_m - t_start - L_i
      frame_offset[i]    <- qstart_orig - t_start_oriented
      which.rc[i]        <- TRUE
    } else {
      oriented_mats[[i]] <- mot.mats.disp[[i]]
      frame_offset[i]    <- qstart_oriented - t_start
    }
  }

  ## Common frame: smallest start to largest end across all motifs.
  frame_starts <- frame_offset
  frame_ends   <- frame_offset + vapply(oriented_mats, ncol, integer(1))
  frame_lo <- min(frame_starts)
  frame_hi <- max(frame_ends)

  padded  <- vector("list", n)
  offsets <- integer(n)
  for (i in seq_len(n)) {
    left_pad  <- frame_offset[i] - frame_lo
    right_pad <- frame_hi - frame_ends[i]
    motm <- oriented_mats[[i]]
    if (left_pad > 0L) {
      motm <- cbind(matrix(0, nrow = nrow(motm), ncol = left_pad,
                           dimnames = list(rownames(motm), NULL)),
                    motm)
    }
    if (right_pad > 0L) {
      motm <- cbind(motm,
                    matrix(0, nrow = nrow(motm), ncol = right_pad,
                           dimnames = list(rownames(motm), NULL)))
    }
    padded[[i]]  <- motm
    offsets[i]   <- left_pad
  }

  ## motIsRC convention from view_motifs_prep is length N-1 (anchor
  ## excluded). Match it.
  list(motifs  = padded,
       motIsRC = which.rc[-anchor],
       offsets = offsets)
}
