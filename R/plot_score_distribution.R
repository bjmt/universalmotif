#' Plot the distribution of motif scores.
#'
#' Draw the distribution of a motif's scores in random sequences, with the
#' score cutoffs that P-values give, to help choose a cutoff for
#' [scan_sequences()].
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param pvalues `numeric` One or several P-value cutoffs to mark, or `NULL`
#'   for none. A motif too short or too degenerate to reach a cutoff gets a
#'   warning instead of a line.
#' @param bkg.probs `numeric`, `list` Background probabilities, in any of the
#'   forms [motif_pvalue()] takes: a single vector for all motifs (named by
#'   letter or in alphabet order) or a list with one vector per motif. If
#'   missing, each motif's own `bkg` slot is used.
#' @param log.y `logical(1)` Draw probabilities on a log scale, which shows the
#'   upper tail where cutoffs usually fall.
#' @param hits Scores of real sites to mark as ticks along the bottom of each
#'   panel: the output of [scan_sequences()] or [scan_sequences_lite()] (a
#'   `DataFrame`, `data.frame` or `GRanges` with `motif` and `score` columns),
#'   a named list of scores with one entry per motif name, or for a single
#'   motif a `numeric` vector of scores.
#' @param bins `numeric(1)` Number of bins the distribution is gathered into.
#' @param ncol `numeric(1)` Number of panels per row when drawing several
#'   motifs.
#'
#' @return A `ggplot` object, with one panel per motif.
#'
#' @details
#' The distribution is the exact one that [motif_pvalue()] works from with
#' `method = "dynamic"`, gathered into `bins` bins of equal width between the
#' lowest and highest possible scores, so that the height of each bin is the
#' probability that a random site scores within it. Since scores lie on a grid
#' of steps of 0.001, the bin edges are placed on that grid and the probability
#' of each bin is the difference between the P-values of its two edges; the
#' probabilities thus add up to one. Scores are the same logodds scores
#' reported by [motif_pvalue()] and [scan_sequences()].
#'
#' Each cutoff in `pvalues` is drawn as a vertical line at the lowest score
#' with at most that P-value, labelled with the P-value, the score and the
#' score as a percentage of the maximum possible score (as in the `score.pct`
#' column of [scan_sequences()]). When a single motif is drawn, the top axis
#' also shows scores as a percentage of the maximum; since this percentage
#' differs between motifs, the top axis is left out when several motifs are
#' drawn as panels.
#'
#' @examples
#' data(ArabidopsisMotif)
#' data(ArabidopsisPromoters)
#' plot_score_distribution(ArabidopsisMotif, pvalues = c(1e-3, 1e-4))
#'
#' ## Mark the scores of hits along the bottom, on a log scale:
#' hits <- scan_sequences(ArabidopsisMotif, ArabidopsisPromoters,
#'   threshold = 1e-4, threshold.type = "pvalue", nthreads = 1)
#' plot_score_distribution(ArabidopsisMotif, hits = hits, log.y = TRUE)
#'
#' ## Several motifs, scored against an AT-rich background:
#' motifs <- read_meme(system.file("extdata", "meme_full.txt",
#'   package = "universalmotif"))
#' plot_score_distribution(motifs[c(1, 3)], pvalues = c(1e-3, 1e-5),
#'   bkg.probs = c(0.32, 0.18, 0.18, 0.32), log.y = TRUE)
#'
#' @seealso [motif_pvalue()], [motif_range()], [scan_sequences()],
#'   [scan_sequences_lite()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
plot_score_distribution <- function(motifs, pvalues = 1e-4, bkg.probs,
  log.y = FALSE, hits = NULL, bins = 100, ncol = 2) {

  if (!is.null(pvalues) &&
      (!is.numeric(pvalues) || anyNA(pvalues) || any(pvalues <= 0 | pvalues > 1)))
    stop("`pvalues` must lie above 0 and at most 1", call. = FALSE)
  if (!is.numeric(bins) || length(bins) != 1 || is.na(bins) || bins < 1 ||
      bins != round(bins))
    stop("`bins` must be a positive whole number", call. = FALSE)
  if (!is.numeric(ncol) || length(ncol) != 1 || is.na(ncol) || ncol < 1)
    stop("`ncol` must be a positive number", call. = FALSE)
  if (!is.logical(log.y) || length(log.y) != 1 || is.na(log.y))
    stop("`log.y` must be TRUE or FALSE", call. = FALSE)

  if (is.list(motifs) && !length(motifs))
    stop("no motifs given", call. = FALSE)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  names <- make.unique(vapply(motifs, function(x) x@name, character(1)))

  if (missing(bkg.probs) || is.null(bkg.probs)) {
    bkgs <- rep(list(NULL), length(motifs))
  } else if (is.list(bkg.probs)) {
    if (length(bkg.probs) == 1) {
      bkgs <- rep(bkg.probs, length(motifs))
    } else if (length(bkg.probs) == length(motifs)) {
      bkgs <- bkg.probs
    } else {
      stop(wmsg("`bkg.probs` must be a single background for all motifs or ",
          "a list with one per motif"), call. = FALSE)
    }
  } else if (is.numeric(bkg.probs)) {
    bkgs <- rep(list(bkg.probs), length(motifs))
  } else {
    stop("`bkg.probs` must be a numeric vector or a list", call. = FALSE)
  }

  hit.scores <- score_dist_hits(hits, names)
  pvalues <- as.numeric(pvalues)

  panels <- lapply(seq_along(motifs), function(i) {
    score_dist_panel(motifs[[i]], names[i], bkgs[[i]], pvalues, bins,
      hit.scores[[i]])
  })

  bars <- do.call(rbind, lapply(panels, `[[`, "bars"))
  steps <- do.call(rbind, lapply(panels, `[[`, "steps"))
  cutoffs <- do.call(rbind, lapply(panels, `[[`, "cutoffs"))
  rug <- do.call(rbind, lapply(panels, `[[`, "hits"))
  keys <- do.call(rbind, lapply(panels, `[[`, "key"))
  # Invisible points spanning each panel's limits, since the x-axis must end
  # exactly at the distribution's edges and the y-axis leaves room for labels.
  limits <- do.call(rbind, lapply(panels, `[[`, "limits"))

  facet.levels <- function(x) {
    x$motif <- factor(x$motif, levels = names)
    x
  }
  bars <- facet.levels(bars)
  steps <- facet.levels(steps)
  limits <- facet.levels(limits)

  if (log.y) {
    limits$y <- limits$y.log
    floors <- vapply(panels, `[[`, numeric(1), "floor")
    bars$ymin <- floors[as.integer(bars$motif)]
    bars$ymax <- pmax(bars$ymax, bars$ymin)
    steps$y <- pmax(steps$y, floors[as.integer(steps$motif)])
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = limits,
      ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_rect(data = bars,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
        ymin = .data$ymin, ymax = .data$ymax),
      fill = "#C9CDD3", colour = NA) +
    ggplot2::geom_step(data = steps,
      ggplot2::aes(x = .data$x, y = .data$y),
      direction = "hv", colour = "#4A4F57", linewidth = 0.35)

  if (!is.null(rug)) {
    rug <- facet.levels(rug)
    keys <- facet.levels(keys)
    p <- p +
      ggplot2::geom_rug(data = rug,
        ggplot2::aes(x = .data$x, alpha = .data$alpha),
        sides = "b", colour = "#255C99", linewidth = 0.45,
        length = grid::unit(0.06, "npc")) +
      ggplot2::geom_text(data = keys,
        ggplot2::aes(x = -Inf, y = Inf, label = .data$label),
        hjust = -0.08, vjust = 1.6, colour = "#255C99",
        size = 8 / ggplot2::.pt) +
      ggplot2::scale_alpha_identity()
  }

  if (!is.null(cutoffs)) {
    cutoffs <- facet.levels(cutoffs)
    p <- p +
      ggplot2::geom_vline(data = cutoffs,
        ggplot2::aes(xintercept = .data$x),
        colour = "#C0504D", linewidth = 0.5) +
      ggplot2::geom_label(data = cutoffs,
        ggplot2::aes(x = .data$x, y = Inf, label = .data$label),
        angle = 90, hjust = 1.03, vjust = -0.3, size = 7.5 / ggplot2::.pt,
        colour = "#333333", fill = grDevices::adjustcolor("white", 0.85),
        label.r = grid::unit(0, "pt"), label.padding = grid::unit(1.2, "pt"),
        border.colour = NA)
  }

  if (length(motifs) == 1 && panels[[1]]$top > 0) {
    top <- panels[[1]]$top
    sec <- ggplot2::sec_axis(transform = ~ 100 * . / top,
      name = "% of maximum score")
  } else {
    sec <- ggplot2::waiver()
  }

  p <- p +
    ggplot2::scale_x_continuous(expand = c(0, 0), sec.axis = sec) +
    ggplot2::facet_wrap(~ motif, ncol = min(ncol, length(motifs)),
      scales = "free") +
    ggplot2::labs(x = "Score",
      y = if (log.y) "Probability (log10)" else "Probability") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background   = ggplot2::element_rect(fill = NA, colour = NA),
      strip.text         = ggplot2::element_text(size = 11),
      strip.placement    = "outside",
      plot.margin        = ggplot2::margin(5.5, 12, 5.5, 5.5),
      panel.grid         = ggplot2::element_blank(),
      panel.border       = ggplot2::element_blank(),
      axis.line.x.bottom = ggplot2::element_line(colour = "black"),
      axis.line.x.top    = ggplot2::element_line(colour = "black"),
      axis.line.y.left   = ggplot2::element_line(colour = "black")
    )

  if (log.y) {
    p <- p + ggplot2::scale_y_log10(expand = c(0, 0),
      breaks = score_dist_log_breaks, labels = score_dist_log_labels)
  } else {
    p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
  }

  p

}

score_dist_panel <- function(motif, name, bkg, pvalues, bins, hits) {

  scores <- score_dist_pvalue(motif, bkg, pvalue = c(1, pvalues))
  dist <- score_dist_bins(motif, bkg, bins, scores[1])
  edges <- dist$edges
  mass <- dist$mass
  top <- dist$top
  n <- length(mass)

  bars <- data.frame(motif = name, xmin = edges[-(n + 1)], xmax = edges[-1],
    ymin = 0, ymax = mass, stringsAsFactors = FALSE)
  steps <- data.frame(motif = name, x = edges, y = c(mass, mass[n]),
    stringsAsFactors = FALSE)

  positive <- mass[mass > 0]
  floor <- if (length(positive)) min(positive) / 3 else 1e-300
  limits <- data.frame(motif = name, x = edges[c(1, n + 1)],
    y = c(0, max(mass) * 1.32), y.log = c(floor, max(mass) * 4),
    stringsAsFactors = FALSE)

  cutoffs <- NULL
  thresholds <- scores[-1]
  reachable <- is.finite(thresholds)
  if (any(reachable)) {
    pct <- if (top > 0) sprintf(" (%.0f%%)", 100 * thresholds / top) else ""
    cutoffs <- data.frame(motif = name, x = thresholds, pvalue = pvalues,
      label = paste0("P = ", format_pvalue_label(pvalues), ", score ",
        sprintf("%.2f", thresholds), pct), stringsAsFactors = FALSE)
    cutoffs <- cutoffs[reachable, , drop = FALSE]
  }
  if (any(!reachable)) {
    warning(wmsg(name, " cannot reach a P-value of ",
        paste(format_pvalue_label(pvalues[!reachable]), collapse = ", "),
        ", so no cutoff is drawn for it; the highest score is ",
        sprintf("%.3f", top)), call. = FALSE)
  }

  hit.df <- key <- NULL
  if (length(hits)) {
    hit.df <- data.frame(motif = name, x = hits,
      alpha = max(0.15, min(1, 20 / length(hits))), stringsAsFactors = FALSE)
    key <- data.frame(motif = name,
      label = paste0("hits (n = ", length(hits), ")"), stringsAsFactors = FALSE)
  }

  list(bars = bars, steps = steps, cutoffs = cutoffs, hits = hit.df,
    key = key, limits = limits, floor = floor, top = top)

}

# Bin edges over the motif's scores, the exact probability of a random site
# scoring within each bin, and the highest score. Scores lie on a grid of
# steps of 0.001, so the edges are put on that grid with each bin holding the
# same number of grid points. The first edge has a P-value of 1 and the last
# one, past the highest score, of 0, so the probabilities add up to 1.
score_dist_bins <- function(motif, bkg, bins, lo) {
  lo.g <- round(lo * 1000)
  hi.g <- max(round(score_dist_max(motif, bkg) * 1000), lo.g)
  npts <- hi.g - lo.g + 1
  step <- ceiling(npts / bins)
  nbins <- ceiling(npts / step)
  grid <- lo.g + step * seq(0, nbins)
  grid[nbins + 1] <- hi.g + 1
  tail <- c(1, rep(0, nbins))
  if (nbins > 1) {
    tail[2:nbins] <- score_dist_pvalue(motif, bkg, score = grid[2:nbins] / 1000)
  }
  mass <- pmax(-diff(tail), 0)
  list(edges = grid / 1000, mass = mass, top = hi.g / 1000)
}

# motif_pvalue() against the motif's own background when bkg is NULL.
score_dist_pvalue <- function(motif, bkg, ...) {
  if (is.null(bkg)) motif_pvalue(motif, ...)
  else motif_pvalue(motif, ..., bkg.probs = bkg)
}

# The highest score a motif can give against a background. As in
# motif_pvalue(), a motif is scored against another background by rebasing
# its probabilities onto that background.
score_dist_max <- function(motif, bkg) {
  if (is.null(bkg)) return(unname(motif_range(motif)["max"]))
  bkg <- motif_pvalue_bkg(motif, bkg, use.freq = 1)
  motif <- convert_type_internal(motif, "PPM")
  motif["bkg"] <- bkg
  unname(suppressMessages(motif_range(motif))["max"])
}

# The hit scores for each motif, from a hit table, a named list of scores, or
# a vector of scores for a single motif.
score_dist_hits <- function(hits, names) {
  if (is.null(hits)) return(rep(list(NULL), length(names)))
  if (methods::is(hits, "Vector") && !methods::is(hits, "DataFrame"))
    hits <- S4Vectors::mcols(hits)
  if (is.data.frame(hits) || methods::is(hits, "DataFrame")) {
    cols <- colnames(hits)
    if (!"score" %in% cols)
      stop(wmsg("`hits` has no `score` column; was it made by ",
          "scan_sequences()?"), call. = FALSE)
    score <- as.numeric(hits$score)
    if (!"motif" %in% cols) {
      if (length(names) > 1)
        stop(wmsg("`hits` needs a `motif` column to be split between ",
            "several motifs"), call. = FALSE)
      return(list(score[!is.na(score)]))
    }
    motif <- as.character(hits$motif)
    return(lapply(names, function(x) {
      s <- score[motif == x]
      s[!is.na(s)]
    }))
  }
  if (is.list(hits)) {
    if (is.null(names(hits)) || any(names(hits) == ""))
      stop("a list of `hits` must be named by motif", call. = FALSE)
    unknown <- setdiff(names(hits), names)
    if (length(unknown))
      stop(wmsg("`hits` names motifs that are not plotted: ",
          paste(utils::head(sort(unknown), 5), collapse = ", ")),
        call. = FALSE)
    return(lapply(names, function(x) {
      if (x %in% names(hits)) as.numeric(hits[[x]]) else NULL
    }))
  }
  if (!is.numeric(hits))
    stop(wmsg("`hits` must be the output of scan_sequences(), a named list ",
        "of scores or a numeric vector"), call. = FALSE)
  if (length(names) > 1)
    stop(wmsg("give `hits` for several motifs as a table with a `motif` ",
        "column or a list of scores named by motif"), call. = FALSE)
  list(as.numeric(hits))
}

# Breaks at powers of ten on the log scale, thinned to at most about six,
# labelled as 10 to a power.
score_dist_log_breaks <- function(limits) {
  k <- seq(ceiling(log10(limits[1])), floor(log10(limits[2])))
  step <- ceiling(length(k) / 6)
  if (step > 1 && any(k %% step == 0)) k <- k[k %% step == 0]
  10^k
}

score_dist_log_labels <- function(breaks) {
  out <- paste0("10^", round(log10(breaks)))
  out[is.na(breaks)] <- "''"
  parse(text = out)
}

# P-values for labels: plain above 1e-3, otherwise as a short mantissa and
# exponent such as 2.5e-6.
format_pvalue_label <- function(p) {
  vapply(p, function(x) {
    if (x >= 1e-3) return(sprintf("%g", x))
    ex <- floor(log10(x))
    mant <- round(x / 10^ex, 2)
    if (mant >= 10) {
      mant <- mant / 10
      ex <- ex + 1
    }
    if (mant == 1) paste0("1e", ex) else paste0(sprintf("%g", mant), "e", ex)
  }, character(1))
}
