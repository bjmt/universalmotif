#' Plot motif logos.
#'
#' Show sequence logo. If given a list of more than one motif, then the motifs
#' are aligned with the first in the list.
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param use.type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#' @param return.raw `logical(1)` Instead of returning a plot, return the
#'    aligned named matrices used to generate the plot. This can be useful
#'    if you wish to use [view_motifs()] alignment capabilities for custom
#'    plotting uses. Alignment is performed by adding empty columns to the
#'    left or right of motifs to generate matrices of equal length.
#' @param dedup.names `logical(1)` Plotting motifs with duplicated names is
#'    not allowed. Setting this to `TRUE` allows the names to be modified
#'    for plotting.
#' @param ... Additional options for [ggseqlogo::geom_logo()].
#' @param show.positions `logical(1)` Show x-axis position tick labels.
#' @param show.positions.once `logical(1)` When plotting multiple motifs,
#'    show x-axis position tick labels only once.
#' @param show.names `logical(1)` Add motif names when plotting multiple
#'    motifs.
#'
#' @return A ggplot object. If `return.raw = TRUE`, a list.
#'
#' @details
#' Since the \pkg{ggseqlogo} package can only plot individual
#' characters and not strings, plotting the `multifreq` slot is not
#' supported. See the `examples` section for plotting the `multifreq`
#' slot using the \pkg{Logolas} package.
#'
#' See [compare_motifs()] for more info on comparison parameters.
#'
#' Note: `score.strat = "a.mean"` is NOT recommended, as [view_motifs()] will
#' not discriminate between two alignments with equal mean scores, even if one
#' alignment is longer than the other.
#'
#' @examples
#' ## plotting multifreq motifs:
#' \dontrun{
#'   motif <- create_motif()
#'   motif <- add_multifreq(motif, sample_sites(motif))
#'   Logolas::logomaker(motif["multifreq"][["2"]], type = "Logo",
#'                      color_type = "per_symbol")
#' }
#'
#' @references
#' \insertRef{logolas}{universalmotif}
#'
#' \insertRef{ggseqlogo}{universalmotif}
#'
#' @seealso [compare_motifs()], [add_multifreq()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams compare_motifs
#' @export
view_motifs <- function(motifs, use.type = "ICM", method = "ALLR",
  tryRC = TRUE, min.overlap = 6, min.mean.ic = 0.25,
  relative_entropy = FALSE, normalise.scores = FALSE, min.position.ic = 0,
  score.strat = "sum", return.raw = FALSE, dedup.names = FALSE,
  show.positions = TRUE, show.positions.once = TRUE, show.names = TRUE, ...) {

# view_motifs(ArabidopsisMotif) +
#   theme(axis.line.y = element_line(size = 0.25),
#     axis.ticks.y = element_line(size = 0.25),
#     axis.text.y = element_text(margin = margin(r = 1))) +
#   scale_y_continuous(breaks=c(0,1,2), limits = c(0, 2), expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0.02, 0), breaks = 1:15)
#
# view_motifs(c(ArabidopsisMotif, create_motif())) +
#   theme(axis.line.y = element_line(size = 0.25),
#     axis.ticks.y = element_line(size = 0.25), strip.text = element_blank(),
#     panel.spacing = unit(1, "lines"),
#     axis.text.y = element_text(margin = margin(r = 1))) +
#   scale_y_continuous(breaks=c(0,1,2), limits = c(0, 2), expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0.02, 0), breaks = 1:15) +  # OR breaks = NULL
#   facet_wrap(~seq_group, scales = "free_x", ncol = 1)
#
# - adding to motif_tree(): first find motif which is most similar to all input
#   motifs, and make sure all the motifs are aligned to that motif!
# - ylim: round(c(0, log2(nrow(motif)) / 2, log2(nrow(motif))))

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!method %in% COMPARE_METRICS) {
    method_check <- " * Incorrect 'method'"
    all_checks <- c(all_checks, method_check)
  }
  if (!use.type %in% c("PPM", "ICM", "PWM", "PCM")) {
    use.type_check <- paste0(" * Incorrect 'use.type': expected `PPM`, `PCM`, ",
                             "`PWM` or `ICM`; got `",
                             use.type, "`")
    use.type_check <- wmsg2(use.type_check, 4, 2)
    all_checks <- c(all_checks, use.type_check)
  }
  char_check <- check_fun_params(list(use.type = args$use.type,
                                      method = args$method,
                                      score.strat = args$score.strat),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     min.position.ic = args$min.position.ic),
                                numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores,
                                      return.raw = args$return.raw,
                                      dedup.names = args$dedup.names),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(all_checks, char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (!score.strat %in% c("sum", "a.mean", "g.mean", "median", "wa.mean",
                          "wg.mean", "fzt"))
    stop("'score.strat' must be one of 'sum', 'a.mean', 'g.mean', 'median', ",
         "'wa.mean', 'wg.mean', 'fzt'")

  if (score.strat %in% c("g.mean", "wg.mean") && method %in%
      c("ALLR", "ALLR_LL", "PCC"))
    stop(wmsg("'g.mean'/'wg.mean' is not allowed for methods which can generate negative values: ",
              "ALLR, ALLR_LL, PCC"))

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PPM")
  if (!is.list(motifs)) motifs <- list(motifs)

  ylim2 <- NULL
  breaks <- NULL
  if (use.type == "ICM" && !relative_entropy) {
    plot.method <- "bits"
    yname <- "bits"
    ylim2 <- log2(nrow(motifs[[1]]@motif))
    breaks <- unique(round(c(0, ylim2 * 0.5, ylim2)))
  } else {
    switch(use.type,
      "PPM" = {
        plot.method <- "prob"
        yname <- "probability"
        breaks <- c(0, 0.5 ,1)
        ylim2 <- 1
      },
      "ICM" = {
        plot.method <- "custom"
        yname <- "bits"
        ylim2 <- log2(nrow(motifs[[1]]@motif))
        breaks <- unique(round(c(0, ylim2 * 0.5, ylim2)))
      },
      "PWM" = {
        plot.method <- "custom"
        yname <- "logodds"
        if (length(motifs) > 1 && method == "KL") {
          stop("cannot use method 'KL' with 'PWM' matrices")
        }
      },
      "PCM" = {
        plot.method <- "custom"
        yname <- "counts"
      },
      stop("'use.type' must be one of 'PCM', 'PPM', 'PWM', 'ICM'")
    )
  }

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  if (length(mot.names) != length(unique(mot.names))) {
    if (!dedup.names) stop(wmsg(
      "All motifs must have unique names. Alternatively, set dedup.names=TRUE."
    ), call. = FALSE)
    tofix <- duplicated(mot.names)
    mot.names[tofix] <- paste0(
      mot.names[tofix], " (duplicated #", seq_len(sum(tofix)), ")"
    )
  }

  mot.mats <- lapply(motifs, function(x) x@motif)

  mot.alph <- unique(vapply(motifs, function(x) x@alphabet, character(1)))
  if (length(mot.alph) > 1) stop("can only have one alphabet")
  use.custom <- FALSE

  switch(mot.alph,
    "DNA" = {
      alph <- DNA_BASES
      seq_type <- "dna"
    },
    "RNA" = {
      alph <- RNA_BASES
      seq_type <- "rna"
    },
    "AA" = {
      alph <- AA_STANDARD2
      seq_type <- "aa"
    },
    {
      if (mot.alph != "custom") {
        alph <- sort_unique_cpp(safeExplode(mot.alph))
        use.custom <- TRUE
      } else {
        alph <- rownames(mot.mats[[1]])
        use.custom <- TRUE
      }
    }
  )

  mot.bkgs <- get_bkgs(motifs)
  mot.nsites <- lapply(motifs, function(x) x@nsites)
  mot.pseudo <- lapply(motifs, function(x) x@pseudocount)

  if (use.type == "ICM") {
    for (i in seq_along(mot.mats)) {
      if (all(mot.mats[[i]] == 0.25)) {
        stop(wmsg("Zero IC motifs cannot be plotted with ggseqlogo if `use.type = \"ICM\"`"))
      }
    }
  }

  if (length(motifs) == 1) {
    mot.mats[[1]] <- convert_mat_type_from_ppm(mot.mats[[1]], use.type, mot.nsites[[1]],
                                               mot.bkgs[[1]], mot.pseudo[[1]],
                                               relative_entropy)
    if (return.raw) {
      colnames(mot.mats[[1]]) <- NULL
      names(mot.mats) <- mot.names
      return(mot.mats)
    }
    if (use.custom) {
      p <- ggseqlogo(mot.mats[[1]], method = plot.method,
                     seq_type = seq_type, namespace = alph, ...) +
             ylab(yname)
    } else {
      p <- ggseqlogo(mot.mats[[1]], method = plot.method,
                     seq_type = seq_type, ...) +
             ylab(yname)
    }
    p <- suppressMessages(
      p +
        scale_y_continuous(breaks = if (!is.null(breaks)) breaks else NULL,
          limits = if (!is.null(ylim2)) c(0, ylim2) else NULL,
          expand = c(0, 0))
    )
    if (!show.positions) p <- p + xlab(element_blank())
    p <- suppressMessages(
      p +
        theme(axis.line.y = element_line(size = 0.25),
          axis.ticks.y = element_line(size = 0.25),
          axis.text.y = element_text(margin = margin(r = 1))) +
        scale_x_continuous(breaks = seq_len(ncol(mot.mats[[1]])),
          expand = c(0.02, 0))
    )
    return(p)
  }

  res <- view_motifs_prep(mot.mats, method, tryRC, min.overlap, min.mean.ic,
                          min.position.ic, mot.bkgs, relative_entropy,
                          normalise.scores, alph, get_nsites(motifs),
                          score.strat)
  which.rc <- res$motIsRC
  mots <- res$motifs
  mots <- check_mot_sizes(mots)

  if (method %in% c("KL", "ALLR", "ALLR_LL")) {
    for (i in seq_along(mots)) {
      mots[[i]][mots[[i]] > 0] <- mots[[i]][mots[[i]] > 0] - 0.01
    }
  }

  mots <- mapply(function(x1, x2, x3, x4)
                   convert_mat_type_from_ppm(x1, use.type, x2, x3, x4, relative_entropy),
                 mots, mot.nsites, mot.bkgs, mot.pseudo, SIMPLIFY = FALSE)

  for (i in seq_along(which.rc)) {
    if (which.rc[i]) mot.names[i + 1] <- paste(mot.names[i + 1], "[RC]")
  }
  names(mots) <- mot.names

  if (return.raw) return(mots)

  if (use.custom) {
    p <- ggplot() + geom_logo(mots, method = plot.method, seq_type = seq_type, 
                         namespace = alph, ...) +
      theme_logo() +
      facet_wrap(~seq_group, ncol = 1,
        scales = if (!show.positions.once) "free_x" else "fixed") +
      ylab(yname)
  } else {
    p <- ggplot() + geom_logo(mots, method = plot.method, seq_type = seq_type, ...) +
      theme_logo() +
      facet_wrap(~seq_group, ncol = 1,
        scales = if (!show.positions.once) "free_x" else "fixed") +
      ylab(yname)
  }

  p <- suppressMessages(
    p +
      scale_y_continuous(breaks = if (!is.null(breaks)) breaks else NULL,
        limits = if (!is.null(ylim2)) c(0, ylim2) else NULL,
        expand = c(0, 0))
  )
  p <- suppressMessages(
    p +
      theme(axis.line.y = element_line(size = 0.25),
        axis.ticks.y = element_line(size = 0.25),
        panel.spacing = unit(1, "lines"),
        axis.text.y = element_text(margin = margin(r = 1))) +
      scale_x_continuous(expand = c(0.02, 0),
        breaks = if (show.positions) seq_len(ncol(mots[[1]])) else NULL)
  )
  if (!show.names) p <- p + theme(strip.text = element_blank())

  p

}

check_mot_sizes <- function(mots) {

  sizes <- vapply(mots, ncol, integer(1))
  msize <- max(sizes)

  if (length(unique(sizes)) == 1) {
    mots <- check_right_side(mots, msize)
    return(mots)
  }

  for (i in seq_along(sizes)) {
    if (sizes[i] < msize) {
      mots[[i]] <- cbind(mots[[i]], matrix(0, nrow = nrow(mots[[i]]),
                                           ncol = msize - sizes[i]))
    }
  }

  check_right_side(mots, msize)

}

check_right_side <- function(mots, msize) {

  ok <- FALSE
  thiscol <- msize

  while (!ok) {

    for (i in seq_along(mots)) {
      if (any(mots[[i]][, thiscol] > 0)) ok <- TRUE
    }

    if (!ok) {
      for (i in seq_along(mots)) {
        mots[[i]] <- mots[[i]][, -thiscol]
      }
      thiscol <- thiscol - 1
    }

  }

  mots

}

convert_mat_type_from_ppm <- function(mot.mat, type, nsites, bkg, pseudocount,
                                      relative_entropy) {

  which.zero <- apply(mot.mat, 2, function(x) any(x != 0))

  if (length(nsites) == 0 || nsites == 1) nsites <- 100
  if (length(pseudocount) == 0) pseudocount <- 1

  mot.mat2 <- mot.mat[, which.zero, drop = FALSE]
  mot.mat2 <- switch(type,
                       "PCM" = apply(mot.mat2,  2,
                                     ppm_to_pcmC, nsites = nsites),
                       "PWM" = apply(mot.mat2, 2,
                                     ppm_to_pwmC, bkg = bkg,
                                     pseudocount = pseudocount,
                                     nsites = nsites),
                       "ICM" = apply(mot.mat2, 2,
                                     ppm_to_icmC, bkg = bkg,
                                     relative_entropy = relative_entropy),
                        mot.mat2
                     )

  mot.mat[, which.zero] <- mot.mat2

  mot.mat

}
