#' Plot motif logos.
#'
#' Show sequence logo. If given a list of more than one motif, then the motifs
#' are aligned with the first in the list.
#'
#' @param motifs See [convert_motifs()] for acceptable motif formats.
#' @param use.type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#' @param method `character(1)` One of `c('PCC', 'MPCC', 'EUCL', 'MEUCL',
#'    'SW', 'MSW', 'KL', 'MKL')`. See [compare_motifs()].
#' @param tryRC `logical(1)` Check if motif reverse complement leads to a
#'    better alignment. See [compare_motifs()].
#' @param min.overlap `numeric(1)` Minimum alignment overlap between
#'    motifs. If `min.overlap < 1`, this represents the minimum fraction
#'    between the two motifs during alignment. See [compare_motifs()].
#' @param min.mean.ic `numeric(1)` Minimum information content between the
#'    two motifs for an alignment to be scored. This helps prevent scoring
#'    alignments between low information content regions of two motifs. See
#'    [compare_motifs()].
#' @param relative_entropy `logical(1)` For ICM calculation. See
#'    [convert_type()].
#' @param normalise.scores `logical(1)` Favour alignments which leave fewer
#'    unaligned positions. See [compare_motifs()].
#' @param min.position.ic `numeric(1)` Minimum information content required between
#'    individual alignment positions for it to be counted in the final alignment
#'    score. It is recommended to use this together with `normalise.scores = TRUE`,
#'    as this will help punish scores resulting from only a fraction of an
#'    alignment.
#' @param ... Additional options for [ggseqlogo::geom_logo()].
#'
#' @return A ggplot object.
#'
#' @details
#' Since the \pkg{ggseqlogo} package can only plot individual
#' characters and not strings, plotting the `multifreq` slot is not
#' supported. See the `examples` section for plotting the `multifreq`
#' slot using the \pkg{Logolas} package.
#'
#' @examples
#' ## plotting multifreq motifs:
#' \dontrun{
#'   motif <- create_motif()
#'   motif <- add_multifreq(motif, sample_sites(motif))
#'   Logolas::logomaker(motif["multifreq"][[as.character(2)]], type = "Logo",
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
#' @export
view_motifs <- function(motifs, use.type = "ICM", method = "MPCC",
                        tryRC = TRUE, min.overlap = 6, min.mean.ic = 0.25,
                        relative_entropy = FALSE, normalise.scores = FALSE,
                        min.position.ic = 0, ...) {

  # Idea for multifreq plotting: just manually assign heights for all letters
  # in a multifreq position. Perhaps also manually add spacing between
  # individual multifreq positions (which are really just multiple positions
  # in sync).

  # Possible bug: min.overlap not being respected?

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!method %in% c("PCC", "MPCC", "EUCL", "MEUCL", "SW", "MSW", "KL",
                     "MKL")) {
    method_check <- paste0(" * Incorrect 'method': expected `PCC`, `MPCC`, ",
                           "`EUCL`, `MEUCL`, `SW`, `MSW`, `KL` or `MKL`; got `",
                           method, "`")
    method_check <- wmsg2(method_check, 4, 2)
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
                                      method = args$method),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic,
                                     min.position.ic = args$min.position.ic),
                                numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(all_checks, char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PPM")
  if (!is.list(motifs)) motifs <- list(motifs)

  if (use.type == "ICM" && !relative_entropy) {
    plot.method <- "bits"
    yname <- "bits"
  } else {
    switch(use.type,
      "PPM" = {
        plot.method <- "prob"
        yname <- "probability"
      },
      "ICM" = {
        plot.method <- "custom"
        yname <- "bits"
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
  if (length(mot.names) != length(unique(mot.names)))
    stop("All motifs must have unique names")

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
      alph <- AA_STANDARD
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

  mot.bkgs <- lapply(motifs, function(x) x@bkg[seq_along(alph)])
  mot.nsites <- lapply(motifs, function(x) x@nsites)
  mot.pseudo <- lapply(motifs, function(x) x@pseudocount)

  if (length(motifs) == 1) {
    mot.mats[[1]] <- convert_mat_type_from_ppm(mot.mats[[1]], use.type, mot.nsites[[1]],
                                               mot.bkgs[[1]], mot.pseudo[[1]],
                                               relative_entropy)
    if (use.custom) {
      p <- ggseqlogo(mot.mats[[1]], method = plot.method,
                     seq_type = seq_type, namespace = alph, ...) +
             ylab(yname)
    } else {
      p <- ggseqlogo(mot.mats[[1]], method = plot.method,
                     seq_type = seq_type, ...) + ylab(yname)
    }
    return(p)
  }

  res <- view_motifs_prep(mot.mats, method, tryRC, min.overlap, min.mean.ic,
                          min.position.ic, mot.bkgs, relative_entropy,
                          normalise.scores, alph)
  which.rc <- res[[1]]
  mots <- res[[2]]

  mots <- mapply(function(x1, x2, x3, x4)
                   convert_mat_type_from_ppm(x1, use.type, x2, x3, x4, relative_entropy),
                 mots, mot.nsites, mot.bkgs, mot.pseudo, SIMPLIFY = FALSE)

  for (i in seq_along(which.rc)) {
    if (which.rc[i]) mot.names[i] <- paste(mot.names[i], "[RC]")
  }
  names(mots) <- mot.names

  if (use.custom) {
    ggplot() + geom_logo(mots, method = plot.method, seq_type = seq_type, 
                         namespace = alph, ...) +
      theme_logo() +
      facet_wrap(~seq_group, ncol = 1) + ylab(yname)
  } else {
    ggplot() + geom_logo(mots, method = plot.method, seq_type = seq_type, ...) +
      theme_logo() +
      facet_wrap(~seq_group, ncol = 1) + ylab(yname)
  }

}

convert_mat_type_from_ppm <- function(mot.mat, type, nsites, bkg, pseudocount,
                                      relative_entropy) {

  which.zero <- apply(mot.mat, 2, function(x) all(x == 0))

  if (length(nsites) == 0) nsites <- 100
  if (length(pseudocount) == 0) pseudocount <- 1

  mot.mat[, which.zero] <- switch(type,
                             "PCM" = apply(mot.mat[, which.zero, drop = FALSE], 2,
                                           ppm_to_pcmC, nsites = nsites),
                             "PWM" = apply(mot.mat[, which.zero, drop = FALSE], 2,
                                           ppm_to_pwmC, bkg = bkg,
                                           pseudocount = pseudocount,
                                           nsites = nsites),
                             "ICM" = apply(mot.mat[, which.zero, drop = FALSE], 2,
                                           ppm_to_icmC, bkg = bkg,
                                           relative_entropy = relative_entropy)
                           )

  mot.mat

}
