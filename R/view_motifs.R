#' Plot motif logos.
#'
#' Show sequence logo. If given a list of more than one motif, then the motifs
#' are aligned with the first in the list.
#'
#' @param motifs See \code{\link{convert_motifs}} for acceptable motif formats.
#' @param use.type \code{character(1)} One of \code{c('PCM', 'PPM', 'PWM', 'ICM')}.
#' @param method \code{character(1)} One of \code{c('Pearson', 'Euclidean', 'KL')}.
#' @param tryRC \code{logical(1)} Check if motif reverse complement leads to a
#'    better alignment.
#' @param min.overlap \code{numeric(1)} Minimum alignment overlap between
#'    motifs. If \code{min.overlap < 1}, this represents the minimum fraction
#'    between the two motifs during alignment.
#' @param min.mean.ic \code{numeric(1)} Minimum information content between the
#'    two motifs for an alignment to be scored. This helps prevent scoring
#'    alignments between low information content regions of two motifs.
#' @param relative_entropy \code{logical(1)} For ICM calculation. See
#'    \code{\link{convert_type}}.
#' @param normalise.scores \code{logical(1)} Favour alignments which leave fewer
#'    unaligned positions.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#' @param ... Addtional options for \code{\link[ggseqlogo]{geom_logo}}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Since the \code{\link[ggseqlogo]{ggseqlogo}} package can only plot individual
#' characters and not strings, plotting the \code{multifreq} slot is not
#' supported. See the \code{examples} section for plotting the \code{multifreq}
#' slot using the 'Logolas' package.
#'
#' @examples
#' ## plotting multifreq motifs:
#' if (requireNamespace("Logolas", quietly = TRUE)) {
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
#' @seealso \code{\link{compare_motifs}}
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
view_motifs <- function(motifs, use.type = "ICM", method = "Pearson",
                        tryRC = TRUE, min.overlap = 6,
                        min.mean.ic = 0.5, relative_entropy = FALSE,
                        normalise.scores = TRUE,
                        BPPARAM = SerialParam(), ...) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(use.type = args$use.type,
                                      method = args$method),
                                 numeric(), logical(), "character")
  num_check <- check_fun_params(list(min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic),
                                numeric(), logical(), "numeric")
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), "logical")
  s4_check <- check_fun_params(list(BPPARAM = args$BPPARAM),
                               numeric(), FALSE, "S4")
  all_checks <- c(char_check, num_check, logi_check, s4_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, use.type, relative_entropy = relative_entropy)
  if (!is.list(motifs)) motifs <- list(motifs)

  if (use.type == "ICM" && !relative_entropy) {
    plot.method <- "bits"
    yname <- "bits"
  } else if (use.type == "PPM") {
    plot.method <- "prob"
    yname <- "probability"
  } else if (use.type == "ICM") {
    plot.method <- "custom"
    yname <- "bits"
  } else if (use.type == "PWM") {
    plot.method <- "custom"
    yname <- "logodds"
    if (length(motifs) > 1 && method == "KL") {
      stop("cannot use method 'KL' with 'PWM' matrices")
    }
  } else if (use.type == "PCM") {
    plot.method <- "custom"
    yname <- "counts"
  } else stop("'use.type' must be one of 'PCM', 'PPM', 'PWM', 'ICM'")

  mot.names <- vapply(motifs, function(x) x["name"], character(1))

  mot.mats <- lapply(motifs, function(x) x["motif"])

  mot.alph <- unique(vapply(motifs, function(x) x["alphabet"], character(1)))
  if (length(mot.alph) > 1) stop("can only have one alphabet")
  use.custom <- FALSE
  if (mot.alph == "DNA") {
    alph <- DNA_BASES
    seq_type <- "dna"
  } else if (mot.alph == "RNA") {
    alph <- RNA_BASES
    seq_type <- "rna"
  } else if (mot.alph == "AA") {
    alph <- AA_STANDARD
    seq_type <- "aa"
  } else if (mot.alph != "custom"){
    alph <- sort(strsplit(mot.alph, "")[[1]])
    use.custom <- TRUE
  } else{
    alph <- rownames(mot.mats[[1]])
    use.custom <- TRUE
  }

  if (length(motifs) == 1) {
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

  motifs.rc <- motif_rc(motifs)
  mot.mats.rc <- lapply(motifs.rc, function(x) x["motif"])

  mot.ics <- bpmapply(function(x, y) .pos_iscscores(x, y, relative_entropy),
                      motifs, mot.mats, BPPARAM = BPPARAM)
  mot.ics.rc <- lapply(mot.ics, rev)

  mot.alns <- bplapply(seq_along(mot.mats)[-1],
                       function(x) {
                         y <- merge_motifs_get_offset(mot.mats[[1]], mot.mats[[x]],
                                                      method, min.overlap,
                                                      mot.ics[[1]], mot.ics[[x]],
                                                      min.mean.ic, normalise.scores)
                         merge_add_cols(y)
                         y},
                       BPPARAM = BPPARAM)
  if (tryRC) {
    mot.alns.rc <- bplapply(seq_along(mot.mats)[-1],
                         function(x) {
                           y <- merge_motifs_get_offset(mot.mats[[1]], mot.mats.rc[[x]],
                                                        method, min.overlap,
                                                        mot.ics[[1]], mot.ics.rc[[x]],
                                                        min.mean.ic, normalise.scores)
                           merge_add_cols(y)
                           y},
                         BPPARAM = BPPARAM)
  } else mot.alns.rc <- NULL

  mot.alns <- lapply(mot.alns, function(x) {
                       x[1:2] <- fix_blank_pos(x[[1]], x[[2]]); x
                         })
  mot.alns.rc <- lapply(mot.alns.rc, function(x) {
                          x[1:2] <- fix_blank_pos(x[[1]], x[[2]]); x
                            })

  mots <- realign_all_mots(alph, mot.alns, mot.alns.rc)
  if (!is.null(mot.alns.rc)) {
    which.rc <- mots[[2]]
    mots <- mots[[1]]
    mot.names[which(which.rc) + 1] <- paste(mot.names[which(which.rc) + 1],
                                            "[RC]")
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

fix_blank_pos <- function(motif1, motif2) {

  na1 <- vapply(motif1[1, ], is.na, logical(1))
  na2 <- vapply(motif2[1, ], is.na, logical(1))

  tokeep <- rep(TRUE, ncol(motif1))

  for (i in seq_along(na1)) {
    if (na1[i] && na2[i]) {
      tokeep[i] <- FALSE
    } else break
  }
  for (i in rev(seq_along(na1))) {
    if (na1[i] && na2[i]) {
      tokeep[i] <- FALSE
    } else break
  }

  list(motif1[, tokeep], motif2[, tokeep])

}

realign_all_mots <- function(alph, mot.alns, mot.alns.rc) {

  mots.1 <- lapply(mot.alns, function(x) x$mot1_new)
  mots.2 <- lapply(mot.alns, function(x) x$mot2_new)
  num.rows <- nrow(mots.1[[1]])

  if (!is.null(mot.alns.rc)) {
    mots.1.rc <- lapply(mot.alns.rc, function(x) x$mot1_new)
    mots.2.rc <- lapply(mot.alns.rc, function(x) x$mot2_new)
    scores <- vapply(mot.alns, function(x) x$score, numeric(1))
    scores.rc <- vapply(mot.alns.rc, function(x) x$score, numeric(1))
    which.rc <- rep(FALSE, length(scores))
    for (i in seq_along(scores)) {
      if (scores.rc[i] > scores[i]) {
        mots.1[[i]] <- mots.1.rc[[i]]
        mots.2[[i]] <- mots.2.rc[[i]]
        which.rc[i] <- TRUE
      }
    }
  }

  mots.1.nas <- vapply(mots.1, count_leading_na, numeric(1))
  mots.2.nas <- vapply(mots.2, count_leading_na, numeric(1))
  mot1.highest <- order(mots.1.nas, decreasing = TRUE)[1]

  mots.2 <- lapply(seq_along(mots.2),
                   function(x) {
                     if (x == mot1.highest) return(mots.2[[x]])
                     cbind(matrix(rep(NA, num.rows * mots.1.nas[mot1.highest]),
                                  nrow = num.rows),
                           mots.2[[x]])
                   })

  mots <- c(mots.1[mot1.highest], mots.2)
  mot.lens <- vapply(mots, ncol, numeric(1))
  mots <- lapply(seq_along(mots),
                 function(x) {
                   if (ncol(mots[[x]]) == max(mot.lens)) return(mots[[x]])
                   cbind(mots[[x]],
                         matrix(rep(NA, num.rows * (max(mot.lens) - mot.lens[x])),
                                nrow = num.rows))
                 })
  mots <- lapply(mots, function(x) {rownames(x) <- alph; x})

  if (!is.null(mot.alns.rc)) {
    list(mots, which.rc)
  } else {
    mots
  }

}

count_leading_na <- function(mat) {

  mat.nas <- vapply(mat[1, ], is.na, logical(1))
  count.na <- 0
  for (i in seq_along(mat.nas)) {
    if (!mat.nas[i]) break
    count.na <- count.na + 1
  }

  count.na

}
