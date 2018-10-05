#' Merge motifs.
#'
#' Aligns the motifs using [compare_motifs()], then averages the
#' motif PPMs. Currently the `multifreq` slot, if filled in any of the motifs,
#' will be dropped.
#'
#' @return A single motif object. See [convert_motifs()] for
#'    available formats.
#'
#' @examples
#' if (requireNamespace("MotifDb", quietly = TRUE)) {
#'   library(MotifDb)
#'   merged.motif <- merge_motifs(MotifDb[1:5])
#' }
#'
#' @seealso [compare_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams compare_motifs
#' @export
merge_motifs <- function(motifs, method = "MPCC", use.type = "PPM",
                         min.overlap = 6, min.mean.ic = 0.5, tryRC = TRUE,
                         relative_entropy = FALSE, normalise.scores = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(method = args$method, use.type = args$use.type),
                                 numeric(), logical(), "character")
  num_check <- check_fun_params(list(min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic),
                                numeric(), logical(), "numeric")
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (use.type %in% c("PCM", "PWM") && method %in% c("EUCL", "NEUCL", "KL")) {
    stop("Method '", method, "' is not supported for type '", use.type, "'")
  }

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  motifs <- convert_type(motifs, use.type, relative_entropy = relative_entropy)

  mot <- merge_mot_list(motifs, tryRC, min.overlap, min.mean.ic, method,
                        relative_entropy, normalise.scores)

  mot <- .internal_convert(mot, unique(CLASS_IN))
  mot

}

#' merge_mot_pair
#'
#' Merge two motifs.
#'
#' @param mot1 Motif matrix 1.
#' @param mot2 Motif matrix 2.
#' @param weight1 Weighing factor for motif 1.
#' @param weight2 Weighing factor for motif 2.
#' @param ic1 Positional ICs for motif 1.
#' @param ic2 Positional ICs for motif 2.
#' @param tryRC Logical.
#' @param min.overlap Minimum overlap.
#' @param min.mean.ic Minimum mean IC of alignment.
#' @param method Comparison metric.
#' @param relative_entropy Logical.
#' @param normalise.scores Logical.
#'
#' @noRd
merge_mot_pair <- function(mot1, mot2, weight1, weight2, ic1, ic2, tryRC,
                           min.overlap, min.mean.ic, method, relative_entropy,
                           normalise.scores) {

  out <- merge_motifs_internal(mot1, mot2, method, min.overlap, tryRC, ic1, ic2,
                               min.mean.ic, weight1, weight2, normalise.scores)

  matrix(out[!is.na(out)], nrow = nrow(out))


}

#' merge_mot_list
#'
#' Merge a list of motifs, pairwise.
#'
#' @param motifs List of universalmotif motifs.
#' @param tryRC Logical.
#' @param min.overlap Minimum overlap.
#' @param min.mean.ic Minimum mean IC for alignment.
#' @param method Comparison metric.
#' @param relative_entropy Logical.
#' @param normalise.scores Logical.
#'
#' @noRd
merge_mot_list <- function(motifs, tryRC, min.overlap, min.mean.ic, method,
                           relative_entropy, normalise.scores) {

  mot.names <- vapply(motifs, function(x) x["name"], character(1))
  mot.altnames <- do.call(c, sapply(motifs, function(x) x["altname"], simplify = FALSE))
  mot.families <- do.call(c, sapply(motifs, function(x) x["family"], simplify = FALSE))
  mot.orgs <- do.call(c, sapply(motifs, function(x) x["organism"], simplify = FALSE))
  mot.bkgsites <- do.call(c, sapply(motifs, function(x) x["bkgsites"], simplify = FALSE))
  mot.strands <- vapply(motifs, function(x) x["strand"], character(1))
  mot.extrainfo <- lapply(motifs, function(x) x["extrainfo"])

  mot.mats <- lapply(motifs, function(x) x["motif"])

  alph <- motifs[[1]]["alphabet"]

  mot.ic1 <- .pos_iscscores(motifs[[1]], mot.mats[[1]], relative_entropy)
  mot.ic2 <- .pos_iscscores(motifs[[2]], mot.mats[[2]], relative_entropy)
  new.mat <- merge_mot_pair(mot.mats[[1]], mot.mats[[2]], 1, 1, mot.ic1,
                            mot.ic2, tryRC, min.overlap, min.mean.ic, method,
                            relative_entropy, normalise.scores)
  bkg.1 <- motifs[[1]]["bkg"]
  bkg.2 <- motifs[[2]]["bkg"]
  bkg.new <- vapply(seq_along(bkg.1), function(x) mean(c(bkg.1[x], bkg.2[x])),
                    numeric(1))
  nsites.1 <- motifs[[1]]["nsites"]
  nsites.2 <- motifs[[2]]["nsites"]
  if (length(nsites.1) == 0 && length(nsites.2) == 0) {
    nsites.new <- numeric()
  } else nsites.new <- max(c(nsites.1, nsites.2))
  pseudo.1 <- motifs[[1]]["pseudocount"]
  pseudo.2 <- motifs[[2]]["pseudocount"]
  pseudo.new <- mean(c(pseudo.1, pseudo.2))

  mot.new <- create_motif(new.mat, alphabet = alph, pseudocount = pseudo.new,
                          bkg = bkg.new, nsites = nsites.new)

  if (length(motifs) > 2) {
    add.weight <- 2
    for (i in seq(3, length(motifs))) {
      mot.ic <- .pos_iscscores(motifs[[i]], mot.mats[[i]], relative_entropy)
      new.ic <- .pos_iscscores(mot.new, mot.new["motif"], relative_entropy)
      new.mat <- merge_mot_pair(new.mat, mot.mats[[i]], add.weight, 1,
                                new.ic, mot.ic, tryRC, min.overlap,
                                min.mean.ic, method, relative_entropy,
                                normalise.scores)
      bkg.1 <- motifs[[i]]["bkg"]
      bkg.2 <- mot.new["bkg"]
      bkg.new <- vapply(seq_along(bkg.1), function(x) mean(c(bkg.1[x], bkg.2[x])),
                      numeric(1))
      nsites.1 <- motifs[[i]]["nsites"]
      nsites.2 <- mot.new["nsites"]
      if (length(nsites.1) == 0 && length(nsites.2) == 0) {
        nsites.new <- numeric()
      } else nsites.new <- max(c(nsites.1, nsites.2))
      pseudo.1 <- motifs[[i]]["pseudocount"]
      pseudo.2 <- mot.new["pseudocount"]
      pseudo.new <- mean(c(pseudo.1, pseudo.2))
      mot.new <- create_motif(new.mat, alphabet = alph, pseudocount = pseudo.new,
                              bkg = bkg.new, nsites = nsites.new)
      add.weight <- add.weight + 1
    }
  }

  new.name <- paste(mot.names, collapse = "/")
  new.altname <- paste(mot.altnames, collapse = "/")
  if (nchar(new.altname) == 0) new.altname <- character(0)
  new.family <- paste(mot.families, collapse = "/")
  if (nchar(new.family) == 0) new.family <- character(0) 
  new.organism <- paste(mot.orgs, collapse = "/")
  if (nchar(new.organism) == 0) new.organism <- character(0)
  if (length(mot.bkgsites) > 1) {
    new.bkgsites <- max(mot.bkgsites)
  } else new.bkgsites <- numeric(0)
  if (length(unique(mot.strands)) > 1) {
    new.strand <- "+-" 
  } else new.strand <- unique(mot.strands)
  new.extrainfo <- do.call(c, mot.extrainfo)

  mot.new["name"] <- new.name
  mot.new["altname"] <- new.altname
  mot.new["family"] <- new.family
  mot.new["organism"] <- new.organism
  mot.new["bkgsites"] <- new.bkgsites
  mot.new["strand"] <- new.strand
  mot.new["extrainfo"] <- new.extrainfo

  mot.new

}
