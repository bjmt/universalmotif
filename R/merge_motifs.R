#' Merge motifs.
#'
#' Aligns the motifs using [compare_motifs()], then averages the
#' motif PPMs. Currently the `multifreq` slot, if filled in any of the motifs,
#' will be dropped. Only 0-order background probabilities will be kept.
#'
#' @return A single motif object. See [convert_motifs()] for
#'    available formats.
#'
#' @examples
#' \dontrun{
#' library(MotifDb)
#' merged.motif <- merge_motifs(MotifDb[1:5])
#' }
#'
#' @seealso [compare_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams compare_motifs
#' @export
merge_motifs <- function(motifs, method = "MPCC", use.type = "PPM",
                         min.overlap = 6, min.mean.ic = 0.25, tryRC = TRUE,
                         relative_entropy = FALSE, normalise.scores = FALSE,
                         min.position.ic = 0) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!method %in% COMPARE_METRICS) {
    method_check <- paste0(" * Incorrect 'method': expected `PCC`, `MPCC`, `EUCL`,",
                           " `MEUCL`, `SW`, `MSW`, `KL` or `MKL`; got `",
                           method, "`")
    method_check <- wmsg2(method_check, 4, 2)
    all_checks <- c(all_checks, method_check)
  }
  char_check <- check_fun_params(list(method = args$method),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic),
                                numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(all_checks, char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (use.type != "PPM")
    stop("currently `use.type = \"PPM\"` is the only acceptable option [use.type=",
         use.type, "]")

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  motifs <- convert_type_internal(motifs, "PPM")

  mot <- merge_motifs_all(motifs, method, tryRC, min.overlap, min.mean.ic,
                          min.position.ic, relative_entropy, normalise.scores)

  mot <- .internal_convert(mot, unique(CLASS_IN))
  mot

}

merge_motifs_all <- function(motifs, method, tryRC, min.overlap, min.mean.ic,
                             min.position.ic, relative_entropy, normalise.scores) {

  alph <- unique(vapply(motifs, function(x) x@alphabet, character(1)))
  if (length(alph) > 1) stop("all motifs must have the same alphabet")

  alph2 <- switch(alph, "DNA" = DNA_BASES, "RNA" = RNA_BASES, "AA" = AA_STANDARD,
                  sort_unique_cpp(safeExplode(alph)))

  mot.mats <- lapply(motifs, function(x) x@motif)
  mot.bkgs <- lapply(motifs, function(x) x@bkg[seq_along(alph2)])

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  mot.altnames <- do.call(c, sapply(motifs, function(x) x@altname, simplify = FALSE))
  mot.families <- unique(do.call(c, sapply(motifs, function(x) x@family, simplify = FALSE)))
  mot.orgs <- unique(do.call(c, sapply(motifs, function(x) x@organism, simplify = FALSE)))
  mot.bkgsites <- do.call(c, sapply(motifs, function(x) x@bkgsites, simplify = FALSE))
  mot.strands <- vapply(motifs, function(x) x@strand, character(1))
  mot.extrainfo <- lapply(motifs, function(x) x@extrainfo)
  mot.nsites <- do.call(c, sapply(motifs, function(x) x@nsites, simplify = FALSE))
  mot.pseudo <- vapply(motifs, function(x) x@pseudocount, numeric(1))
  mot.pvals <- do.call(c, sapply(motifs, function(x) x@pval, simplify = FALSE))
  mot.qvals <- do.call(c, sapply(motifs, function(x) x@qval, simplify = FALSE))
  mot.evals <- do.call(c, sapply(motifs, function(x) x@eval, simplify = FALSE))

  ans <- merge_motifs_cpp(mot.mats, method, tryRC, min.overlap, min.mean.ic,
                          min.position.ic, mot.bkgs, relative_entropy,
                          normalise.scores)

  new.name <- paste0(mot.names, collapse = "/")
  new.altname <- paste0(mot.altnames, collapse = "/")
  if (nchar(new.altname) == 0) new.altname <- character(0)
  new.family <- paste0(mot.families, collapse = "/")
  if (nchar(new.family) == 0) new.family <- character(0)
  new.organism <- paste0(mot.orgs, collapse = "/")
  if (nchar(new.organism) == 0) new.organism <- character(0)
  if (length(mot.bkgsites) > 1) {
    new.bkgsites <- max(mot.bkgsites)
  } else new.bkgsites <- numeric(0)
  if (length(unique(mot.strands)) > 1) {
    new.strand <- "+-"
  } else new.strand <- unique(mot.strands)
  new.extrainfo <- unique(do.call(c, mot.extrainfo))
  new.pseudo <- ifelse(length(mot.pseudo) > 0, mean(mot.pseudo, na.rm = TRUE), numeric())
  new.nsites <- ifelse(length(mot.nsites) > 0, max(mot.nsites, na.rm = TRUE), numeric())
  new.pvals <- ifelse(length(mot.pvals) > 0, max(mot.pvals, na.rm = TRUE), numeric())
  new.qvals <- ifelse(length(mot.qvals) > 0, max(mot.qvals, na.rm = TRUE), numeric())
  new.evals <- ifelse(length(mot.evals) > 0, max(mot.evals, na.rm = TRUE), numeric())

  mot <- universalmotif_cpp(motif = ans[[1]], name = new.name, altname = new.altname, 
                            family = new.family, organism = new.organism,
                            alphabet = alph, type = "PPM", nsites = new.nsites,
                            pseudocount = new.pseudo, bkg = ans[[2]],
                            bkgsites = new.bkgsites, strand = new.strand,
                            pval = new.pvals, qval = new.qvals, eval = new.evals,
                            extrainfo = new.extrainfo)

  validObject_universalmotif(mot)

  mot
 
}

#-------------------------------------------------------------------------------

merge_motifs_v1 <- function(motifs, method = "MPCC", use.type = "PPM",
                            min.overlap = 6, min.mean.ic = 0.5, tryRC = TRUE,
                            relative_entropy = FALSE, normalise.scores = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!method %in% COMPARE_METRICS) {
    method_check <- paste0(" * Incorrect 'method': expected `PCC`, `MPCC`, `EUCL`,",
                           " `MEUCL`, `SW`, `MSW`, `KL` or `MKL`; got `",
                           method, "`")
    method_check <- wmsg2(method_check, 4, 2)
    all_checks <- c(all_checks, method_check)
  }
  if (!use.type %in% c("PPM", "ICM")) {
    use.type_check <- paste0(" * Incorrect 'use.type': expected `PPM` or `ICM`; got `",
                             use.type, "`")
    use.type_check <- wmsg2(use.type_check, 4, 2)
    all_checks <- c(all_checks, use.type_check)
  }
  char_check <- check_fun_params(list(method = args$method, use.type = args$use.type),
                                 numeric(), logical(), TYPE_CHAR)
  num_check <- check_fun_params(list(min.overlap = args$min.overlap,
                                     min.mean.ic = args$min.mean.ic),
                                numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(tryRC = args$tryRC,
                                      relative_entropy = args$relative_entropy,
                                      normalise.scores = args$normalise.scores),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(all_checks, char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  motifs <- convert_type_internal(motifs, use.type, relative_entropy = relative_entropy)

  mot <- merge_mot_list(motifs, tryRC, min.overlap, min.mean.ic, method,
                        relative_entropy, normalise.scores)

  mot <- .internal_convert(mot, unique(CLASS_IN))
  mot

}
merge_mot_pair <- function(mot1, mot2, weight1, weight2, ic1, ic2, tryRC,
                           min.overlap, min.mean.ic, method, relative_entropy,
                           normalise.scores) {

  out <- merge_motifs_internal(mot1, mot2, method, min.overlap, tryRC, ic1, ic2,
                               min.mean.ic, weight1, weight2, normalise.scores)

  matrix(out[!is.na(out)], nrow = nrow(out))


}

merge_mot_list <- function(motifs, tryRC, min.overlap, min.mean.ic, method,
                           relative_entropy, normalise.scores) {

  mot.names <- vapply(motifs, function(x) x@name, character(1))
  mot.altnames <- do.call(c, sapply(motifs, function(x) x@altname, simplify = FALSE))
  mot.families <- unique(do.call(c, sapply(motifs, function(x) x@family, simplify = FALSE)))
  mot.orgs <- unique(do.call(c, sapply(motifs, function(x) x@organism, simplify = FALSE)))
  mot.bkgsites <- do.call(c, sapply(motifs, function(x) x@bkgsites, simplify = FALSE))
  mot.strands <- vapply(motifs, function(x) x@strand, character(1))
  mot.extrainfo <- lapply(motifs, function(x) x@extrainfo)

  mot.mats <- lapply(motifs, function(x) x@motif)

  alph <- motifs[[1]]@alphabet

  mot.ic1 <- .pos_iscscores(motifs[[1]], mot.mats[[1]], relative_entropy)
  mot.ic2 <- .pos_iscscores(motifs[[2]], mot.mats[[2]], relative_entropy)
  new.mat <- merge_mot_pair(mot.mats[[1]], mot.mats[[2]], 1, 1, mot.ic1,
                            mot.ic2, tryRC, min.overlap, min.mean.ic, method,
                            relative_entropy, normalise.scores)
  bkg.1 <- motifs[[1]]@bkg[rownames(motifs[[1]]@motif)]
  bkg.2 <- motifs[[2]]@bkg[rownames(motifs[[2]]@motif)]
  bkg.new <- vapply(seq_along(bkg.1), function(x) mean(c(bkg.1[x], bkg.2[x])),
                    numeric(1))
  nsites.1 <- motifs[[1]]@nsites
  nsites.2 <- motifs[[2]]@nsites
  if (length(nsites.1) == 0 && length(nsites.2) == 0) {
    nsites.new <- numeric()
  } else nsites.new <- max(c(nsites.1, nsites.2))
  pseudo.1 <- motifs[[1]]@pseudocount
  pseudo.2 <- motifs[[2]]@pseudocount
  pseudo.new <- mean(c(pseudo.1, pseudo.2))

  mot.new <- create_motif(new.mat, alphabet = alph, pseudocount = pseudo.new,
                          bkg = bkg.new, nsites = nsites.new)

  if (length(motifs) > 2) {
    add.weight <- 2
    for (i in seq(3, length(motifs))) {
      mot.ic <- .pos_iscscores(motifs[[i]], mot.mats[[i]], relative_entropy)
      new.ic <- .pos_iscscores(mot.new, mot.new@motif, relative_entropy)
      new.mat <- merge_mot_pair(new.mat, mot.mats[[i]], add.weight, 1,
                                new.ic, mot.ic, tryRC, min.overlap,
                                min.mean.ic, method, relative_entropy,
                                normalise.scores)
      bkg.1 <- motifs[[i]]@bkg[rownames(motifs[[i]]@motif)]
      bkg.2 <- mot.new@bkg[rownames(mot.new@motif)]
      bkg.new <- vapply(seq_along(bkg.1), function(x) mean(c(bkg.1[x], bkg.2[x])),
                      numeric(1))
      nsites.1 <- motifs[[i]]@nsites
      nsites.2 <- mot.new@nsites
      if (length(nsites.1) == 0 && length(nsites.2) == 0) {
        nsites.new <- numeric()
      } else nsites.new <- max(c(nsites.1, nsites.2))
      pseudo.1 <- motifs[[i]]@pseudocount
      pseudo.2 <- mot.new@pseudocount
      pseudo.new <- mean(c(pseudo.1, pseudo.2))
      mot.new <- create_motif(new.mat, alphabet = alph, pseudocount = pseudo.new,
                              bkg = bkg.new, nsites = nsites.new)
      add.weight <- add.weight + 1
    }
  }

  new.name <- paste0(mot.names, collapse = "/")
  new.altname <- paste0(mot.altnames, collapse = "/")
  if (nchar(new.altname) == 0) new.altname <- character(0)
  new.family <- paste0(mot.families, collapse = "/")
  if (nchar(new.family) == 0) new.family <- character(0) 
  new.organism <- paste0(mot.orgs, collapse = "/")
  if (nchar(new.organism) == 0) new.organism <- character(0)
  if (length(mot.bkgsites) > 1) {
    new.bkgsites <- max(mot.bkgsites)
  } else new.bkgsites <- numeric(0)
  if (length(unique(mot.strands)) > 1) {
    new.strand <- "+-" 
  } else new.strand <- unique(mot.strands)
  new.extrainfo <- do.call(c, mot.extrainfo)

  mot.new@name <- new.name
  mot.new@altname <- new.altname
  mot.new@family <- new.family
  mot.new@organism <- new.organism
  mot.new@bkgsites <- new.bkgsites
  mot.new@strand <- new.strand
  mot.new@extrainfo <- new.extrainfo

  validObject_universalmotif(mot.new)
  mot.new

}

.pos_iscscores <- function(motif, mot.mats, relative = FALSE) {

  bkg <- motif@bkg[rownames(motif@motif)]
  pseudo <- motif@pseudocount
  nsites <- motif@nsites
  if (length(nsites) == 0) nsites <- 100
  apply(mot.mats, 2, function(x) position_icscoreC(x, bkg, "PPM", pseudo,
                                                   nsites, relative))

}
