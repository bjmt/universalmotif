#' Merge motifs.
#'
#' Aligns the motifs using [compare_motifs()], then averages the
#' motif PPMs. Currently the `multifreq` slot, if filled in any of the motifs,
#' will be dropped. Only 0-order background probabilities will be kept.
#' Motifs are merged one at a time, starting with the first entry in the
#' list.
#'
#' @return A single motif object. See [convert_motifs()] for
#'    available formats.
#'
#' @details
#'    See [compare_motifs()] for more info on comparison parameters.
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
    method_check <- paste0(" * Incorrect 'method': expected `PCC`, `MPCC`, `EUCL`, `MEUCL`",
                           ", `SW`, `MSW`, `KL`, `MKL`, `ALLR`, or `MALLR`; got `",
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
    stop(wmsg("deprecated, as `use.type = \"PPM\"` is now the only acceptable ",
              "option [use.type=", use.type, "]"))

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
  mot.bkgs <- get_bkgs(motifs)

  mot.names <- unique(vapply(motifs, function(x) x@name, character(1)))
  mot.altnames <- unique(do.call(c, sapply(motifs, function(x) x@altname, simplify = FALSE)))
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
                          normalise.scores, get_nsites(motifs))

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
  new.extrainfo <- new.extrainfo[!duplicated(new.extrainfo)]
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
