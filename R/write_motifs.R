#' Export motifs in universalmotif format.
#'
#' Write motifs as universalmotif objects to file. For optimal storage of
#' `universalmotif` class motifs, consider using [saveRDS()] and
#' [readRDS()]. The `universalmotif` format will not be documented,
#' as realistically the need to generate these manually/elsewhere should
#' be non-existent.
#'
#' @param minimal `logical(1)` Only write essential motif information.
#' @param multifreq `logical(1)` Write `multifreq` slot, if present.
#' @param progress `logical(1)` Show progress.
#'
#' @return `NULL`, invisibly.
#'
#' @family write_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams write_jaspar
#' @export
write_motifs <- function(motifs, file, minimal = FALSE, multifreq = TRUE,
                         progress = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file), 1, FALSE, "character")
  logi_check <- check_fun_params(list(minimal = args$minimal,
                                      multifreq = args$multifreq,
                                      progress = args$progress),
                                 numeric(), logical(), "logical")
  all_checks <- c(char_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  mots <- lapply_(motifs, write_motifs2_single, minimal = minimal,
                  multifreq = multifreq, PB = progress)
  names(mots) <- paste0("MOTIF", seq_along(mots))

  mots <- as.yaml(mots, indent = 4, precision = 12)

  mots <- collapse_cpp(c("# universalmotif version ",
                         packageDescription("universalmotif")$Version,
                         "\n#\n", mots))

  writeLines(mots, con <- file(file))
  close(con)

  invisible(NULL)

}

write_motifs2_single <- function(motif, minimal, multifreq) {

  if (minimal) {

    lmot <- list(name = motif@name,
                 alphabet = motif@alphabet,
                 type = motif@type,
                 pseudocount = motif@pseudocount,
                 bkg = as.list(motif@bkg),
                 strand = motif@strand)

  } else {

    lmot <- list(name = motif@name,
                 altname = motif@altname,
                 family = motif@family,
                 organism = motif@organism,
                 alphabet = motif@alphabet,
                 type = motif@type,
                 icscore = round(motif@icscore, 2),
                 nsites = motif@nsites,
                 pseudocount = motif@pseudocount,
                 bkg = as.list(motif@bkg),
                 bkgsites = motif@bkgsites,
                 consensus = motif@consensus,
                 strand = motif@strand,
                 pval = motif@pval,
                 qval = motif@qval,
                 eval = motif@eval,
                 extrainfo = as.list(motif@extrainfo))

    to.keep <- vapply(lmot, length, integer(1)) > 0
    lmot <- lmot[to.keep]

  }


  mot.t <- t(motif@motif)
  mot.bycol <- apply(mot.t, 1, format_pos)
  mot.bycol <- as.list(unname(mot.bycol))

  lmot$motif <- mot.bycol

  if (multifreq && length(motif@multifreq) > 0) {

    mult <- vector("list", length(motif@multifreq))
    for (i in seq_along(mult)) {
      mult.t <- t(motif@multifreq[[i]])
      mult.l <- apply(mult.t, 1, format_pos)
      mult.l <- as.list(unname(mult.l))
      mult[[i]] <- mult.l
    }

    names(mult) <- names(motif@multifreq)
    lmot$multifreq <- mult

  }

  lmot

}

format_pos <- function(x) {
  paste0(formatC(x, format = "e", digits = 3), collapse = " ")
}
