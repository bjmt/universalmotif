#' Export motifs in universalmotif format.
#'
#' Write motifs as universalmotif objects to file. For optimal storage of
#' `universalmotif` class motifs, consider using [saveRDS()] and
#' [readRDS()]. Currently the `universalmotif` format is YAML-based, but
#' this is subject to change.
#'
#' @param minimal `logical(1)` Only write essential motif information.
#' @param multifreq `logical(1)` Write `multifreq` slot, if present.
#' @param progress `logical(1)` Show progress.
#' @param overwrite `logical(1)` Overwrite existing file.
#' @param append `logical(1)` Add to an existing motif file. Package version
#'    in existing motif file must be greater than 1.2.0.
#' @param BP `logical(1)` Allows for the use of \pkg{BiocParallel} within
#'    [write_motifs()]. See [BiocParallel::register()] to change the
#'    default backend.
#'
#' @return `NULL`, invisibly.
#'
#' @family write_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams write_jaspar
#' @export
write_motifs <- function(motifs, file, minimal = FALSE, multifreq = TRUE,
                         progress = FALSE, overwrite = FALSE,
                         append = FALSE, BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file), 1, FALSE, TYPE_CHAR)
  logi_check <- check_fun_params(list(minimal = args$minimal,
                                      multifreq = args$multifreq,
                                      progress = args$progress,
                                      overwrite = args$overwrite,
                                      BP = args$BP),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(char_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (file.exists(file) && !overwrite && !append)
    stop(wmsg("Existing file found, set `overwrite = TRUE` to continue."))

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  if (file.exists(file) && append) {

    old <- readLines(con <- file(file), n = 1); close(con)

    if (substr(old[1], 1, 24) == "# universalmotif version") {

      file.version <- strsplit(old[1], " ")[[1]][4]

      if (file.version < "1.2.0")
        stop(wmsg("To append motifs to an existing file, file version",
                  " must be at least 1.2.0"))

    } else {

      stop(wmsg("In order to append motifs, existing file must contain package version"))

    }

  }

  mots <- lapply_(motifs, write_motifs2_single, minimal = minimal,
                  multifreq = multifreq, PB = progress, BP = BP)
  mots <- unlist(mots)

  if (append) {

    mots <- collapse_cpp(mots)
    mots <- substr(mots, 2, nchar(mots))
    cat(mots, file = file, append = TRUE)

  } else {

    mots <- collapse_cpp(c("# universalmotif version ",
                           packageDescription("universalmotif")$Version,
                           "\n", mots))
    writeLines(mots, con <- file(file)); close(con)

  }

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

  collapse_cpp(c("\n---\n", as.yaml(lmot)))

}

format_pos <- function(x) {
  paste0(formatC(x, format = "e", digits = 3), collapse = " ")
}
