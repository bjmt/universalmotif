#' Export motifs in MEME format.
#'
#' Convert motifs to minimal MEME format and write to file.
#' See \url{http://meme-suite.org/doc/meme-format.html}.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#' @param version `numeric(1)` MEME version.
#' @param bkg `numeric` Background letter frequencies. If missing, will use background
#'   frequencies from motif objects (if they are identical); else
#'   background frequencies will be set to freq = 1/length(alphabet)
#' @param strand `character` If missing, will use strand from motif objects (if identical);
#'   otherwise will default to "+ -"
#' @param overwrite `logical(1)` Overwrite existing file.
#' @param append `logical(1)` Add to an existing file. Motifs will be written
#'    in minimal format, so it is recommended to only use this if the existing
#'    file is also a minimal MEME format file.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' transfac <- read_transfac(system.file("extdata", "transfac.txt",
#'                                     package = "universalmotif"))
#' write_meme(transfac, tempfile())
#'
#' @references
#'    \insertRef{meme}{universalmotif}
#'
#' @family write_motifs
#' @seealso [read_meme()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_meme <- function(motifs, file, version = 5, bkg, strand,
                       overwrite = FALSE, append = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file, strand = args$strand),
                                 numeric(), c(FALSE, TRUE), "character")
  num_check <- check_fun_params(list(version = args$version), 1, FALSE, "numeric")
  logi_check <- check_fun_params(list(overwrite = args$overwrite,
                                      append = args$append),
                                 c(1, 1), c(FALSE, FALSE), "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (file.exists(file) && !overwrite && !append)
    stop(wmsg("Existing file found, set `overwrite = TRUE` to continue."))

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PPM")
  if (!is.list(motifs)) motifs <- list(motifs)

  if (missing(strand)) {
    strand <- unique(vapply(motifs, function(x) x@strand, character(1)))
    if (length(strand) > 1) strand <- "+ -"
    if (strand == "+-" || strand == "-+") strand <- "+ -"
  }

  alph <- unique(vapply(motifs, function(x) x@alphabet, character(1)))
  if (length(alph) > 1) stop("all motifs must use the same alphabet")

  switch(alph,
    "DNA" = {
      alph.2 <- "ACGT"
      alph.length <- 4
    },
    "RNA" = {
      alph.2 <- "ACGU"
      alph.length <- 4
    },
    "AA" = {
      alph.2 <- "ACDEFGHIKLMNPQRSTVWY"
      alph.length <- 20
    },
    stop("unknown alphabet")
  )

  if (missing(bkg)) {
    bkg <- lapply(motifs, function(x) x@bkg[safeExplode(alph.2)])
    bkgtest <- lapply(bkg, function(x) all(x == bkg[[1]]))
    if (all(unlist(bkgtest))) {
      bkg <- bkg[[1]]
    } else {
      if (alph %in% c("DNA", "RNA")) bkg <- rep(0.25, 4)
      if (alph == "AA") bkg <- rep(1 / 20, 20)
    }
  }
  bkg.2 <- safeExplode(alph.2)

  lines_out <- c(paste("MEME version", version), "",
                 paste("ALPHABET=", alph.2), "",
                 paste("strands:", strand), "",
                 paste("Background letter frequencies"),
                 paste(bkg.2, bkg, collapse = " "), "")

  .write_meme <- function(motifs) {
  
    motif <- motifs
    lines_out <- paste("MOTIF", motif@name)
    nsites <- motif@nsites
    nsites <- ifelse(length(nsites) == 0, 100, nsites)
    eval <- motif@eval
    eval <- ifelse(length(eval) == 0, 0, eval)
    lines_out <- c(lines_out,
                   paste("letter-probability matrix:", "alength=", alph.length,
                         "w=", ncol(motif@motif), "nsites=", nsites,
                         "E=", eval))
    for (i in seq_len(ncol(motif@motif))) {
      pos <- motif@motif[, i]
      pos <- vapply(pos, function(x) format(x, nsmall = 6), character(1))
      pos <- paste("", pos, "", collapse = "")
      lines_out <- c(lines_out, pos)
    }

    c(lines_out, "")

  }

  if (append) {
    lines_out <- unlist(lapply(motifs, .write_meme))
    cat(lines_out, sep = "\n", file = file, append = TRUE)
  } else {
    lines_out <- c(lines_out, unlist(lapply(motifs, .write_meme)))
    writeLines(lines_out, con <- file(file)); close(con)
  }

  invisible(NULL)


}
