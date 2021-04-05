#' Export motifs in JASPAR format.
#'
#' Convert motifs to JASPAR format and write to file.
#' See \url{http://jaspar.genereg.net/}.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#' @param overwrite `logical(1)` Overwrite existing file.
#' @param append `logical(1)` Add to an existing file.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' transfac <- read_transfac(system.file("extdata", "transfac.txt",
#'                                     package = "universalmotif"))
#' write_jaspar(transfac, tempfile())
#'
#' @references
#'
#' Khan A, Fornes O, Stigliani A, Gheorghe M, Castro-Mondragon JA,
#' van der Lee R, Bessy A, Cheneby J, Kulkarni SR, Tan G, Baranasic
#' D, Arenillas DJ, Sandelin A, Vandepoele K, Lenhard B, Ballester B,
#' Wasserman WW, Parcy F, Mathelier A (2018). “JASPAR 2018: update of
#' the open-access database of transcription factor binding profiles
#' and its web framework.” *Nucleic Acids Research*, **46**, D260-D266.
#'
#' @family write_motifs
#' @seealso [read_jaspar()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
write_jaspar <- function(motifs, file, overwrite = FALSE, append = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file),
                                 1, FALSE, TYPE_CHAR)
  logi_check <- check_fun_params(list(overwrite = args$overwrite,
                                      append = args$append),
                                 c(1, 1), c(FALSE, FALSE), TYPE_LOGI)
  all_checks <- c(char_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (file.exists(file) && !overwrite && !append)
    stop(wmsg("Existing file found, set `overwrite = TRUE` to continue."))

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PCM")
  if (!is.list(motifs)) motifs <- list(motifs)

  .write_jaspar <- function(motifs) {

    lines_out <- vector()

    motif <- motifs
    if (length(motif@altname) > 0) {
      lines_out <- c(lines_out, paste0(">", motif@name, " ",
                                       motif@altname))
    } else {
      lines_out <- c(lines_out, paste0(">", motif@name))
    }

    alph <- switch(motif@alphabet, "DNA" = DNA_BASES, "RNA" = RNA_BASES,
                   "AA" = AA_STANDARD2, stop("unknown alphabet"))

    nsites <- motif@nsites
    if (length(nsites) == 0) nsites <- 100
    if (nsites > 99999) {
      width <- 6
    } else if (nsites > 9999) {
      width <- 5
    } else if (nsites > 999) {
      width <- 4
    } else if (nsites > 99) {
      width <- 3
    } else width <- 2

    for (j in seq_along(alph)) {
      p1 <- alph[j]
      p2 <- "["
      p3 <- as.numeric(motif@motif[j, ])
      p3 <- formatC(p3, width = width, format = "d")
      p3 <- paste0(p3, collapse = " ")
      p4 <- "]"
      lines_out <- c(lines_out, paste(p1, p2, p3, p4))
    }

    lines_out <- c(lines_out, "")

  }

  lines_out <- lapply(motifs, .write_jaspar)
  lines_out <- unlist(lines_out)

  if (append) {
    cat(lines_out, sep = "\n", file = file, append = TRUE)
  } else {
    writeLines(lines_out, con <- file(file)); close(con)
  }

  invisible(NULL)

}
