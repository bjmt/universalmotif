#' Export motifs as raw matrices.
#'
#' Write motifs as simple matrices with optional headers to file.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#' @param positions `character(1)` One of `c('columns', 'rows')`. Partial
#'   matching allowed.
#' @param rownames `logical(1)` Include alphabet letters as rownames.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM',
#'   'CWM')`. If missing, the motif is written with whatever type it
#'   is currently stored as, so a CWM motif (signed real values, no
#'   shared column sums ) round-trips through
#'   `write_matrix()` / `read_matrix(type = "CWM")` verbatim. Note
#'   that `type = "CWM"` here only works with
#'   CWM motifs; you cannot ask `write_matrix()` to convert a
#'   PPM (or any other type) into a CWM (the underlying
#'   [convert_type()] step refuses any `*` -> `CWM` direction).
#' @param sep `character(1)` Indicates how to separate individual motifs. Set as
#'    `NULL` to have no seperating lines between motifs (the default is to
#'    use a blank line).
#' @param headers `logical(1)`, `character(1)` Indicating if and how to write names.
#' @param overwrite `logical(1)` Overwrite existing file.
#' @param append `logical(1)` Add to an existing file.
#' @param digits `numeric(1)` Number of digits to use for motif positions.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' motif <- create_motif()
#' write_matrix(motif, tempfile(), headers = ">")
#'
#' ## CWM round-trip: write a CWM with its raw signed values, then
#' ## read it back with type = "CWM" preserved.
#' cwm.mat <- matrix(c( 0.1, -0.2,  0.8, -0.1,
#'                     -0.3,  0.9, -0.2,  0.1,
#'                      0.7, -0.1, -0.1, -0.4,
#'                     -0.2, -0.1, -0.1,  0.9),
#'                   nrow = 4, byrow = FALSE,
#'                   dimnames = list(c("A","C","G","T"), NULL))
#' cwm <- create_motif(cwm.mat, type = "CWM", name = "modisco_demo")
#' f <- tempfile()
#' write_matrix(cwm, f, headers = ">", rownames = TRUE)
#' re <- read_matrix(f, headers = ">", rownames = TRUE, type = "CWM")
#' re["type"]
#'
#' @family write_motifs
#' @seealso [read_matrix()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
write_matrix <- function(motifs, file, positions = "columns", rownames = FALSE,
                         type, sep = "", headers = TRUE, overwrite = FALSE,
                         append = FALSE, digits = 6) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file, positions = args$positions,
                                      type = args$type, sep = args$sep),
                                 numeric(), c(FALSE, FALSE, TRUE, TRUE),
                                 TYPE_CHAR)
  num_check <- check_fun_params(list(digits = args$digits),
                                 numeric(), logical(), TYPE_NUM)
  logi_check <- check_fun_params(list(rownames = args$rownames,
                                      overwrite = args$overwrite,
                                      append = args$append),
                                 c(1, 1, 1), c(FALSE, FALSE, FALSE), TYPE_LOGI)
  header_check <- character()
  if (!is.logical(headers) && !is.character(headers)) {
    header_check <- paste0(" * Incorrect type for 'headers': ",
                           "expected `logical` or `character`; got `",
                           class(headers), "`")
  } else if (length(headers) != 1) {
    header_check <- paste0(" * Incorrect vector length for 'headers': ",
                           "expected 1; got ", length(headers))
  }
  all_checks <- c(char_check, logi_check, header_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  positions <- match.arg(positions, c("columns", "rows"))

  if (file.exists(file) && !overwrite && !append)
    stop(wmsg("Existing file found, set `overwrite = TRUE` to continue."))

  motifs <- convert_motifs(motifs)
  if (!missing(type)) motifs <- convert_type_internal(motifs, type)
  if (!is.list(motifs)) motifs <- list(motifs)

  .write_matrix <- function(motifs, positions, rownames, sep, headers) {

    motif <- motifs

    lines_out <- vector()

    if (!isFALSE(headers)) {
      if (isTRUE(headers)) lines_out <- c(lines_out, motif@name) else {
        lines_out <- c(lines_out, paste0(headers, motif@name))
      }
    }

    switch(positions,
      "columns" = {
        for (i in seq_len(nrow(motif@motif))) {
          pos <- motif@motif[i, ]
          pos <- formatC(pos, format = "f", digits = digits)
          if (rownames) pos <- c(rownames(motif@motif)[i], pos)
          lines_out <- c(lines_out, paste0(pos, collapse = "\t"))
        }
      },
      "rows" = {
        if (rownames) lines_out <- c(lines_out, paste0(rownames(motif@motif),
                                                       collapse = "\t"))
        for (i in seq_len(ncol(motif@motif))) {
          pos <- motif@motif[, i]
          pos <- formatC(pos, format = "f", digits = digits)
          pos <- as.character(pos)
          lines_out <- c(lines_out, paste0(pos, collapse = "\t"))
        }
      }
    )

    c(lines_out, sep)

  }

  lines_final <- lapply(motifs, .write_matrix, positions = positions,
                          rownames = rownames, sep = sep,
                          headers = headers)  # not working??
  lines_final <- unlist(lines_final)

  if (append) {
    cat(lines_final, sep = "\n", file = file, append = TRUE)
  } else {
    con <- file(file)
    on.exit(close(con), add = TRUE)
    writeLines(lines_final, con)
  }

  invisible(NULL)

}
