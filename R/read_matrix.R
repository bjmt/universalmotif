#' Import motifs from raw matrices.
#'
#' Import simply formatted motifs.
#'
#' @param file `character(1)` File name.
#' @param skip `numeric(1)` If not zero, will skip however many desired lines in the
#'    file before starting to read.
#' @param positions `character(1)` One of `c('columns', 'rows')`. Partial matching
#'    allowed. Indicate whether each
#'    position within a motif is represented as a row or a column in the file.
#' @param alphabet `character(1)` One of `c('DNA', 'RNA', 'AA')`,
#'    or a string of letters.
#' @param type `character(1)` One of `c('PCM', 'PPM', 'PWM', 'ICM')`.
#'    If missing will try and guess which one.
#' @param sep `character(1)` Indicates how individual motifs are separated. Set as
#'    `NULL` if there are no seperating lines between motifs (the default is to
#'    assume a blank line).
#' @param headers `logical(1)`, `character(1)` Indicating if and how to read names.
#' @param rownames `logical(1)` Are there alphabet letters present as rownames?
#' @param comment `NULL`, `character(1)` Character denoting lines to be considered
#'    comments.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @examples
#'    hocomoco <- system.file("extdata", "hocomoco.txt", package = "universalmotif")
#'    hocomoco <- read_matrix(hocomoco, headers = ">", positions = "rows")
#'
#' @family read_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
read_matrix <- function(file, skip = 0, type, positions = "columns",
                        alphabet = "DNA", sep = "", headers = TRUE,
                        rownames = FALSE, comment = NULL) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file, type = args$type,
                                      positions = args$positions,
                                      alphabet = args$alphabet, sep = args$sep),
                                 numeric(), c(FALSE, TRUE, FALSE, FALSE, TRUE),
                                 TYPE_CHAR)
  num_check <- check_fun_params(list(skip = args$skip), 1, FALSE, TYPE_NUM)
  logi_check <- check_fun_params(list(rownames = args$rownames),
                                 1, FALSE, TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  positions <- match.arg(positions, c("columns", "rows"))

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]

  if (!is.null(comment))
    raw_lines <- raw_lines[!grepl(collapse_cpp(c("^", comment)), raw_lines)]

  if (!is.null(sep)) {
    seperators <- which(raw_lines == sep)
    if (length(seperators) != 0) {
      if (seperators[length(seperators)] != length(raw_lines)) {
        seperators <- c(seperators, length(raw_lines) + 1)
      }
      if (seperators[1] == 1) {
        motif_stops <- seperators[-1] - 1
      } else motif_stops <- seperators - 1
    } else motif_stops <- length(raw_lines)
  }

  if (!isFALSE(headers)) {
    if (!isTRUE(headers)) {
      preheader <- headers
      headers <- grep(headers, raw_lines)
      motif_starts <- headers + 1
      headers <- raw_lines[headers]
      headers <- vapply(headers, function(x) sub(preheader, "", x),
                        character(1))
    } else {
      headers <- c(1, seperators + 1)
      if (headers[length(headers)] >= length(raw_lines))
        headers <- headers[-length(headers)]
      motif_starts <- headers + 1
      headers <- raw_lines[headers]
    }
  } else {
    if (!is.null(sep)) {
      motif_starts <- c(1, seperators + 1)
      if (motif_starts[length(motif_starts)] >= length(raw_lines))
        motif_starts <- motif_starts[-length(motif_starts)]
    } else {
      motif_starts <- 1
    }
  }

  if (is.null(sep)) {
    if (!isFALSE(headers)) {
      motif_stops <- c(motif_starts[-1] - 2, length(raw_lines))
    } else {
      motif_stops <- length(raw_lines)
    }
  }

  if (length(motif_starts) != length(motif_stops))
    stop("parsing error; try setting 'header' and 'sep'")
  motifs <- mapply(function(x, y) raw_lines[x:y],
                     motif_starts, motif_stops, SIMPLIFY = FALSE)
  if (rownames && positions == "columns") {
    motifs <- lapply(motifs, function(x) as.matrix(read.table(text = x,
                                                                row.names = 1)))
  } else if (rownames && positions == "rows") {
    motifs <- lapply(motifs, function(x) as.matrix(read.table(text = x,
                                                                header = TRUE)))
  } else {
    motifs <- lapply(motifs, function(x) as.matrix(read.table(text = x)))
  }

  if (positions == "rows") {
    motifs <- lapply(motifs, t)
  }

  motifs <- lapply(motifs, function(x) {rownames(x) <- NULL; x})
  if (!missing(type)) {
    motifs <- lapply(motifs, function(x) {
                      universalmotif_cpp(motif = x, type = type, alphabet = alphabet)
                     })
  } else {
    motifs <- lapply(motifs, function(x) universalmotif_cpp(motif = x,
                                                          alphabet = alphabet))
  }

  if (!isFALSE(headers)) {
    motifs <- mapply(function(x, y) {x@name <- y; x}, motifs, headers,
                       SIMPLIFY = FALSE)
  }

  motifs <- lapply(motifs, function(x) {validObject_universalmotif(x); x})

  if (length(motifs) == 1) motifs <- motifs[[1]]

  motifs

}
