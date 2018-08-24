#' Import motifs from raw matrices.
#'
#' Import motifs formatted simply, as matrices.
#'
#' @param file \code{character(1)} File name.
#' @param skip \code{numeric(1)} If not zero, will skip however many desired lines in the
#'    file before starting to read.
#' @param positions \code{character(1)} One of \code{c('columns', 'rows')}.
#'    Indicate whether each
#'    position within a motif is represented as a row or a column in the file.
#' @param alphabet \code{character(1)} One of \code{c('DNA', 'RNA', 'AA')},
#'    or a string of letters.
#' @param type \code{character(1)} One of \code{c('PCM', 'PPM', 'PWM', 'ICM')}.
#'    If missing will try and guess which one.
#' @param sep \code{character(1)} Indicates how individual motifs are seperated.
#' @param headers \code{logical(1)}, \code{character(1)} Indicating if and how to read names. 
#' @param rownames \code{logical(1)} Are there alphabet letters present as rownames?
#'
#' @return \code{list} \linkS4class{universalmotif} objects.
#'
#' @examples
#'    hocomoco <- system.file("extdata", "hocomoco.txt", package = "universalmotif")
#'    hocomoco <- read_matrix(hocomoco, headers = ">", positions = "rows")
#'
#' @family read_motifs
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_matrix <- function(file, skip = 0, type, positions = "columns",
                        alphabet = "DNA", sep = "", headers = FALSE,
                        rownames = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file, type = args$type,
                                      positions = args$positions,
                                      alphabet = args$alphabet, sep = args$sep),
                                 numeric(), c(FALSE, TRUE, FALSE, FALSE, FALSE),
                                 "character")
  num_check <- check_fun_params(list(skip = args$skip), 1, FALSE, "numeric")
  logi_check <- check_fun_params(list(rownames = args$rownames),
                                 1, FALSE, "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]

  seperators <- which(raw_lines == sep)
  if (length(seperators) != 0) {
    if (seperators[length(seperators)] != length(raw_lines)) {
      seperators <- c(seperators, length(raw_lines) + 1)
    }
    if (seperators[1] == 1) {
      motif_stops <- seperators[-1] - 1
    } else motif_stops <- seperators - 1
  } else motif_stops <- length(raw_lines)

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
      motif_starts <- headers + 1
      headers <- raw_lines[headers]
    }
  } else {
    motif_starts <- c(1, seperators + 1)
  }

  motifs <- bpmapply(function(x, y) raw_lines[x:y],
                     motif_starts, motif_stops, SIMPLIFY = FALSE)
  if (rownames && positions == "columns") {
    motifs <- bplapply(motifs, function(x) as.matrix(read.table(text = x,
                                                                row.names = 1)))
  } else if (rownames && positions == "rows") {
    motifs <- bplapply(motifs, function(x) as.matrix(read.table(text = x,
                                                                header = TRUE)))
  } else {
    motifs <- bplapply(motifs, function(x) as.matrix(read.table(text = x)))
  }

  if (positions == "rows") {
    motifs <- bplapply(motifs, t)
  }

  motifs <- bplapply(motifs, function(x) {rownames(x) <- NULL; x})
  if (!missing(type)) {
    motifs <- bplapply(motifs, function(x) {
                      universalmotif_cpp(motif = x, type = type, alphabet = alphabet)
                     })
  } else {
    motifs <- bplapply(motifs, function(x) universalmotif_cpp(motif = x,
                                                          alphabet = alphabet))
  }

  if (!isFALSE(headers)) {
    motifs <- bpmapply(function(x, y) {x@name <- y; x}, motifs, headers,
                       SIMPLIFY = FALSE)
  }

  motifs <- bplapply(motifs, function(x) {
                         msg <- validObject_universalmotif(x)
                         if (length(msg) > 0) stop(msg) else x
                       })

  if (length(motifs) == 1) motifs <- motifs[[1]]

  motifs

}
