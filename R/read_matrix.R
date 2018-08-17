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
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
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
                        rownames = FALSE, BPPARAM = SerialParam()) {

  args <- as.list(environment())
  check_input_params(char = list(positions = args$positions, file = args$file,
                                 alphabet = args$alphabet, sep = args$sep),
                     num = list(skip = args$skip),
                     logi = list(rownames = args$rownames))

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
                     motif_starts, motif_stops, SIMPLIFY = FALSE,
                     BPPARAM = BPPARAM)
  if (rownames && positions == "columns") {
    motifs <- bplapply(motifs, function(x) as.matrix(read.table(text = x,
                                                                row.names = 1)),
                      BPPARAM = BPPARAM)
  } else if (rownames && positions == "rows") {
    motifs <- bplapply(motifs, function(x) as.matrix(read.table(text = x,
                                                                header = TRUE)),
                      BPPARAM = BPPARAM)
  } else {
    motifs <- bplapply(motifs, function(x) as.matrix(read.table(text = x)),
                      BPPARAM = BPPARAM)
  }

  if (positions == "rows") {
    motifs <- bplapply(motifs, t, BPPARAM = BPPARAM)
  }

  motifs <- bplapply(motifs, function(x) {rownames(x) <- NULL; x}, BPPARAM = BPPARAM)
  if (!missing(type)) {
    motifs <- bplapply(motifs, function(x) {
                      universalmotif_cpp(motif = x, type = type, alphabet = alphabet)
                     }, BPPARAM = BPPARAM)
  } else {
    motifs <- bplapply(motifs, function(x) universalmotif_cpp(motif = x,
                                                          alphabet = alphabet),
                       BPPARAM = BPPARAM)
  }

  if (!isFALSE(headers)) {
    motifs <- bpmapply(function(x, y) {x@name <- y; x}, motifs, headers,
                       BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  }

  motifs <- bplapply(motifs, function(x) {
                         msg <- validObject_universalmotif(x)
                         if (length(msg) > 0) stop(msg) else x
                       }, BPPARAM = BPPARAM)

  if (length(motifs) == 1) motifs <- motifs[[1]]

  motifs

}
