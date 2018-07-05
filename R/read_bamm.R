#' Import BaMM motifs.
#'
#' @param file Character.
#' @param alphabet Character.
#' @param name Character.
#' @param skip Numeric.
#'
#' @return A universalmotif object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_bamm <- function(file, alphabet = "DNA", name = "motif", skip = 0) {

  raw_lines <- readLines(con <- file(file)); close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]
  numbers <- lapply(raw_lines, function(x) scan(text = x, quiet = TRUE))
  line_lengths <- vapply(numbers, length, numeric(1))
  max_order <- length(unique(line_lengths))
  zero_order <- matrix(nrow = line_lengths[1],
                       ncol = length(numbers) / max_order)
  first_order <- NULL
  if (max_order >= 2) {
    first_order <- matrix(nrow = nrow(zero_order) * 4, ncol = ncol(zero_order))
  }
  second_order <- NULL
  if (max_order >= 3) {
    second_order <- matrix(nrow = nrow(first_order) * 4, ncol = ncol(zero_order))
  }
  
  j <- 1
  for (i in seq_len(ncol(zero_order))) {
    if (j > length(numbers)) break
    zero_order[, i] <- numbers[[j]]
    j <- j + max_order
  }

  if (!is.null(first_order)) {
    j <- 2
    for (i in seq_len(ncol(first_order))) {
      if (j > length(numbers)) break
      first_order[, i] <- numbers[[j]] / sum(numbers[[j]])
      j <- j + max_order
    }
  } else first_order <- matrix()

  if (!is.null(second_order)) {
    j <- 3
    for (i in seq_len(ncol(second_order))) {
      if (j > length(numbers)) break
      second_order[, i] <- numbers[[j]] / sum(numbers[[j]])
      j <- j + max_order
    }
  } else second_order <- matrix()

  motif <- universalmotif(alphabet = alphabet, name = name, motif = zero_order,
                          difreq = first_order, trifreq = second_order)
  motif

}
