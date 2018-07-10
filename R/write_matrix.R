#' Export motifs as raw matrices.
#'
#' @param motifs List of motifs or a single motif object.
#' @param file Character.
#' @param positions Character. One of 'columns' or 'rows'.
#' @param rownames Logical. Include alphabet letters as rownames.
#' @param type Character. One of 'PCM', 'PPM', 'PWM', and 'ICM'. If missing
#'             will use whatever type the motif is currently stored as.
#' @param sep Character. Indicates how to separate individual motifs.
#' @param headers Logical or character, indicating if and how to write names.
#' @param BPPARAM See \code{\link[BiocParallel]{SerialParam}}.
#'
#' @return NULL, invisibly.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_matrix <- function(motifs, file, positions = "columns", rownames = FALSE,
                         type, sep = "", headers = TRUE, BPPARAM = SerialParam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  if (!missing(type)) motifs <- convert_type(motifs, type, BPPARAM = BPPARAM)
  if (!is.list(motifs)) motifs <- list(motifs)

  .write_matrix <- function(motifs, positions, rownames, sep, headers) {

    motif <- motifs

    lines_out <- vector()

    if (!isFALSE(headers)) {
      if (isTRUE(headers)) lines_out <- c(lines_out, motif["name"]) else {
        lines_out <- c(lines_out, paste0(headers, motif["name"]))
      }
    }

    if (positions == "columns") {
      for (i in seq_len(nrow(motif["motif"]))) {
        pos <- motif["motif"][i, ]
        # pos <- vapply(pos, function(x) formatC(x, format = "f", digits = 3),
                      # character(1))
        pos <- as.character(pos)
        if (rownames) pos <- c(rownames(motif["motif"])[i], pos)
        lines_out <- c(lines_out, paste(pos, collapse = "\t"))
      }
    } else if (positions == "rows") {
      if (rownames) lines_out <- c(lines_out, paste(rownames(motif["motif"]),
                                                    collapse = "\t"))
      for (i in seq_len(ncol(motif["motif"]))) {
        pos <- motif["motif"][, i]
        # pos <- vapply(pos, function(x) formatC(x, format = "f", digits = 3),
                      # character(1))
        pos <- as.character(pos)
        lines_out <- c(lines_out, paste(pos, collapse = "\t"))
      }
    }

    c(lines_out, sep)

  }

  lines_final <- bplapply(motifs, .write_matrix, positions = positions,
                          rownames = rownames, sep = sep,
                          headers = headers, BPPARAM = BPPARAM)  # not working??
  lines_final <- unlist(lines_final)

  writeLines(lines_final, con <- file(file))
  close(con)

  invisible(NULL)

}
