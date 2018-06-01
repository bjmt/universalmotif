#' Export motifs in HOMER format.
#'
#' @param motifs List of motifs or a motif object.
#' @param file Character.
#' @param logodds_threshold Stringency required for HOMER to match a motif.
#' @param BPPARAM Param for bplapply.
#'
#' @return NULL, invisibly.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_homer <- function(motifs, file, logodds_threshold = 0.6,
                        BPPARAM = bpparam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, "PWM", BPPARAM = BPPARAM)
  if (!is.list(motifs)) motifs <- list(motifs)

  get_max_logodds <- function(x) {
    x <- x["motif"]
    pos_max <- apply(x, 2, function(x) sort(x, decreasing = TRUE)[1])
    sum(pos_max)
  }
  max_logodds <- vapply(motifs, get_max_logodds, numeric(1))
  logodds_thresholds <- max_logodds * logodds_threshold

  motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)

  .write_homer <- function(motifs, logodds_thresholds) {
    motif <- motifs
    threshold <- logodds_thresholds
    lines_out <- vector()
    header <- c(paste0(">", motif["consensus"]), motif["name"], threshold)
    header <- paste(header, collapse = "\t")
    lines_out <- header
    mat <- t(motif["motif"])
    for (i in seq_len(nrow(mat))) {
      pos <- mat[i, ]
      pos <- vapply(pos, function(x) formatC(x, format = "f", digits = 3),
                    character(1))
      lines_out <- c(lines_out, paste(pos, collapse = "\t"))
    }
    lines_out
  }

  lines_out <- bpmapply(.write_homer, motifs, logodds_thresholds,
                        BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  lines_out <- unlist(lines_out)

  writeLines(lines_out, con <- file(file))
  close(con)

  invisible(NULL)

}
