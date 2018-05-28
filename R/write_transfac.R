#' Export motifs in TRANSFAC format.
#'
#' @param motifs List of motifs or a motif object.
#' @param file Character.
#' @param BPPARAM Param for bplapply.
#'
#' @return NULL, invisibly.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                     package = "universalmotif"))
#' write_transfac(jaspar, tempfile())
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_transfac <- function(motifs, file, BPPARAM = bpparam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, "PCM", BPPARAM = BPPARAM)
  if (!is.list(motifs)) motifs <- list(motifs)

  .write_transfac <- function(motifs) {
    motif <- motifs
    lines_out <- vector()
    lines_out <- c(lines_out, paste("ID", motif["name"]))
    if (length(motif["altname"]) > 0) {
      lines_out <- c(lines_out, paste("NA", motif["altname"]))
    }
    if (length(motif["family"]) > 0) {
      lines_out <- c(lines_out, paste("HC", motif["family"]))
    }
    if (length(motif["organism"]) > 0) {
      lines_out <- c(lines_out, paste("OS", motif["organism"]))
    }
    lines_out <- c(lines_out, "P0\tA\tC\tG\tT")
    consensus <- strsplit(motif["consensus"], "")[[1]]
    for (j in seq_along(consensus)) {
      p1 <- formatC(j, width = 2, format = "d", flag = "0")
      p2 <- paste(as.numeric(motif["motif"][, j]), collapse = "\t")
      p3 <- consensus[j]
      lines_out <- c(lines_out, paste(p1, p2, p3, sep = "\t"))
    }
    lines_out <- c(lines_out, "XX", "//")
  }

  lines_out <- bplapply(motifs, .write_transfac, BPPARAM = BPPARAM)
  lines_out <- unlist(lines_out)

  writeLines(lines_out, con <- file(file))
  close(con)

  invisible(NULL)

}
