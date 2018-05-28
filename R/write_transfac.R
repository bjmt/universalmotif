#' Export motifs in TRANSFAC format.
#'
#' @param motifs List of motifs or a motif object.
#' @param file Character.
#'
#' @return NULL, invisibly.
#'
#' @examples
#' \dontrun{
#' jaspar <- read_transfac(system.file("extdata", "jaspar.txt",
#'                                     package = "universalmotif"))
#' write_transfac(jaspar, "motifs.transfac")
#' }
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_transfac <- function(motifs, file) {

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, "PCM")
  if (!is.list(motifs)) motifs <- list(motifs)

  lines_out <- vector()
  for (i in seq_along(motifs)) {
    motif <- motifs[[i]]
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

  writeLines(lines_out, con <- file(file))
  close(con)

  invisible(NULL)

}
