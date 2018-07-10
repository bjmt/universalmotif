#' Export motifs in JASPAR format.
#'
#' @param motifs List of motifs or a motif object.
#' @param file Character.
#' @param BPPARAM See \code{\link[BiocParallel]{SerialParam}}.
#'
#' @return NULL, invisibly.
#'
#' @examples
#' transfac <- read_transfac(system.file("extdata", "transfac.txt",
#'                                     package = "universalmotif"))
#' write_jaspar(transfac, tempfile())
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_jaspar <- function(motifs, file, BPPARAM = SerialParam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, "PCM", BPPARAM = BPPARAM)
  if (!is.list(motifs)) motifs <- list(motifs)

  .write_jaspar <- function(motifs) {

    lines_out <- vector()

    motif <- motifs
    if (length(motif["altname"]) > 0) {
      lines_out <- c(lines_out, paste0(">", motif["name"], " ",
                                       motif["altname"]))
    } else {
      lines_out <- c(lines_out, paste0(">", motif["name"]))
    }

    if (motif["alphabet"] == "DNA") {
      alph <- c("A", "C", "G", "T") 
    } else if (motif["alphabet"] == "RNA") {
      alph <- c("A", "C", "G", "U")
    } else if (motif["alphabet"] == "AA") {
      alph <- c("A", "C", "D", "E", "F", "G", "H", "I", "K",
              "L", "M", "N", "P", "Q", "R", "S", "T", "V",
              "W", "Y")
    } else stop("unknown alphabet")

    nsites <- motif["nsites"]
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
      p3 <- as.numeric(motif["motif"][j, ])
      p3 <- vapply(p3, function(x) formatC(x, width = width, format = "d"),
                   character(1))
      p3 <- paste(p3, collapse = " ")
      p4 <- "]"
      lines_out <- c(lines_out, paste(p1, p2, p3, p4))
    }

    lines_out <- c(lines_out, "")

  }

  lines_out <- bplapply(motifs, .write_jaspar, BPPARAM = BPPARAM)
  lines_out <- unlist(lines_out)

  writeLines(lines_out, con <- file(file))
  close(con)

  invisible(NULL)

}
