#' Export motifs in JASPAR format.
#'
#' Convert motifs to JASPAR format and write to file.
#' See \url{http://jaspar.genereg.net/}.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' transfac <- read_transfac(system.file("extdata", "transfac.txt",
#'                                     package = "universalmotif"))
#' write_jaspar(transfac, tempfile())
#'
#' @references
#'    \insertRef{jaspar}{universalmotif}
#'
#' @family write_motifs
#' @seealso [read_jaspar()]
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_jaspar <- function(motifs, file) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file),
                                 1, FALSE, "character")
  all_checks <- c(char_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, "PCM")
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
    if (length(nsites) == 0) nsites <- 100
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

  lines_out <- lapply(motifs, .write_jaspar)
  lines_out <- unlist(lines_out)

  writeLines(lines_out, con <- file(file))
  close(con)

  invisible(NULL)

}
