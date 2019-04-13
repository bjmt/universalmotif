#' Export motifs in TRANSFAC format.
#'
#' Convert motifs to TRANSFAC format and write to file.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#' @param overwrite `logical(1)` Overwrite existing file.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                     package = "universalmotif"))
#' write_transfac(jaspar, tempfile())
#'
#' @references
#'    \insertRef{transfac}{universalmotif}
#'
#' @family write_motifs
#' @seealso [read_transfac()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_transfac <- function(motifs, file, overwrite = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file), 1, FALSE, "character")
  logi_check <- check_fun_params(list(overwrite = args$overwrite),
                                 1, FALSE, "logical")
  all_checks <- c(char_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (file.exists(file) && !overwrite)
    stop(wmsg("Existing file found, set `overwrite = TRUE` to continue."))

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PCM")
  if (!is.list(motifs)) motifs <- list(motifs)

  .write_transfac <- function(motifs) {
    motif <- motifs
    lines_out <- vector()
    lines_out <- c(lines_out, paste("ID", motif@name))
    if (length(motif@altname) > 0) {
      lines_out <- c(lines_out, paste("NA", motif@altname))
    }
    if (length(motif@family) > 0) {
      lines_out <- c(lines_out, paste("HC", motif@family))
    }
    if (length(motif@organism) > 0) {
      lines_out <- c(lines_out, paste("OS", motif@organism))
    }
    lines_out <- c(lines_out, "P0\tA\tC\tG\tT")
    consensus <- safeExplode(motif@consensus)
    for (j in seq_along(consensus)) {
      p1 <- formatC(j, width = 2, format = "d", flag = "0")
      p2 <- paste0(as.numeric(motif@motif[, j]), collapse = "\t")
      p3 <- consensus[j]
      lines_out <- c(lines_out, paste(p1, p2, p3, sep = "\t"))
    }
    lines_out <- c(lines_out, "XX", "//")
  }

  lines_out <- lapply(motifs, .write_transfac)
  lines_out <- unlist(lines_out)

  writeLines(lines_out, con <- file(file))
  close(con)

  invisible(NULL)

}
