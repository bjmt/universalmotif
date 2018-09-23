#' Export motifs in HOMER format.
#'
#' Convert motifs to HOMER format and write to file.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#' @param logodds_threshold `numeric` Stringency required for HOMER to match a motif.
#'    See [scan_sequences()].
#'
#' @return NULL, invisibly.
#'
#' @references
#'    \insertRef{homer}{universalmotif}
#'
#' @examples
#' motif <- create_motif()
#' write_homer(motif, tempfile())
#'
#' @family write_motifs
#' @seealso [read_homer()]
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_homer <- function(motifs, file, logodds_threshold = 0.6) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file), 1, FALSE, "character")
  num_check <- check_fun_params(list(logodds_threshold = args$logodds_threshold),
                                1, FALSE, "numeric")
  all_checks <- c(char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, "PWM")
  if (!is.list(motifs)) motifs <- list(motifs)

  max_logodds <- vapply(motifs, function(x) sum(apply(x["motif"], 2, max)),
                        numeric(1))
  logodds_thresholds <- max_logodds * logodds_threshold

  motifs <- convert_type(motifs, "PPM")

  .write_homer <- function(motifs, logodds_thresholds) {
    motif <- motifs
    threshold <- logodds_thresholds
    header <- c(paste0(">", motif["consensus"]), motif["name"], threshold)
    header <- paste(header, collapse = "\t")
    mat <- t(motif["motif"])
    lines_out <- vector("character", 1 + nrow(mat))
    lines_out[1] <- header
    for (i in seq_len(nrow(mat))) {
      pos <- mat[i, ]
      pos <- vapply(pos, function(x) formatC(x, format = "f", digits = 3),
                    character(1))
      lines_out[1 + i] <- paste(pos, collapse = "\t")
    }
    lines_out
  }

  lines_out <- mapply(.write_homer, motifs, logodds_thresholds,
                        SIMPLIFY = FALSE)
  lines_out <- unlist(lines_out)

  writeLines(lines_out, con <- file(file))
  close(con)

  invisible(NULL)

}
