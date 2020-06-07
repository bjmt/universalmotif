#' Export motifs in HOMER format.
#'
#' Convert DNA motifs to HOMER format and write to file.
#' See \url{http://homer.ucsd.edu/homer/motif/}.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#' @param logodds_threshold `numeric` Stringency required for HOMER to match a motif.
#'    See [scan_sequences()].
#' @param overwrite `logical(1)` Overwrite existing file.
#' @param append `logical(1)` Add to an existing file.
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
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_homer <- function(motifs, file, logodds_threshold = 0.6,
                        overwrite = FALSE, append = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file), 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(logodds_threshold = args$logodds_threshold),
                                1, FALSE, TYPE_NUM)
  logi_check <- check_fun_params(list(overwrite = args$overwrite,
                                      append = args$append),
                                 c(1, 1), c(FALSE, FALSE), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (file.exists(file) && !overwrite && !append)
    stop(wmsg("Existing file found, set `overwrite = TRUE` to continue."))

  motifs <- convert_motifs(motifs)
  motifs <- convert_type_internal(motifs, "PWM")
  if (!is.list(motifs)) motifs <- list(motifs)

  max_logodds <- vapply(motifs, function(x) sum(apply(x@motif, 2, max)),
                        numeric(1))
  logodds_thresholds <- max_logodds * logodds_threshold

  motifs <- convert_type(motifs, "PPM")

  .write_homer <- function(motifs, logodds_thresholds) {
    motif <- motifs
    threshold <- logodds_thresholds
    header <- c(paste0(">", motif["consensus"]), motif@name, threshold)
    header <- paste0(header, collapse = "\t")
    mat <- t(motif@motif)
    lines_out <- vector("character", 1 + nrow(mat))
    lines_out[1] <- header
    for (i in seq_len(nrow(mat))) {
      pos <- mat[i, ]
      pos <- formatC(pos, format = "f", digits = 3)
      lines_out[1 + i] <- paste0(pos, collapse = "\t")
    }
    lines_out
  }

  lines_out <- mapply(.write_homer, motifs, logodds_thresholds,
                        SIMPLIFY = FALSE)
  lines_out <- unlist(lines_out)

  if (append) {
    cat(lines_out, sep = "\n", file = file, append = TRUE)
  } else {
    writeLines(lines_out, con <- file(file)); close(con)
  }

  invisible(NULL)

}
