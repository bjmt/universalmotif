#' Export motifs in HOMER format.
#'
#' Convert DNA motifs to HOMER format and write to file.
#' See \url{http://homer.ucsd.edu/homer/motif/}.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#' @param logodds_threshold Deprecated. If set, [read_homer()] will behave like
#'    pre-version 1.12.0 of the `universalmotif` package for backwards
#'    compatibility (though a warning will be printed).
#' @param overwrite `logical(1)` Overwrite existing file.
#' @param append `logical(1)` Add to an existing file.
#' @param threshold `numeric(1)` Stringency required for HOMER to match a motif.
#'    See [scan_sequences()] for how to use this argument. Can be a single value to
#'    be recycled for all motifs, or a vector of equal length to the number of motifs.
#' @param threshold.type `character(1)` How the `threshold` value should be used
#'    to obtain the final threshold value in the written motif. See [scan_sequences()]
#'    for how to use this.
#'
#' @return `NULL`, invisibly.
#'
#' @references
#'
#' Heinz S, Benner C, Spann N, Bertolino E, Lin YC, Laslo P, Cheng
#' JX, Murre C, Singh H, Glass CK (2010). “Simple combinations of
#' lineage-determining transcription factors prime cis-regulatory
#' elements required for macrophage and B cell identities.”
#' *Molecular Cell*, **38**, 576-589.
#'
#' @examples
#' motif <- create_motif()
#' write_homer(motif, tempfile())
#'
#' @family write_motifs
#' @seealso [read_homer()]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
write_homer <- function(motifs, file, logodds_threshold = NULL,
  overwrite = FALSE, append = FALSE, threshold = 0.8,
  threshold.type = c("logodds", "logodds.abs", "pvalue")) {

  threshold.type <- match.arg(threshold.type)

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file), 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(threshold = args$threshold,
                                     logodds_threshold = args$logodds_threshold),
                                c(0, 1), c(FALSE, TRUE), TYPE_NUM)
  logi_check <- check_fun_params(list(overwrite = args$overwrite,
                                      append = args$append),
                                 c(1, 1), c(FALSE, FALSE), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (file.exists(file) && !overwrite && !append)
    stop(wmsg("Existing file found, set `overwrite = TRUE` to continue."))

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  motifsPWM <- convert_type_internal(motifs, "PWM")

  if (!is.null(logodds_threshold)) {
    warning(wmsg("Deprecated `logodds_threshold` argument has been ",
        "set. In future please use `threshold` and `threshold.type`."),
      call. = FALSE)
    max_logodds <- vapply(motifsPWM, function(x) sum(apply(x@motif, 2, max)),
                          numeric(1))
    logodds_thresholds <- max_logodds * logodds_threshold
  } else {
    if (length(threshold) != 1 && length(threshold) != length(motifs))
      stop("A single value or one per motif should be provided for `threshold`")
    if (threshold.type == "logodds") {
      max_logodds <- vapply(motifsPWM, function(x) sum(apply(x@motif, 2, max)),
                            numeric(1))
      logodds_thresholds <- max_logodds * threshold
    } else if (threshold.type == "logodds.abs") {
      logodds_thresholds <- threshold
    } else if (threshold.type == "pvalue") {
      logodds_thresholds <- motif_pvalue(motifs, pvalue = threshold)
    }
  }

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
