#' Import CIS-BP motifs.
#'
#' @param file Character.
#' @param skip Numeric.
#'
#' @return List of universalmotif objects.
#'
#' @examples
#' cisbp <- read_cisbp(system.file("extdata", "cisbp.txt",
#'                                 package = "universalmotif"))
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_cisbp <- function(file, skip = 0) {

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]

  meta_starts <- which(grepl("^TF", raw_lines))
  for (i in seq_along(meta_starts)) {
    if ((meta_starts[i] + 1) %in% meta_starts) {
      meta_starts <- meta_starts[-(i + 1)]
    }
  }

  if (length(meta_starts) == 1) motif_stops <- length(raw_lines) else {
    motif_stops <- c(meta_starts[-1] - 1, length(raw_lines))
  }
  motif_starts <- which(grepl("^Pos", raw_lines))
  meta_stops <- motif_starts - 1

  meta_list <- bpmapply(function(x, y) raw_lines[x:y],
                        meta_starts, meta_stops,
                        SIMPLIFY = FALSE)
  motif_list <- bpmapply(function(x, y) raw_lines[x:y],
                         motif_starts, motif_stops,
                         SIMPLIFY = FALSE)

  parse_meta <- function(x) {
    metas <- lapply(x, function(x) strsplit(x, "\\s+")[[1]])
    metas_correct <- c(name = metas[[1]][2],
                       altname = metas[[2]][3],
                       family = metas[[5]][2],
                       organism = metas[[6]][2])
    metas_correct
  }

  parse_motifs <- function(x) {
    alph <- strsplit(x[1], "\\s+")[[1]][-1]
    x <- vapply(x[-1], function(x) as.numeric(strsplit(x, "\\s+")[[1]][-1]),
                numeric(length(alph)))
    x <- matrix(x, ncol = 4, byrow = TRUE)
    colnames(x) <- alph
    x
  }

  meta_list <- bplapply(meta_list, parse_meta)
  motif_list <- bplapply(motif_list, parse_motifs)

  motifs <- bpmapply(function(x, y) {
                    if (all(colnames(x) %in% c("A", "C", "G", "U"))) {
                      alph <- "RNA"
                    } else if (all(colnames(x) %in% c("A", "C", "G", "T"))) {
                      alph <- "DNA"
                    }
                    universalmotif(name = y[1],
                                   altname = y[2],
                                   family = y[3],
                                   organism = y[4],
                                   motif = t(x),
                                   alphabet = alph,
                                   type = "PPM")
                  }, motif_list, meta_list)

  motifs <- motifs[vapply(motifs, function(x) ncol(x["motif"]) > 0,
                          logical(1))]

  motifs

}
