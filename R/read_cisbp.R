#' Import CIS-BP motifs.
#' 
#' Import CIS-BP formatted motifs. See \url{http://cisbp.ccbr.utoronto.ca/index.php}.
#' Assumed to be DNA motifs.
#'
#' @param file `character(1)` File name.
#' @param skip `numeric(1)` If not zero, will skip however many desired lines in the
#'    file before starting to read.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @references
#'    \insertRef{cisbp}{universalmotif}
#'
#' @family read_motifs
#' @examples
#' cisbp <- read_cisbp(system.file("extdata", "cisbp.txt",
#'                                 package = "universalmotif"))
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_cisbp <- function(file, skip = 0) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file),
                                 1, FALSE, "character")
  num_check <- check_fun_params(list(skip = args$skip), 1, FALSE, "numeric")
  all_checks <- c(char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

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

  meta_list <- mapply(function(x, y) raw_lines[x:y],
                        meta_starts, meta_stops,
                        SIMPLIFY = FALSE)
  motif_list <- mapply(function(x, y) raw_lines[x:y],
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

  meta_list <- lapply(meta_list, parse_meta)
  motif_list <- lapply(motif_list, parse_motifs)

  motifs <- mapply(function(x, y) {
                    if (all(colnames(x) %in% c("A", "C", "G", "U"))) {
                      alph <- "RNA"
                    } else if (all(colnames(x) %in% c("A", "C", "G", "T"))) {
                      alph <- "DNA"
                    }
                    mot <- universalmotif_cpp(name = y[1],
                                   altname = y[2],
                                   family = y[3],
                                   organism = y[4],
                                   motif = t(x),
                                   alphabet = alph,
                                   type = "PPM")
                    msg <- validObject_universalmotif(mot)
                    if (length(msg) > 0) stop(msg) else mot
                  }, motif_list, meta_list)

  motifs <- motifs[vapply(motifs, function(x) ncol(x["motif"]) > 0,
                          logical(1))]

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs

}
