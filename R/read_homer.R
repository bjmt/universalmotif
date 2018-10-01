#' Import HOMER motifs.
#'
#' Import HOMER formatted motifs. See \url{http://homer.ucsd.edu/homer/motif/}.
#' Assumed to be DNA motifs.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @examples
#' homer <- read_homer(system.file("extdata", "homer.txt",
#'                                 package = "universalmotif"))
#'
#' @references
#'    \insertRef{homer}{universalmotif}
#'
#' @family read_motifs
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams read_cisbp
#' @export
read_homer <- function(file, skip = 0) {

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

  headers <- which(grepl("^>", raw_lines))
  motif_starts <- headers + 1
  if (length(headers) == 1) motif_stops <- length(raw_lines) else {
    motif_stops <- c(headers[-1] - 1, length(raw_lines))
  }

  parse_motifs <- function(x, y) {
    motif <- raw_lines[x:y]
    motif <- vapply(motif,
                    function(x) as.numeric(strsplit(x, "\\s+")[[1]]),
                    numeric(4))
    matrix(motif, ncol = 4, byrow = TRUE)
  }

  parse_meta <- function(x) {
    x <- strsplit(x, "\\s+")[[1]]
    if (grepl("/", x[2])) {
      y <- strsplit(x[2], "/")[[1]]
      y <- strsplit(y[1], "\\(")[[1]]
      x[2] <- y[1]
      family <- strsplit(y[2], "\\)")[[1]][1]
    } else family <- character(0)
    x2 <- strsplit(x[6], ",")[[1]]
    nsites <- strsplit(strsplit(x2[1], "T:")[[1]][2], "\\(")[[1]][1]
    bkgsites <- strsplit(strsplit(x2[2], "B:")[[1]][2], "\\(")[[1]][1]
    pval <- strsplit(x2[3], "P:")[[1]][2]
    c(name = x[2], nsites = nsites, bkgsites = bkgsites, pval = pval,
      threshold = x[3], family = family)
  }

  motif_list <- mapply(parse_motifs, motif_starts, motif_stops,
                         SIMPLIFY = FALSE)

  motif_meta <- lapply(raw_lines[headers], parse_meta)

  homer2umot <- function(x, y) {
    mot <- universalmotif_cpp(name = x[1],
                   nsites = ifelse(is.na(as.numeric(x[2])), numeric(0),
                                   as.numeric(x[2])),
                   bkgsites = ifelse(is.na(as.numeric(x[3])), numeric(0),
                                     as.numeric(x[3])),
                   pval = ifelse(is.na(as.numeric(x[4])), numeric(0),
                                 as.numeric(x[4])),
                   motif = t(y),
                   alphabet = "DNA",
                   type = "PPM",
                   family = x[6],
                   extrainfo = c(logodds = x[5]))
    msg <- validObject_universalmotif(mot)
    if (length(msg) > 0) stop(msg) else mot
  }

  motifs <- mapply(homer2umot, motif_meta, motif_list,
                     SIMPLIFY = FALSE)

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs

}
