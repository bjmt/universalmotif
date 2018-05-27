#' Import HOMER motifs.
#'
#' Assumed DNA.
#'
#' @param file Character.
#' @param skip Numeric.
#'
#' @return List of universalmotif objects.
#'
#' @examples
#' homer <- read_homer(system.file("extdata", "homer.txt",
#'                                 package = "universalmotif"))
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_homer <- function(file, skip = 0) {

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
    universalmotif(name = x[1],
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
  }

  motifs <- mapply(homer2umot, motif_meta, motif_list,
                   SIMPLIFY = FALSE)

  motifs

}
