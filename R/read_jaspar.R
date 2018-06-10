#' Import JASPAR motifs.
#'
#' @param file Character.
#' @param skip Numeric.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return List of universalmotif objects.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_jaspar <- function(file, skip = 0, BPPARAM = bpparam()) {

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]

  motif_names <- which(grepl("^>", raw_lines))
  motif_starts <- motif_names + 1
  if (length(motif_starts) == 0) motif_stops <- length(raw_lines) else {
    motif_stops <- c(motif_names[-1] - 1, length(raw_lines))
  }

  if (length(unique(c(length(motif_names), length(motif_starts),
                      length(motif_stops)))) != 1) {
    stop("motifs incorrectly formatted")
  }

  motif_names <- raw_lines[motif_names]
  motif_names <- sub(">", "", motif_names)
  motif_names <- lapply(motif_names, function(x) strsplit(x, "\\s+")[[1]])

  motifs <- bpmapply(function(x, y) raw_lines[x:y],
                     motif_starts, motif_stops,
                     SIMPLIFY = FALSE, BPPARAM = BPPARAM)

  get_matrix <- function(x) {
    x <- sub("\\[", "", x)
    x <- sub("\\]", "", x)
    per_line1 <- function(x) {
      x <- strsplit(x, "\\s+")[[1]]
      x <- x[x != ""]
      as.numeric(x[-1])
    }
    per_line2 <- function(x) {
      x <- strsplit(x, "\\s+")[[1]]
      x <- x[x != ""]
      as.numeric(x)
    }
    alphabet <- vapply(x, function(x) strsplit(x, "\\s+")[[1]][1],
                       character(1))
    if (any(alphabet %in% LETTERS)) {
      x2 <- sapply(x, per_line1)
      x2 <- matrix(x2, nrow = length(x), byrow = TRUE)
      rownames(x2) <- alphabet
      x2
    } else {
      x2 <- sapply(x, per_line2)
      matrix(x2, nrow = length(x), byrow = TRUE)
    }
  }

  motifs <- bplapply(motifs, get_matrix, BPPARAM = BPPARAM)

  jaspar2umot <- function(motif, name) {
    alphabet <- rownames(motif)
    if (all(c("A", "C", "D", "E", "F", "G", "H", "I", "K",
              "L", "M", "N", "P", "Q", "R", "S", "T", "V",
              "W", "Y") %in% alphabet)) {
      alphabet <- "AA" 
    } else if (all(c("A", "C", "G", "U") %in% alphabet)) {
      alphabet <- "RNA" 
    } else if (all(c("A", "C", "G", "T") %in% alphabet)) {
      alphabet <- "DNA"
    } else alphabet <- "DNA"
    universalmotif(name = name[1], altname = name[2],
                   type = "PCM", alphabet = alphabet,
                   motif = motif)
  }

  motifs <- bpmapply(jaspar2umot, motifs, motif_names, 
                     SIMPLIFY = FALSE, BPPARAM = BPPARAM)

  motifs

}
