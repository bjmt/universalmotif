#' Import TRANSFAC motifs.
#'
#' Assumed to be DNA motifs, type PCM.
#'
#' @param file Character.
#' @param skip Numeric.
#'
#' @return List of universalmotif objects.
#'
#' @examples
#' transfac <- read_transfac(system.file("extdata", "transfac.txt",
#'                                       package = "universalmotif"))
#' 
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_transfac <- function(file, skip = 0) {

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[!grepl("^XX", raw_lines)]

  motif_starts <- which(grepl("^P0", raw_lines) | grepl("^PO", raw_lines))
  motif_stops <- which(grepl("^//", raw_lines))

  for (i in seq_along(motif_stops)) {
    if ((motif_stops[i] + 1) %in% motif_stops) {
      motif_stops <- motif_stops[-(i + 1)]
    }
  }

  if (length(motif_starts) != length(motif_stops)) {
    stop("motifs incorrectly formatted")
  }

  motifs <- mapply(function(x, y) raw_lines[x:(y - 1)],
                   motif_starts, motif_stops,
                   SIMPLIFY = FALSE)

  motif_stops <- c(1, motif_stops[-length(motif_stops)])

  motif_meta <- mapply(function(x, y) raw_lines[(x):(y - 1)],
                       motif_stops, motif_starts,
                       SIMPLIFY = FALSE)

  get_matrix <- function(x) {
    x <- x[-1]
    per_line <- function(x) {
      x <- strsplit(x, "\\s+")[[1]]
      as.numeric(x[-c(1, 6)])
    }
    x <- vapply(x, per_line, numeric(4))
    matrix(x, ncol = 4, byrow = TRUE)
  }

  motifs <- lapply(motifs, get_matrix)

  parse_meta <- function(x) {
    metas <- lapply(x, function(x) strsplit(x, "\\s+")[[1]])
    metas_correct <- vector()
    for (i in seq_along(metas)) {
      if (length(metas[[i]]) == 0) next
      if (metas[[i]][1] == "AC") {
        metas_correct <- c(metas_correct, AC = metas[[i]][2])
      }
      if (metas[[i]][1] == "ID") {
        metas_correct <- c(metas_correct, ID = metas[[i]][2])
      } 
      if (metas[[i]][1] == "NA") {
        metas_correct <- c(metas_correct, N.A = metas[[i]][2])
      }
      if (metas[[i]][1] == "HC") {
        metas_correct <- c(metas_correct, family = metas[[i]][2])
      }
      if (metas[[i]][1] == "OS") {
        metas_correct <- c(metas_correct, organism = metas[[i]][2])
      }
      if (all(c("AC", "ID") %in% names(metas_correct))) {
        metas_correct <- c(metas_correct,
                           name = metas_correct[names(metas_correct) ==
                                                "AC"],
                           altname = metas_correct[names(metas_correct) == 
                                                   "ID"])
      } else if (all(c("AC", "N.A") %in% names(metas_correct))) {
        metas_correct <- c(metas_correct,
                           name = metas_correct[names(metas_correct) == 
                                                "AC"],
                           altname = metas_correct[names(metas_correct) == 
                                                   "N.A"])
      } else if (all(c("ID", "N.A") %in% names(metas_correct))) {
        metas_correct <- c(metas_correct,
                           name = metas_correct[names(metas_correct) == 
                                                "ID"],
                           altname = metas_correct[names(metas_correct) == 
                                                   "N.A"])
      } else {
        metas_correct <- c(metas_correct,
                           name = metas_correct[names(metas_correct) %in%
                                                c("ID", "AC", "N.A")])
      }
    }
    names(metas_correct) <- vapply(names(metas_correct),
                                   function(x) strsplit(x, "[.]")[[1]][1],
                                   character(1))
    metas_correct <- metas_correct[!duplicated(names(metas_correct))]
    metas_correct
  }

  motif_meta <- lapply(motif_meta, parse_meta)

  motifs <- mapply(function(x, y) {
                    universalmotif(name = as.character(y[names(y) == 
                                                       "name"]),
                                   altname = as.character(y[names(y) == 
                                                          "altname"]),
                                   family = as.character(y[names(y) == 
                                                         "family"]),
                                   organism = as.character(y[names(y) == 
                                                           "organism"]),
                                   motif = t(x),
                                   alphabet = "DNA",
                                   type = "PCM")
                   }, motifs, motif_meta)

  motifs

}
