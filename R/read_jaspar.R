#' Import JASPAR motifs.
#'
#' Import JASPAR formatted motifs. See \url{http://jaspar.genereg.net/}.
#' Can be either DNA, RNA, or AA motifs.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#'
#' @references
#'    \insertRef{jaspar}{universalmotif}
#'
#' @family read_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams read_cisbp
#' @export
read_jaspar <- function(file, skip = 0) {

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

  motifs <- mapply(function(x, y) raw_lines[x:y],
                     motif_starts, motif_stops,
                     SIMPLIFY = FALSE)

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

  motifs <- lapply(motifs, get_matrix)

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
    mot <- universalmotif_cpp(name = name[1], altname = name[2],
                   type = "PCM", alphabet = alphabet,
                   motif = motif)
    msg <- validObject_universalmotif(mot)
    if (length(msg) > 0) stop(msg) else mot
  }

  motifs <- mapply(jaspar2umot, motifs, motif_names, 
                     SIMPLIFY = FALSE)

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs

}
