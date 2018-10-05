#' Import TRANSFAC motifs.
#'
#' Import TRANSFAC formatted motifs. Assumed to be DNA motifs, type PCM.
#' See `system.file("extdata", "transfac.txt", pacakge="universalmotif")`
#' for an example motif.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @examples
#' transfac <- read_transfac(system.file("extdata", "transfac.txt",
#'                                       package = "universalmotif"))
#' 
#' @references
#'    \insertRef{transfac}{universalmotif}
#'
#' @family read_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams read_cisbp
#' @export
read_transfac <- function(file, skip = 0) {

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
  raw_lines <- raw_lines[!grepl("^XX", raw_lines)]
  raw_lines <- raw_lines[raw_lines != ""]

  motif_stops <- grep("^//", raw_lines)

  if (1 == motif_stops[1]) {
    motif_stops <- motif_stops[-1]
    raw_lines <- raw_lines[-1]
  }
  motif_starts <- c(1, motif_stops[-length(motif_stops)] + 1)
  motif_stops <- motif_stops - 1

  motifs <- mapply(function(x, y) raw_lines[x:y],
                     motif_starts, motif_stops,
                     SIMPLIFY = FALSE)

  get_matrix <- function(x) {
    mot_start <- which(grepl("^P0", x) | grepl("^PO", x)) + 1
    x <- x[-c(1:(mot_start - 1))]
    mot_i <- vapply(x, function(x) {
                         x <- strsplit(x, "\\s+")[[1]][1]
                         x <- suppressWarnings(as.numeric(x))
                         !is.na(x)
                       }, logical(1))
    mot <- x[mot_i]
    per_line <- function(x) {
      x <- strsplit(x, "\\s+")[[1]]
      as.numeric(x[-c(1, 6)])
    }
    mot <- vapply(mot, per_line, numeric(4))
    matrix(mot, ncol = 4, byrow = TRUE)
  }

  motif_matrix <- lapply(motifs, get_matrix)

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

  motif_meta <- lapply(motifs, parse_meta)

  motifs <- mapply(function(x, y) {
                      mot <- universalmotif_cpp(name = as.character(y[names(y) == 
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
                      msg <- validObject_universalmotif(mot)
                      if (length(msg) > 0) stop(msg) else mot
                     }, motif_matrix, motif_meta)

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs

}
