#' Import UNIPROBE motifs.
#'
#' Import UNIPROBE formatted motifs. Assumed DNA.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @examples
#' uniprobe.minimal <- read_uniprobe(system.file("extdata", "uniprobe_minimal.txt",
#'                                               package = "universalmotif"))
#' uniprobe.full <- read_uniprobe(system.file("extdata", "uniprobe_full.txt",
#'                                            package = "universalmotif"))
#'
#' @references
#'    \insertRef{uniprobe}{universalmotif}
#'
#' @family read_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams read_cisbp
#' @export
read_uniprobe <- function(file, skip = 0) {

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

  if (any(grepl("#\\s+Motif:", raw_lines))) {

    # full
    motif_meta <- grep("#\\s+Motif:", raw_lines) + 1
    motif_starts <- grep("^Probability matrix", raw_lines) + 1
    motif_stops <- motif_starts + 3

    motif_meta <- lapply(raw_lines[motif_meta],
                         function(x) strsplit(x, "\\s+")[[1]][-1])
    parse_motifs <- function(mstart, mstop) {
      motif_A <- raw_lines[mstart]
      motif_C <- raw_lines[mstart + 1]
      motif_G <- raw_lines[mstart + 2]
      motif_T <- raw_lines[mstop]
      motif_A <- sapply(motif_A, function(x) strsplit(x, "\\s+")[[1]][-1])
      motif_C <- sapply(motif_C, function(x) strsplit(x, "\\s+")[[1]][-1])
      motif_G <- sapply(motif_G, function(x) strsplit(x, "\\s+")[[1]][-1])
      motif_T <- sapply(motif_T, function(x) strsplit(x, "\\s+")[[1]][-1])
      as.numeric(c(motif_A, motif_C, motif_G, motif_T))
    }

    motif_list <- mapply(parse_motifs, motif_starts, motif_stops,
                           SIMPLIFY = FALSE)
    motif_list <- mapply(function(meta, motif) {
                            universalmotif_cpp(name = meta[1],
                                           motif = matrix(motif, nrow = 4,
                                                          byrow = TRUE),
                                           alphabet = "DNA",
                                           type = "PPM",
                                           extrainfo = c(enrichment.score = meta[2]))
                           }, motif_meta, motif_list, SIMPLIFY = FALSE)

  } else {

    # minimal
    motif_A <- raw_lines[grepl("^A:", raw_lines)]
    motif_A <- lapply(motif_A, function(x) as.numeric(strsplit(x, "\\s+")[[1]][-1]))
    motif_C <- raw_lines[grepl("^C:", raw_lines)]
    motif_C <- lapply(motif_C, function(x) as.numeric(strsplit(x, "\\s+")[[1]][-1]))
    motif_G <- raw_lines[grepl("^G:", raw_lines)]
    motif_G <- lapply(motif_G, function(x) as.numeric(strsplit(x, "\\s+")[[1]][-1]))
    motif_T <- raw_lines[grepl("^T:", raw_lines)]
    motif_T <- lapply(motif_T, function(x) as.numeric(strsplit(x, "\\s+")[[1]][-1]))
    motif_names <- grep("^A:", raw_lines) - 1
    motif_names <- raw_lines[motif_names]

    motif_list <- mapply(function(name, a, c, g, t) {
                             motif <- matrix(c(a, c, g, t), nrow = 4,
                                             byrow = TRUE)
                             mot <- universalmotif_cpp(name = name,
                                            motif = motif,
                                            alphabet = "DNA",
                                            type = "PPM")
                             validObject(mot)
                             mot
                           }, motif_names, motif_A, motif_C, motif_G, motif_T,
                           SIMPLIFY = FALSE)

  }

  if (length(motif_list) == 1) motif_list <- motif_list[[1]]
  motif_list

}
