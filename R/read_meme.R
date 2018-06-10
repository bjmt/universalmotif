#' Import MEME motifs.
#'
#' @param file Character.
#' @param skip Numeric.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return List of universalmotif objects.
#'
#' @examples
#' meme.minimal <- read_meme(system.file("extdata", "meme_minimal.txt",
#'                                       package = "universalmotif"))
#' meme.full <- read_meme(system.file("extdata", "meme_full.txt",
#'                                    package = "universalmotif"))
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_meme <- function(file, skip = 0, BPPARAM = bpparam()) {

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]
  raw_lines <- raw_lines[!grepl("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", raw_lines)]
  raw_lines <- raw_lines[!grepl("------------------------", raw_lines)]

  alph <- raw_lines[grepl("^ALPHABET=", raw_lines)]
  alph <- strsplit(alph, "\\s+")[[1]][2]
  if (alph == "ACGT") {
    alph <- "DNA"
  } else if (alph == "ACGU") {
    alph <- "RNA"
  } else if (alph == "ACDEFGHIKLMNPQRSTVWY") {
    alph <- "AA"
  } else alph <- "custom"
  strands <- raw_lines[grepl("^strands:", raw_lines)]
  strands <- strsplit(strands, "\\s+")[[1]][-1]
  if (all(c("+", "-") %in% strands)) {
    strands <- "+-"
  } 
  bkg <- raw_lines[grep("^Background letter frequencies", raw_lines) + 1]
  bkg <- strsplit(bkg, "\\s+")[[1]]
  bkg <- as.numeric(bkg[seq_len(length(bkg)) %% 2 == 0])

  motif_meta <- grep("^letter-probability matrix:", raw_lines)
  motif_names <- motif_meta - 1
  motif_names <- bplapply(raw_lines[motif_names], function(x) {
                            x <- strsplit(x, "\\s+")[[1]]
                            if (x[1] == "") x[3] else x[2]
                          }, BPPARAM = BPPARAM)
  motif_starts <- motif_meta + 1
  motif_stops <- sapply(raw_lines[motif_meta],
                        function(x) strsplit(x, "\\s+")[[1]][6])
  motif_stops <- motif_meta + as.numeric(motif_stops)

  motif_meta <- bplapply(raw_lines[motif_meta],
                         function(x) {
                           x <- strsplit(x, "\\s+")[[1]]
                           c(nsites = as.numeric(x[8]),
                             eval = as.numeric(x[10]))
                         }, BPPARAM = BPPARAM)
  motif_list <- bpmapply(function(x, y) {
                           z <- raw_lines[x:y]
                           z <- sapply(z, function(x) strsplit(x, "\\s+")[[1]])
                           z <- suppressWarnings(as.numeric(z))
                           z <- z[!is.na(z)]
                         }, motif_starts, motif_stops, SIMPLIFY = FALSE,
                         BPPARAM = BPPARAM)

  motif_list <- bpmapply(function(x, y, z) {
                          universalmotif(name = x,
                                         type = "PPM",
                                         nsites = y[1],
                                         eval = y[2],
                                         bkg = bkg,
                                         alphabet = alph,
                                         strand = strands,
                                         motif = t(matrix(z, ncol = 4,
                                                          byrow = TRUE)))
                         }, motif_names, motif_meta, motif_list,
                         SIMPLIFY = FALSE, BPPARAM = BPPARAM)

  motif_list

}
