#' Import MEME motifs.
#'
#' Import MEME formatted motifs, as well as original motif sequences. See
#' \url{http://meme-suite.org/doc/meme-format.html}. Both 'full' and 'minimal'
#' formats are supported.
#'
#' @param file `character(1)` File name.
#' @param skip `numeric(1)` If not zero, will skip however many desired lines in the
#'    file before starting to read.
#' @param readsites `logical(1)` If `TRUE`, the motif sites will be read as well.
#'
#' @return `list` [universalmotif-class] objects. If `readsites = TRUE`, a list
#'    comprising of a sub-list of motif objects and a sub-list of 
#'    motif sites will be returned.
#'
#' @examples
#' meme.minimal <- read_meme(system.file("extdata", "meme_minimal.txt",
#'                                       package = "universalmotif"))
#' meme.full <- read_meme(system.file("extdata", "meme_full.txt",
#'                                    package = "universalmotif"))
#'
#' @references
#'    \insertRef{meme}{universalmotif}
#'
#' @family read_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
read_meme <- function(file, skip = 0, readsites = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file),
                                 1, FALSE, "character")
  num_check <- check_fun_params(list(skip = args$skip), 1, FALSE, "numeric")
  logi_check <- check_fun_params(list(readsites = args$readsites),
                                 1, FALSE, "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]
  raw_lines <- raw_lines[!grepl("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", raw_lines)]
  raw_lines <- raw_lines[!grepl("------------------------", raw_lines)]

  alph <- raw_lines[grepl("^ALPHABET=", raw_lines)]
  alph <- strsplit(alph, "\\s+")[[1]][2]
  alph.len <- nchar(alph)
  if (alph == "ACGT") {
    alph <- "DNA"
  } else if (alph == "ACGU") {
    alph <- "RNA"
  } else if (alph == "ACDEFGHIKLMNPQRSTVWY") {
    alph <- "AA"
  } else alph <- "custom"
  strands <- raw_lines[grepl("^strands:", raw_lines)]
  if (length(strands) > 0) {
    strands <- strsplit(strands, "\\s+")[[1]][-1]
  } else {
    strands <- c("+", "-")
  }
  if (all(c("+", "-") %in% strands)) {
    strands <- "+-"
  } 
  bkg <- raw_lines[grep("^Background letter frequencies", raw_lines) + 1]
  bkg <- strsplit(bkg, "\\s+")[[1]]
  bkg <- as.numeric(bkg[seq_len(length(bkg)) %% 2 == 0])

  motif_meta <- grep("^letter-probability matrix:", raw_lines)
  motif_names <- grep("^MOTIF ", raw_lines)
  # motif_names <- motif_meta - 1
  motif_names <- lapply(raw_lines[motif_names], function(x) {
                            x <- strsplit(x, "\\s+")[[1]]
                            if (x[1] == "") x[3] else x[2]
                          })
  motif_starts <- motif_meta + 1
  motif_stops <- sapply(raw_lines[motif_meta],
                        function(x) strsplit(x, "\\s+")[[1]][6])
  motif_stops <- motif_meta + as.numeric(motif_stops)

  motif_meta <- lapply(raw_lines[motif_meta],
                         function(x) {
                           x <- strsplit(x, "\\s+")[[1]]
                           c(nsites = as.numeric(x[8]),
                             eval = as.numeric(x[10]))
                         })
  motif_list <- mapply(function(x, y) {
                           z <- raw_lines[x:y]
                           z <- sapply(z, function(x) strsplit(x, "\\s+")[[1]])
                           z <- suppressWarnings(as.numeric(z))
                           z <- z[!is.na(z)]
                         }, motif_starts, motif_stops, SIMPLIFY = FALSE)

  motif_list <- mapply(function(x, y, z) {
                          mot <- universalmotif_cpp(name = x,
                                         type = "PPM",
                                         nsites = y[1],
                                         eval = y[2],
                                         bkg = bkg,
                                         alphabet = alph,
                                         strand = strands,
                                         motif = t(matrix(z, ncol = alph.len,
                                                          byrow = TRUE)))
                          msg <- validObject_universalmotif(mot)
                          if (length(msg) > 0) stop(msg) else mot
                         }, motif_names, motif_meta, motif_list,
                         SIMPLIFY = FALSE)

  if (length(motif_list) == 1) motif_list <- motif_list[[1]]

  if (readsites) {
    if (is.list(motif_list)) 
      mot.names <- vapply(motif_list, function(x) x["name"], character(1))
    else
      mot.names <- motif_list["name"]
    block.starts <- sapply(mot.names,
                           function(x) grep("in BLOCKS format", raw_lines))
    if (length(block.starts) == 0) {
      warning("could not find BLOCKS formatted motifs in MEME file")
      motif_list <- list(motifs = motif_list, sites = NULL)
    } else {
      block.len <- vapply(block.starts,
                          function(x) strsplit(raw_lines[x + 1], "seqs=")[[1]][2],
                          character(1))
      block.len <- as.numeric(block.len)
      block.starts <- block.starts + 2
      block.stops <- block.starts +  block.len - 1
      blocks <- mapply(function(x, y) read.table(text = raw_lines[x:y],
                                                  stringsAsFactors = FALSE),
                        block.starts, block.stops,
                        SIMPLIFY = FALSE)
      sites <- lapply(blocks, function(x) x$V4)
      site.names <- lapply(blocks, function(x) x$V1)
      if (alph == "DNA") {
        sites <- lapply(sites, DNAStringSet)
      } else if (alph == "RNA") {
        sites <- lapply(sites, RNAStringSet)
      } else if (alph == "AA") {
        sites <- lapply(sites, AAStringSet)
      } else {
        sites <- lapply(sites, BStringSet)
      }
      sites <- mapply(function(x, y) {names(x) <- y; x},
                        sites, site.names,
                        SIMPLIFY = FALSE)
      if (length(sites) == 1) sites <- sites[[1]]
      if (is.list(sites) && is.list(motif_list))  # TODO: this is a bug..
        if (length(sites) != length(motif_list))
          sites <- sites[seq_len(length(motif_list))]
      motif_list <- list(motifs = motif_list, sites = sites)
    }
  }

  motif_list

}
