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
#' @param readsites.meta `logical(1)` If `readsites = TRUE`, then additionally
#'    read site positions and P-values.
#'
#' @return `list` [universalmotif-class] objects. If `readsites = TRUE`, a list
#'    comprising of a sub-list of motif objects and a sub-list of 
#'    motif sites will be returned. If `readsites.meta = TRUE`, then two
#'    additional list items will be present, one containing site positions
#'    and P-values, and another containing combined sequence p-values.
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
read_meme <- function(file, skip = 0, readsites = FALSE,
                      readsites.meta = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file),
                                 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(skip = args$skip), 1, FALSE, TYPE_NUM)
  logi_check <- check_fun_params(list(readsites = args$readsites,
                                      readsites.meta = args$readsites.meta),
                                 numeric(), logical(), TYPE_LOGI)
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]
  raw_lines <- raw_lines[!grepl("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", raw_lines)]
  raw_lines <- raw_lines[!grepl("------------", raw_lines)]

  alph <- raw_lines[grepl("^ALPHABET=", raw_lines)]
  alph <- strsplit(alph, "\\s+")[[1]][2]
  alph.len <- nchar(alph)
  alph <- switch(alph, "ACGT" = "DNA", "ACGU" = "RNA",
                 "ACDEFGHIKLMNPQRSTVWY" = "AA", alph)
  strands <- raw_lines[grepl("^strands:", raw_lines)]
  if (length(strands) > 0) {
    strands <- strsplit(strands, "\\s+")[[1]][-1]
  } else {
    strands <- "+"
  }
  if (all(c("+", "-") %in% strands)) {
    strands <- "+-"
  } 
  bkg <- raw_lines[grep("^Background letter frequencies", raw_lines) + 1]
  bkg <- strsplit(bkg, "\\s+")[[1]]
  bkg <- as.numeric(bkg[seq_len(length(bkg)) %% 2 == 0])

  motif_meta <- grep("^letter-probability matrix:", raw_lines)
  motif_names_i <- grep("^MOTIF ", raw_lines)
  # motif_names <- motif_meta - 1
  motif_names <- lapply(raw_lines[motif_names_i], function(x) {
                            x <- strsplit(x, "\\s+")[[1]]
                            if (x[1] == "") x[3] else x[2]
                          })
  motif_altnames <- lapply(raw_lines[motif_names_i], function(x) {
                            x <- strsplit(x, "\\s+")[[1]]
                            if (x[1] == "") x[4] else x[3]
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

  motif_list <- mapply(function(x, y, z, x2) {
                          mot <- universalmotif_cpp(name = x,
                                           type = "PPM",
                                           altname = x2,
                                           nsites = y[1],
                                           eval = y[2],
                                           bkg = bkg,
                                           alphabet = alph,
                                           strand = strands,
                                           motif = t(matrix(z, ncol = alph.len,
                                                            byrow = TRUE)))
                          validObject_universalmotif(mot)
                          mot
                         }, motif_names, motif_meta, motif_list, motif_altnames,
                         SIMPLIFY = FALSE)

  if (length(motif_list) == 1) motif_list <- motif_list[[1]]

  if (readsites) {
    if (is.list(motif_list)) 
      mot.names <- vapply(motif_list, function(x) x@name, character(1))
    else
      mot.names <- motif_list@name
    block.starts <- grep("in BLOCKS format", raw_lines)
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
      sites <- switch(alph, "DNA" = lapply(sites, DNAStringSet),
                      "RNA" = lapply(sites, RNAStringSet),
                      "AA" = lapply(sites, AAStringSet),
                      lapply(sites, BStringSet))
      sites <- mapply(function(x, y) {names(x) <- y; x},
                        sites, site.names,
                        SIMPLIFY = FALSE)
      names(sites) <- mot.names
      if (length(sites) == 1) sites <- sites[[1]]
      if (is.list(sites) && is.list(motif_list))  # TODO: this is a bug..
        if (length(sites) != length(motif_list))
          sites <- sites[seq_len(length(motif_list))]
      motif_list <- list(motifs = motif_list, sites = sites)
    }

    if (readsites.meta) {
      site.starts <- grep("sites sorted by position p-value", raw_lines)
      site.stops <- grep("block diagrams$", raw_lines)
      if (length(site.starts) == 0 || length(site.stops) == 0) {
        warning("Could not find site P-values in MEME file")
      } else {
        site.starts <- site.starts + 2
        site.stops <- site.stops - 1
        site.tables <- mapply(function(x, y) {
                                z <- raw_lines[x:y]
                                lapply(z, function(x) strsplit(x, "\\s+")[[1]])
                              }, site.starts, site.stops, SIMPLIFY = FALSE)
        col.seqname <- 1
        if (all(grepl("Strand", raw_lines[site.starts - 1]))) {
          col.pos <- 3
          col.pval <- 4
          col.seq <- 6
        } else {
          col.pos <- 2
          col.pval <- 3
          col.seq <- 5
        }
        site.tables <- lapply(site.tables,
                              function(x) {
                                z1 <- vapply(x, function(x) x[col.seqname], character(1))
                                z2 <- vapply(x, function(x) x[col.pos], character(1))
                                z3 <- vapply(x, function(x) x[col.pval], character(1))
                                z4 <- vapply(x, function(x) x[col.seq], character(1))
                                data.frame(Sequence = z1,
                                           Position = as.numeric(z2),
                                           Pvalue = as.numeric(z3),
                                           Site = z4,
                                           stringsAsFactors = FALSE)
                              })
        names(site.tables) <- mot.names
        if (length(site.tables) == 1) site.tables <- site.tables[[1]]
        motif_list <- list(motifs = motif_list$motifs,
                           sites = motif_list$sites,
                           sites.meta = site.tables)

        summ.start <- grep("Combined block diagrams: non-overlapping sites",
                           raw_lines)
        if (length(summ.start) == 0) {
          warning("Could not find combined P-values in MEME file")
        } else {
          summ.start <- summ.start + 2
          summ.raw <- raw_lines[summ.start:(length(raw_lines) - 2)]
          need.fix <- grep("\\", summ.raw, fixed = TRUE)
          if (length(need.fix) > 0) {
            need.fix2 <- need.fix + 1
            summ.raw[need.fix] <- vapply(summ.raw[need.fix],
                                         function(x) strsplit(x, "\\",
                                                              fixed = TRUE)[[1]],
                                         character(1))
            summ.raw[need.fix2] <- gsub("\\s+", "", summ.raw[need.fix2])
            summ.raw[need.fix] <- mapply(function(x, y) paste0(x, y),
                                         summ.raw[need.fix], summ.raw[need.fix2])
            summ.raw <- summ.raw[-need.fix2]
          }
          summ.tab <- read.table(text = summ.raw, stringsAsFactors = FALSE)
          colnames(summ.tab) <- c("Sequence", "Combined.Pvalue", "Diagram")
          motif_list <- list(motifs = motif_list$motifs,
                             sites = motif_list$sites,
                             sites.meta = motif_list$sites.meta,
                             sites.meta.combined = summ.tab)
        }
      }
    }

  } else if (readsites.meta) {
    warning("'readsites.meta' is not valid if 'readsites = FALSE'")
  }

  motif_list

}
