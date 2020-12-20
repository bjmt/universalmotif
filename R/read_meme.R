#' Import MEME motifs.
#'
#' Import MEME formatted motifs, as well as original motif sequences. See
#' \url{http://meme-suite.org/doc/meme-format.html}. Both 'full' and 'minimal'
#' formats are supported. DREME and STREME motifs can also be imported, but note
#' that `readsites` and `readsites.meta` arguments do nothing.
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
#' @details
#' Please note that the typical number precision limit in R is around 1e-308.
#' This means that motif P-values in MEME files below this limit are rounded
#' automatically to 0. To get around this, the E-value is also stored as a
#' string in the `extrainfo` slot. If you require a numeric value for analysis,
#' use the [log_string_pval()] function to get the log of the string-formatted
#' p-value.
#'
#' @examples
#' meme.minimal <- read_meme(system.file("extdata", "meme_minimal.txt",
#'                                       package = "universalmotif"))
#' meme.full <- read_meme(system.file("extdata", "meme_full.txt",
#'                                    package = "universalmotif"))
#' ## Get numeric p-value:
#' log_string_pval(meme.minimal[[1]]["extrainfo"]["eval.string"])
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
  raw_lines <- raw_lines[!grepl("^#", raw_lines)]
  raw_lines <- raw_lines[!grepl("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", raw_lines)]
  raw_lines <- raw_lines[!grepl("------------", raw_lines)]

  check_meme_version(raw_lines)

  alph <- get_meme_alph(raw_lines)
  alph.len <- get_meme_alph_len(alph)
  alph.split <- switch(alph, DNA = DNA_BASES, RNA = RNA_BASES, AA = AA_STANDARD2,
      safeExplode(alph))

  strands <- raw_lines[grepl("^strands:", raw_lines)]
  if (length(strands) > 0) {
    strands <- strsplit(strands, "\\s+")[[1]][-1]
  } else {
    message("Could not find strand info, assuming +.")
    strands <- "+"
  }
  if (all(c("+", "-") %in% strands)) {
    strands <- "+-"
  }
  bkg.start <- grep("^Background letter frequencies", raw_lines)
  if (length(bkg.start)) {
    bkg.offset <- 1
    bkg <- raw_lines[bkg.start + bkg.offset]
    bkg <- strsplit(bkg, "\\s+")[[1]]
    bkg <- as.numeric(bkg[seq_len(length(bkg)) %% 2 == 0])
    while (length(bkg) < alph.len) {
      bkg.offset <- bkg.offset + 1
      bkg.tmp <- raw_lines[bkg.start + bkg.offset]
      bkg.tmp <- strsplit(bkg.tmp, "\\s+")[[1]]
      bkg.tmp <- as.numeric(bkg.tmp[seq_along(bkg.tmp) %% 2 == 0])
      bkg <- c(bkg, bkg.tmp)
    }
    if (anyNA(bkg))
      stop("Could not parse background frequencies, check that they match alphabet")
  } else {
    message("Could not find background, assuming uniform frequencies.")
    bkg <- rep(1 / length(alph.split), length(alph.split))
  }

  motif_meta <- grep("^letter-probability matrix:", raw_lines)
  motif_names_i <- grep("^MOTIF ", raw_lines)
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
  #                          c(nsites = as.numeric(x[8]),
  # # Add a C++ function here to convert long double to log10 values?
  #                            eval = as.numeric(x[10]))
                           c(nsites = x[8], eval = x[10])
                         })
  motif_list <- mapply(function(x, y) {
                           z <- trimws(raw_lines[x:y])
                           z <- sapply(z, function(x) strsplit(x, "\\s+")[[1]])
                           if (nrow(z) != alph.len)
                             stop("Alphabet length does not match motif length")
                           z <- z[order(alph.split, method = "radix"), ]
                           as.numeric(z)
                         }, motif_starts, motif_stops, SIMPLIFY = FALSE)

  motif_list <- mapply(function(x, y, z, x2) {
                          mot <- universalmotif_cpp(name = x,
                                           type = "PPM",
                                           altname = x2,
                                           nsites = as.numeric(y[1]),
                                           eval = as.numeric(y[2]),
                                           bkg = bkg,
                                           alphabet = alph,
                                           strand = strands,
                                           extrainfo = c(eval.string = unname(y[2])),
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

#' Returns type of MEME alphabet to use
#'
#' Used to deploy a switch() statement in get_meme_alph()
#'
#' @param raw_lines raw lines from .meme file
#'
#' @return
#'
#' @noRd
check_meme_alph_type <- function(raw_lines){
  if (any(grepl("^ALPHABET=", raw_lines))) {
    "default"
  } else if (any(grepl("^ALPHABET", raw_lines))) {
    "custom"
  } else {
    message("Could not find alphabet, assuming DNA.")
    "none"
  }
}

#' Grabs alphabet type for defined alphabet entries
#'
#' @param raw_lines raw lines from .meme file
#'
#' @return
#'
#' @noRd
get_default_meme_alph <- function(raw_lines) {
  alph <- raw_lines[grepl("^ALPHABET=", raw_lines)]
  alph <- strsplit(alph, "\\s+")[[1]][2]
  alph <- switch(alph, "ACGT" = "DNA", "ACGU" = "RNA",
                 "ACDEFGHIKLMNPQRSTVWY" = "AA", alph)
  return(alph)
}

#' Get custom alphabet defintion from .meme file
#'
#' If this ever runs, import will probably break downstream because custom
#' alphabets typically define a lookup table  which doesn't correspond to
#' alphabet entries in the final matrix. (ie so alph.len will no longer be the
#' correct matrix width)
#'
#' This is kind of a placeholder. So I short-circuit it with an informative
#' error about custom alphabets not being supported.
#'
#' @param raw_lines raw lines from .meme file
#'
#' @return
#'
#' @noRd
get_custom_meme_alph <- function(raw_lines) {

  i <- grep("^ALPHABET", raw_lines) + 1
  n <- grep("^END ALPHABET", raw_lines) - 1

  alph_lines <- raw_lines[i:n]

  alph_lines <- alph_lines[!grepl("=", alph_lines, fixed = TRUE)]
  alph_lines <- trimws(unlist(strsplit(alph_lines, "~", fixed = TRUE)))
  alph <- vapply(strsplit(alph_lines, " ", fixed = TRUE), function(x) x[1], character(1))

  alph.order <- c(LETTERS, letters, 0:9, "*", "-", ".")
  alph <- alph[order(match(alph, alph.order))]
  alph <- collapse_cpp(alph)

  switch(alph, "ACGT" = "DNA", "ACGU" = "RNA", "ACDEFGHIKLMNPQRSTVWY" = "AA", alph)

}

#' Return MEME alphabet string
#'
#' @param raw_lines raw lines from .meme file
#'
#' @return `character(1)` alphabet string
#'
#' @noRd
get_meme_alph <- function(raw_lines) {

  switch(check_meme_alph_type(raw_lines),
         default = get_default_meme_alph(raw_lines),
         custom = get_custom_meme_alph(raw_lines),
         none = "DNA")

}

#' Return alphabet length
#'
#' Uses lookup table for default alphabets, otherwise just counts number of
#' letters
#'
#' @param alph alphabet string from get_meme_alph
#'
#' @return `integer(1)` alphabet length
#'
#' @noRd
get_meme_alph_len <- function(alph) {

  switch(alph, "DNA" = 4L, "RNA" = 4L, "AA" = 20L, nchar(alph))

}

check_meme_version <- function(raw_lines) {
  ver <- raw_lines[grep("^MEME version", raw_lines)]
  if (!length(ver)) {
    warning(wmsg("Could not find MEME version. Note that this is considered",
        " a requirement by MEME programs."),
      call. = FALSE)
  } else{
    ver <- strsplit(ver, "\\s+")[[1]][3]
    ver <- strsplit(ver, "[.]")[[1]][1]
    ver <- as.numeric(ver)
    if (!ver %in% 4:5) {
      warning(wmsg("read_meme() is configured for MEME versions 4-5."),
        call. = FALSE)
    }
  }
}
