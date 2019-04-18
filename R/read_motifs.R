#' Import universalmotif formatted motifs.
#'
#' Import motifs created from [write_motifs()]. For optimal storage of
#' `universalmotif` class motifs, consider using [saveRDS()] and
#' [readRDS()]. Currently the `universalmotif` format is YAML-based.
#'
#' @param progress `logical(1)` Show progress.
#' @param BP `logical(1)` Allows for the use of \pkg{BiocParallel} within
#'    [read_motifs()]. See [BiocParallel::register()] to change the
#'    default backend.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @family read_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams read_cisbp
#' @export
read_motifs <- function(file, skip = 0, progress = FALSE, BP = FALSE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file),
                                 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(skip = args$skip), 1, FALSE, TYPE_NUM)
  logi_check <- check_fun_params(list(progress = args$progress,
                                      BP = args$BP),
                                 c(1, 1), c(FALSE, FALSE), TYPE_LOGI)
  all_checks <- c(char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  raw_lines <- readLines(con <- file(file))
  close(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]

  if (substr(raw_lines[1], 1, 24) == "# universalmotif version") {
    version <- strsplit(raw_lines[1], " ")[[1]][4]
    raw_lines <- raw_lines[-1]
  } else {
    message(wmsg("Could not detect universalmotif version, assuming most recent format"))
    version <- packageDescription("universalmotif")$Version
  }

  raw_lines <- raw_lines[!grepl("^#", raw_lines)]

  if (version < "1.1.67")
    motifs <- read_motifs_1_0(raw_lines, version)
  else
    motifs <- read_motifs_1_2(raw_lines, progress, BP)

  motifs

}

read_motifs_1_2 <- function(raw_lines, progress, BP) {

  total <- length(raw_lines)
  # starts <- which(raw_lines == "---")
  starts <- grep("^---([^-]|$)", raw_lines)
  if (length(starts) == 0) stop("Incorrectly formatted motif file")
  if (starts[length(starts)] == total) stop("Incorrectly formatted motif file")
  stops <- starts - 1
  stops <- c(stops[-1], total)
  starts <- starts + 1

  motifs <- mapply_(read_motifs_single_pre, starts, stops,
                    MoreArgs = list(raw_lines = raw_lines),
                    SIMPLIFY = FALSE, PB = progress, BP = BP)

  mots.check <- vapply(motifs, is.null, logical(1))
  if (all(mots.check))
    stop("Failed to read any motifs")
  else if (any(mots.check)) {
    failed.mots <- seq_along(motifs)[mots.check]
    warning(wmsg("Some motifs could not be read: ", paste0(failed.mots, collapse = ", ")))
    motifs <- motifs[!mots.check]
  }

  if (length(motifs) == 1) motifs[[1]] else motifs

}

read_motifs_single_pre <- function(start, stop, raw_lines) {
  tryCatch(read_motifs_single(raw_lines[start:stop]),
           error = function(e) return(NULL))
}

read_motifs_single <- function(mot) {

  mot <- yaml.load(mot)

  fields <- names(mot)

  mot[["bkg"]] <- unlist(mot[["bkg"]])

  if ("extrainfo" %in% fields)
    mot[["extrainfo"]] <- unlist(mot[["extrainfo"]])

  mot.mat <- mot[["motif"]]
  mot.mat <- strsplit(mot.mat, " ")
  n_row <- length(mot.mat[[1]])
  mot.mat <- matrix(as.numeric(unlist(mot.mat)), nrow = n_row)
  mot[["motif"]] <- mot.mat

  if ("multifreq" %in% fields) {

    mult.names <- names(mot$multifreq)
    mults <- vector("list", length(mult.names))
    for (i in seq_along(mults)) {
      mults[[i]] <- mot$multifreq[[mult.names[i]]]
      mults[[i]] <- strsplit(mults[[i]], " ")
      n_row <- length(mults[[i]][[1]])
      mults[[i]] <- matrix(as.numeric(unlist(mults[[i]])), nrow = n_row)
      colnames(mults[[i]]) <- seq_len(ncol(mults[[i]]))
      rownames(mults[[i]]) <- get_klets(rownames(motif@motif),
                                        log(nrow(mults[[i]]), nrow(motif@motif)))
    }

    names(mults) <- mult.names
    motif <- do.call(universalmotif_cpp, mot[fields != "multifreq"])
    motif@multifreq <- mults

  } else {

    motif <- do.call(universalmotif_cpp, mot)

  }

  validObject_universalmotif(motif)
  motif

}

#-------------------------------------------------------------------------------
# Pre-1.2.X

read_motifs_1_0 <- function(raw_lines, version) {

  names_search <- substr(raw_lines, 1, 5)
  motif_starts <- which(names_search == "name:")
  motif_stops <- c(motif_starts[-1] - 1, length(raw_lines))

  motifs <- mapply(function(x, y) motifs_to_list(raw_lines, x, y),
                   motif_starts, motif_stops, SIMPLIFY = FALSE)

  motifs <- lapply(motifs, motifs_to_umot, version = version)

  motifs.check <- mapply(check_list, motifs, seq_along(motifs), SIMPLIFY = TRUE)

  motifs <- motifs[motifs.check]

  if (length(motifs) == 0) stop("no motifs were read")
  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs

}

check_list <- function(motif, count) {
  msg <- validObject_universalmotif(motif, FALSE)
  if (length(msg) > 0) {
    warning(paste0("error parsing motif number ", count))
    FALSE
  } else TRUE
}

motifs_to_list <- function(raw_lines, motif_starts, motif_stops) {
  list(raw_lines[motif_starts:motif_stops])
}

motifs_to_umot <- function(motif, version) {

  motif <- motif[[1]]
  motif.split <- strsplit(motif, " ")

  which.name <- vapply(motif.split, function(x) x[1] == "name:", logical(1))
  which.name <- which(which.name)
  motif.name <- collapse2(motif.split[[which.name]][-1])

  which.altname <- vapply(motif.split, function(x) x[1] == "altname:", logical(1))
  which.altname <- which(which.altname)
  if (length(which.altname) == 1) {
    motif.altname <- collapse2(motif.split[[which.altname]][-1])
  } else motif.altname <- FALSE

  which.family <- vapply(motif.split, function(x) x[1] == "family:", logical(1))
  which.family <- which(which.family)
  if (length(which.family) == 1) {
    motif.family <- collapse2(motif.split[[which.family]][-1])
  } else motif.family <- FALSE

  which.organism <- vapply(motif.split, function(x) x[1] == "organism:", logical(1))
  which.organism <- which(which.organism)
  if (length(which.organism) == 1) {
    motif.organism <- collapse2(motif.split[[which.organism]][-1])
  } else motif.organism <- FALSE

  which.alphabet <- vapply(motif.split, function(x) x[1] == "alphabet:", logical(1))
  which.alphabet <- which(which.alphabet)
  if (length(which.alphabet) == 1) {
    motif.alphabet <- motif.split[[which.alphabet]][2]
  } else motif.alphabet <- FALSE

  which.type <- vapply(motif.split, function(x) x[1] == "type:", logical(1))
  which.type <- which(which.type)
  if (length(which.type) == 1) {
    motif.type <- motif.split[[which.type]][2]
  } else motif.type <- FALSE

  which.nsites <- vapply(motif.split, function(x) x[1] == "nsites:", logical(1))
  which.nsites <- which(which.nsites)
  if (length(which.nsites) == 1) {
    motif.nsites <- as.numeric(motif.split[[which.nsites]][2])
  } else motif.nsites <- FALSE

  which.pseudocount <- vapply(motif.split, function(x) x[1] == "pseudocount:", logical(1))
  which.pseudocount <- which(which.pseudocount)
  if (length(which.pseudocount) == 1) {
    motif.pseudocount <- as.numeric(motif.split[[which.pseudocount]][2])
  } else motif.pseudocount <- FALSE

  which.bkg <- vapply(motif.split, function(x) x[1] == "bkg:", logical(1))
  which.bkg <- which(which.bkg)
  if (length(which.bkg) == 1) {
    if (version < "1.1.57") {
      motif.bkg <- as.numeric(motif.split[[which.bkg]][-1])
    } else {
      motif.bkg <- vapply(motif.split[[which.bkg]][-1],
                          function(x) strsplit(x, "=")[[1]][2],
                          character(1))
      motif.bkg <- as.numeric(motif.bkg)
      names(motif.bkg) <- vapply(motif.split[[which.bkg]][-1],
                                 function(x) strsplit(x, "=")[[1]][1],
                                 character(1))
    }
  } else motif.bkg <- FALSE

  which.bkgsites <- vapply(motif.split, function(x) x[1] == "bkgsites:", logical(1))
  which.bkgsites <- which(which.bkgsites)
  if (length(which.bkgsites) == 1) {
    motif.bkgsites <- as.numeric(motif.split[[which.bkgsites]][2])
  } else motif.bkgsites <- FALSE

  which.strand <- vapply(motif.split, function(x) x[1] == "strand:", logical(1))
  which.strand <- which(which.strand)
  if (length(which.strand) == 1) {
    motif.strand <- collapse2(motif.split[[which.strand]][-1])
  } else motif.strand <- FALSE
  
  which.pval <- vapply(motif.split, function(x) x[1] == "pval:", logical(1))
  which.pval <- which(which.pval)
  if (length(which.pval) == 1) {
    motif.pval <- as.numeric(motif.split[[which.pval]][2])
  } else motif.pval <- FALSE

  which.qval <- vapply(motif.split, function(x) x[1] == "qval:", logical(1))
  which.qval <- which(which.qval)
  if (length(which.qval) == 1) {
    motif.qval <- as.numeric(motif.split[[which.qval]][2])
  } else motif.qval <- FALSE

  which.eval <- vapply(motif.split, function(x) x[1] == "eval:", logical(1))
  which.eval <- which(which.eval)
  if (length(which.eval) == 1) {
    motif.eval <- as.numeric(motif.split[[which.eval]][2])
  } else motif.eval <- FALSE

  which.extrainfo <- vapply(motif.split, function(x) x[1] == "extrainfo:", logical(1))
  extrainfo.start <- which(which.extrainfo)

  which.motif <- vapply(motif.split, function(x) x[1] == "motif:", logical(1))
  motif.start <- which(which.motif)
  which.multifreq <- vapply(motif.split, function(x) x[1] == "multifreq:", logical(1))
  multifreq.start <- which(which.multifreq)

  if (length(extrainfo.start) > 0) {
    extrainfo.all <- motif.split[(extrainfo.start + 1):(motif.start - 1)]
    extrainfo.names <- lapply(extrainfo.all, function(x) x[2])
    extrainfo.content <- lapply(extrainfo.all, function(x) collapse2(x[-c(1:2)]))
    motif.extrainfo <- mapply(function(x, y) { names(x) <- y; x },
                              extrainfo.content, extrainfo.names,
                              SIMPLIFY = TRUE)
  } else motif.extrainfo <- FALSE

  if (length(multifreq.start) == 0) {
    motif.motif <- parse_matrix(motif.split[(motif.start + 1):length(motif.split)])
    motif.multifreq <- FALSE
  } else {
    motif.motif <- parse_matrix(motif.split[(motif.start + 1):(multifreq.start - 1)])
    motif.multifreq <- parse_multi(motif.split[(multifreq.start + 1):length(motif.split)])
  }

  make_umot(motif.name, motif.altname, motif.family, motif.organism, motif.alphabet,
            motif.type, motif.nsites, motif.pseudocount, motif.bkg, motif.bkgsites,
            motif.strand, motif.pval, motif.qval, motif.eval, motif.motif,
            motif.multifreq, motif.extrainfo)

}

make_umot <- function(motif.name, motif.altname, motif.family, motif.organism,
                      motif.alphabet, motif.type, motif.nsites, motif.pseudocount,
                      motif.bkg, motif.bkgsites, motif.strand, motif.pval,
                      motif.qval, motif.eval, motif.motif, motif.multifreq,
                      motif.extrainfo) {

  args <- list(name = motif.name, altname = motif.altname, family = motif.family,
               organism = motif.organism, alphabet = motif.alphabet,
               nsites = motif.nsites, pseudocount = motif.pseudocount,
               bkg = motif.bkg, bkgsites = motif.bkgsites, strand = motif.strand,
               pval = motif.pval, qval = motif.qval, eval = motif.eval,
               motif = motif.motif, multifreq = motif.multifreq,
               extrainfo = motif.extrainfo)
  args <- args[!vapply(args, isFALSE, logical(1))]

  do.call(universalmotif_cpp, args)

}

parse_matrix <- function(lines) {
  lines <- lapply(lines, function(x) as.numeric(x[-1]))
  lines <- lines[vapply(lines, length, integer(1)) > 0]
  matrix(do.call(c, lines), byrow = TRUE, nrow = length(lines))
}

parse_multi <- function(lines) {

  mult.starts <- which(vapply(lines, function(x) x[1] == ">", logical(1)))
  mult.ends <- c(mult.starts[-1] - 1, length(lines))

  mults <- mapply(function(x, y) motifs_to_list(lines, x, y),
                  mult.starts + 1, mult.ends, SIMPLIFY = FALSE)
  mults <- lapply(mults, parse_matrix)
  names(mults) <- vapply(mult.starts, function(x) lines[[x]][2])

  mults

}

collapse2 <- function(char) {
  if (length(char) == 0) character(0)
  else paste0(char, collapse = " ")
}
