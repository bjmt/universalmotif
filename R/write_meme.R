#' Export motifs in MEME format.
#'
#' Convert motifs to minimal MEME format and write to file.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param file `character(1)` File name.
#' @param version `numeric(1)` MEME version.
#' @param bkg `numeric` Background letter frequencies. If missing, will use background
#'   frequencies from motif objects (if they are identical); else
#'   background frequencies will be set to freq = 1/length(alphabet)
#' @param strand `character` If missing, will use strand from motif objects (if identical);
#'   otherwise will default to "+ -"
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' transfac <- read_transfac(system.file("extdata", "transfac.txt",
#'                                     package = "universalmotif"))
#' write_meme(transfac, tempfile())
#'
#' @references
#'    \insertRef{meme}{universalmotif}
#'
#' @family write_motifs
#' @seealso [read_meme()]
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
write_meme <- function(motifs, file, version = 4, bkg, strand) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file, strand = args$strand),
                                 numeric(), c(FALSE, TRUE), "character")
  num_check <- check_fun_params(list(version = args$version), 1, FALSE, "numeric")
  all_checks <- c(char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, "PPM")
  if (!is.list(motifs)) motifs <- list(motifs)

  if (missing(strand)) {
    strand <- unique(vapply(motifs, function(x) x["strand"], character(1)))
    if (length(strand) > 1) strand <- "+ -"
    if (strand == "+-" || strand == "-+") strand <- "+ -"
  }

  alph <- unique(vapply(motifs, function(x) x["alphabet"], character(1)))
  if (length(alph) > 1) stop("all motifs must use the same alphabet")
  if (alph == "DNA")  {
    alph.2 <- "ACGT" 
    alph.length <- 4
  } else if (alph == "RNA") {
    alph.2 <- "ACGU"
    alph.length <- 4
  } else if (alph == "AA") {
    alph.2 <- "ACDEFGHIKLMNPQRSTVWY"
    alph.length <- 20
  } else stop("unknown alphabet")

  if (missing(bkg)) {
    bkg <- lapply(motifs, function(x) x["bkg"])
    bkgtest <- lapply(bkg, function(x) all(x == bkg[[1]]))
    if (all(unlist(bkgtest))) {
      bkg <- bkg[[1]]
    } else {
      if (alph == "DNA" || alph == "RNA") bkg <- c(0.25, 0.25, 0.25, 0.25)
      if (alph == "AA") bkg <- rep(1 / 20, 20)
    }
  }
  bkg.2 <- strsplit(alph.2, "")[[1]]

  lines_out <- c(paste("MEME version", version), "",
                 paste("ALPHABET=", alph.2), "",
                 paste("strands:", strand), "",
                 paste("Background letter frequencies"),
                 paste(bkg.2, bkg, collapse = " "), "")

  .write_meme <- function(motifs) {
  
    motif <- motifs
    lines_out <- paste("MOTIF", motif["name"])
    nsites <- motif["nsites"]
    nsites <- ifelse(length(nsites) == 0, 100, nsites)
    eval <- motif["eval"]
    eval <- ifelse(length(eval) == 0, 0, eval)
    lines_out <- c(lines_out,
                   paste("letter-probability matrix:", "alength=", alph.length,
                         "w=", ncol(motif["motif"]), "nsites=", nsites,
                         "E=", eval))
    for (i in seq_len(ncol(motif["motif"]))) {
      pos <- motif["motif"][, i]
      pos <- vapply(pos, function(x) format(x, nsmall = 6), character(1))
      pos <- paste("", pos, "", collapse = "")
      lines_out <- c(lines_out, pos)
    }

    c(lines_out, "")
  
  }

  lines_out <- c(lines_out, unlist(lapply(motifs, .write_meme)))

  writeLines(lines_out, con <- file(file))
  close(con)

  invisible(NULL)


}
