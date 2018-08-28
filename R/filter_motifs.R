#' Filter a list of \linkS4class{universalmotif} motifs.
#'
#' Filter motifs based on the contents of available [universalmotif-class]
#' slots.
#'
#' @param motifs `list` See [convert_motifs()] for acceptable
#'    formats.
#' @param name `character` Keep motifs by names.
#' @param altname `altname` Keep motifs by altnames.
#' @param family `family` Keep motifs by family.
#' @param organism `organism` Keep motifs by organism.
#' @param width `numeric(1)` Keep motifs with minimum width.
#' @param alphabet `character` Keep motifs by alphabet.
#' @param type `character` Keep motifs by type.
#' @param icscore `numeric(1)` Keep motifs with minimum total IC.
#' @param nsites `numeric(1)` Keep motifs with minimum number of target sites.
#' @param strand `character` Keeps motifs by strand.
#' @param pval `numeric(1)` Keep motifs by max P-value.
#' @param qval `numeric(1)` Keep motifs by max Q-value.
#' @param eval `numeric(1)` Keep motifs by max E-val.
#'
#' @return `list` Motifs.
#'
#' @examples
#' ## By minimum IC:
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.ic10 <- filter_motifs(jaspar, icscore = 10)
#'
#' ## By organism:
#' if (requireNamespace("MotifDb", quietly = TRUE)) {
#'   library(MotifDb)
#'   motifs <- convert_motifs(MotifDb)
#'   motifs <- filter_motifs(motifs, organism = c("Athaliana", "Mmusculus"))
#' }
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
filter_motifs <- function(motifs, name, altname, family, organism, width,
                          alphabet, type, icscore, nsites, strand, pval, qval,
                          eval) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(name = args$name, altname = args$altname,
                                      family = args$family, organism = args$organism,
                                      alphabet = args$alphabet, type = args$type,
                                      strand = args$strand),
                                 rep(0, 7), rep(TRUE, 7), "character")
  num_check <- check_fun_params(list(width = args$width, icscore = args$icscore,
                                     nsites = args$nsites, pval = args$pval,
                                     qval = args$qval, eval = args$eval),
                                rep(0, 6), rep(TRUE, 6), "numeric")
  all_checks <- c(char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------
  
  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  if (!missing(name)) {
    motif_names <- vapply(motifs, function(x) x["name"], character(1))
    motifs <- motifs[motif_names %in% name]
  }

  if (!missing(altname)) {
    motif_altnames <- sapply(motifs, function(x) x["altname"])
    motifs <- motifs[motif_altnames %in% altname]
  }

  if (!missing(family)) {
    motif_families <- sapply(motifs, function(x) x["family"])
    motifs <- motifs[motif_families %in% family]
  }

  if (!missing(organism)) {
    motif_organisms <- sapply(motifs, function(x) x["organism"])
    motifs <- motifs[motif_organisms %in% organism]
  }

  if (!missing(width)) {
    motif_widths <- vapply(motifs, function(x) ncol(x["motif"]), numeric(1))
    motifs <- motifs[motif_widths >= width]
  }

  if (!missing(alphabet)) {
    motif_alphabets <- vapply(motifs, function(x) x["alphabet"], character(1))
    motifs <- motifs[motif_alphabets %in% alphabet]
  }

  if (!missing(type)) {
    motif_types <- vapply(motifs, function(x) x["type"], character(1))
    motifs <- motifs[motif_types %in% type]
  }

  if (!missing(icscore)) {
    motif_icscores <- vapply(motifs, function(x) x["icscore"], numeric(1))
    motifs <- motifs[motif_icscores >= icscore]
  }

  if (!missing(nsites)) {
    motif_nsites <- apply(motifs, function(x) x["nsites"])
    motifs <- motifs[motif_nsites >= nsites]
  }

  if (!missing(strand)) {
    motif_strands <- vapply(motifs, function(x) x["strand"], character(1))
    motifs <- motifs[motif_strands %in% strand]
  }

  if (!missing(pval)) {
    motif_pvals <- sapply(motifs, function(x) x["pval"])
    motifs <- motifs[motif_pvals <= pval]
  }

  if (!missing(qval)) {
    motif_qvals <- sapply(motifs, function(x) x["qval"])
    motifs <- motifs[motif_qvals <= qval]
  }

  if (!missing(eval)) {
    motif_evals <- sapply(motifs, function(x) x["eval"])
    motifs <- motifs[motif_evals <= eval]
  }

  motifs <- .internal_convert(motifs, unique(CLASS_IN))

  motifs

}
