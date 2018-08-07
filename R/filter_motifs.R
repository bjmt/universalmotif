#' Filter a list of \linkS4class{universalmotif} motifs.
#'
#' Filter motifs based on the contents of available \linkS4class{universalmotif}
#' slots.
#'
#' @param motifs List of motifs.
#' @param name Keep motifs by names.
#' @param altname Keep motifs by altnames.
#' @param family Keep motifs by family.
#' @param organism Keep motifs by organism.
#' @param width Keep motifs with minimum width.
#' @param alphabet Keep motifs by alphabet.
#' @param type Keep motifs by type.
#' @param icscore Keep motifs with minimum total IC.
#' @param nsites Keep motifs with minimum number of target sites.
#' @param strand Keeps motifs by strand.
#' @param pval Keep motifs by max P-value.
#' @param qval Keep motifs by max Q-value.
#' @param eval Keep motifs by max E-val.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return List of motifs.
#'
#' @examples
#' # By minimum IC:
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.ic10 <- filter_motifs(jaspar, icscore = 10)
#'
#' # By organism:
#' library(MotifDb)
#' motifs <- convert_motifs(MotifDb)
#' motifs <- filter_motifs(motifs, organism = c("Athaliana", "Mmusculus"))
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
filter_motifs <- function(motifs, name, altname, family, organism, width,
                          alphabet, type, icscore, nsites, strand, pval, qval,
                          eval, BPPARAM = SerialParam()) {
  
  CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)

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

  motifs <- .internal_convert(motifs, unique(CLASS_IN), BPPARAM = BPPARAM)

  motifs

}
