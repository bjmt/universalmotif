#' Trim motifs.
#'
#' @param motifs Motif object.
#' @param IC_cutoff Numeric. Minimum allowed information content.
#'
#' @return Motifs Motif object or list.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.trimmed <- trim_motifs(jaspar)
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
trim_motifs <- function(motifs, IC_cutoff = 0.25) {

  if (class(motifs) == "list") {
    motifs <- lapply(motifs, function(x) trim_motifs(x, IC_cutoff = IC_cutoff))
    return(motifs)
  }

  CLASS_IN <- .internal_convert(motifs)
  motif <- convert_motifs(motifs)

  motif_scores <- apply(motif["motif"], 2, position_icscore,
                        bkg = motif["bkg"], type = motif["type"],
                        pseudoweight = motif["pseudoweight"],
                        nsites = motif["nsites"])
  to_cut <- rep(TRUE, length(motif_scores))
  for (i in seq_along(motif_scores)) {
    if (motif_scores[i] < IC_cutoff) to_cut[i] <- FALSE else break
  }
  for (i in rev(seq_along(motif_scores))) {
    if (motif_scores[i] < IC_cutoff) to_cut[i] <- FALSE else break
  }
  motif@motif <- as.matrix(motif@motif[, to_cut])
  if (ncol(motif@motif) == 0) return(NULL)
  motif <- universalmotif(name = motif["name"], altname = motif["altname"],
                          family = motif["family"], motif = motif["motif"],
                          alphabet = motif["alphabet"], type = motif["type"],
                          organism = motif["organism"],
                          nsites = motif["nsites"],
                          pseudoweight = motif["pseudoweight"],
                          bkg = motif["bkg"], bkgsites = motif["bkgsites"],
                          strand = motif["strand"], pval = motif["pval"],
                          qval = motif["qval"], eval = motif["eval"],
                          extrainfo = motif["extrainfo"])

  motif <- .internal_convert(motif, CLASS_IN)
  motif

}
