#' Trim motifs.
#'
#' Remove edges of a motif with low information content.
#'
#' @param motifs Motif object.
#' @param IC_cutoff Numeric. Minimum allowed information content.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Motifs Motif object or list.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.trimmed <- trim_motifs(jaspar)
#'
#' @seealso \code{\link{create_motif}}
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
trim_motifs <- function(motifs, IC_cutoff = 0.25, BPPARAM = SerialParam()) {

  # TODO: too slow. Need to implement in c++.

  if (class(motifs) == "list") {
    motifs <- bplapply(motifs, function(x) trim_motifs(x, IC_cutoff = IC_cutoff),
                       BPPARAM = BPPARAM)
    tokeep <- vapply(motifs, function(x) !is.null(x), logical(1))
    if (FALSE %in% tokeep) {
      num_gone <- length(which(tokeep == FALSE))
      warning(num_gone, " motifs were completely trimmed and subsequently deleted")
    }
    return(motifs[tokeep])
  }

  CLASS_IN <- .internal_convert(motifs)
  motif <- convert_motifs(motifs, BPPARAM = BPPARAM)

  if (length(motif["nsites"]) == 0) nsites <- 100 else nsites <- motif["nsites"]

  motif_scores <- apply(motif["motif"], 2, position_icscoreC,
                        bkg = motif["bkg"], type = motif["type"],
                        pseudocount = motif["pseudocount"],
                        nsites = nsites)
  to_cut <- rep(TRUE, length(motif_scores))
  cut_left <- 0
  for (i in seq_along(motif_scores)) {
    if (motif_scores[i] < IC_cutoff) to_cut[i] <- FALSE else break
    cut_left <- cut_left + 1
  }
  cut_right <- 0
  for (i in rev(seq_along(motif_scores))) {
    if (motif_scores[i] < IC_cutoff) to_cut[i] <- FALSE else break
    cut_right <- cut_right + 1
  }
  motif_mat <- as.matrix(motif@motif[, to_cut])
  if (ncol(motif_mat) == 0) return(NULL)

  if (cut_right != 0) {
    cut_right <- ncol(motif_mat) - cut_right + 1
    if (cut_right != ncol(motif_mat)) cut_right <- cut_right:ncol(motif_mat)
  }
  multifreq <- motif@multifreq
  if (length(multifreq) > 0) {
    for (i in seq_along(multifreq)) {
      multifreq[[i]] <- multifreq[[i]][, -seq_len(cut_left)]
      multifreq[[i]] <- multifreq[[i]][, -cut_right]
      colnames(multifreq[[i]]) <- seq_len(ncol(multifreq[[i]]))
    }
  }

  if (!is.matrix(motif_mat)) motif_mat <- as.matrix(motif_mat)
  motif <- universalmotif_cpp(name = motif["name"], altname = motif["altname"],
                          family = motif["family"], motif = motif["motif"],
                          alphabet = motif["alphabet"], type = motif["type"],
                          organism = motif["organism"],
                          nsites = motif["nsites"],
                          pseudocount = motif["pseudocount"],
                          bkg = motif["bkg"], bkgsites = motif["bkgsites"],
                          strand = motif["strand"], pval = motif["pval"],
                          qval = motif["qval"], eval = motif["eval"],
                          extrainfo = motif["extrainfo"])
  if (length(multifreq) > 0) motif@multifreq <- multifreq

  msg <- validObject_universalmotif(motif)
  if (length(msg) > 0) stop(msg)

  motif <- .internal_convert(motif, CLASS_IN, BPPARAM = BPPARAM)
  motif

}
