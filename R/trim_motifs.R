#' Trim motifs.
#'
#' Remove edges of a motif with low information content.
#'
#' @param motifs Motif object.
#' @param min.ic Numeric. Minimum allowed information content.
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
trim_motifs <- function(motifs, min.ic = 0.25, BPPARAM = SerialParam()) {

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  if (!is.list(motifs)) motifs <- list(motifs)
  mot.names <- vapply(motifs, function(x) x["name"], character(1))

  motifs <- bplapply(motifs,
                     function(x) {
                       y <- x["nsites"]
                       if (length(y) == 0) x["nsites"] <- 100
                       x
                     }, BPPARAM = BPPARAM)

  mot.mats <- bplapply(motifs, function(x) x["motif"], BPPARAM = BPPARAM)

  mot.mats.k <- bplapply(motifs, function(x) x["multifreq"], BPPARAM = BPPARAM)

  mot.scores <- bplapply(motifs,
                         function(x) {
                          apply(x["motif"], 2, position_icscoreC,
                                bkg = x["bkg"], type = x["type"],
                                pseudocount = x["pseudocount"],
                                nsites = x["nsites"])
                         }, BPPARAM = BPPARAM)

  new.mats <- bpmapply(function(x, y) trim_motif_internal(x, y, min.ic),
                       mot.mats, mot.scores, BPPARAM = BPPARAM,
                       SIMPLIFY = FALSE)

  new.mats.k <- bpmapply(function(x, y) {
                        if (length(x) > 0) {
                          lapply(x, function(z) trim_motif_internal(z, y, min.ic))
                        } else list()
                       }, mot.mats.k, mot.scores, BPPARAM = BPPARAM,
                       SIMPLIFY = FALSE)

  motifs <- bpmapply(function(x, y, z) {
                        if (length(x) == 0) return(NULL)
                        z@motif <- x
                        z@multifreq <- y
                        z
                       }, new.mats, new.mats.k, motifs, BPPARAM = BPPARAM,
                       SIMPLIFY = FALSE)

  dont_keep <- vapply(motifs, is.null, logical(1))
  num_bar <- which(dont_keep)
  if (length(num_bar) > 0) {
    if (length(num_bar) == length(mot.names)) {
      stop("All motifs were completely trimmed")
    }
    message("The following motifs were completely trimmed: ",
            mot.names[num_bar])
    return(invisible(NULL))
  }
  
  motifs <- motifs[!dont_keep]

  motifs <- bplapply(motifs,
                     function(x) {
                       alph <- x@alphabet
                       type <- x@type
                       pseudo <- x@pseudocount
                       mat <- x@motif
                       bkg <- x@bkg
                       nsites <- x@nsites
                       ic <- sum(apply(mat, 2, position_icscoreC, bkg = bkg,
                                 type = type, pseudocount = pseudo,
                                 nsites = nsites))
                       x@icscore <- ic
                       if (alph %in% c("DNA", "RNA")) {
                         consensus <- apply(mat, 2, get_consensusC,
                                            alphabet = alph, type = type,
                                            pseudocount = pseudo)
                         colnames(mat) <- consensus
                         x@motif <- mat
                         x@consensus <- paste(consensus, collapse = "")
                       } else if (alph == "AA") {
                         consensus <- apply(mat, 2, get_consensusAAC,
                                            type = type, pseudocount = pseudo)
                         colnames(mat) <- consensus
                         x@motif <- mat
                         x@consensus <- paste(consensus, collapse = "")
                       }
                       msg <- validObject_universalmotif(x)
                       if (length(msg) > 0) stop(msg)
                       x
                     }, BPPARAM = BPPARAM)

  if (length(motifs) == 1) motifs[[1]] else motifs

}
