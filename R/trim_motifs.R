#' Trim motifs.
#'
#' Remove edges of motifs with low information content.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param min.ic `numeric(1)` Minimum allowed information content. See
#'    [convert_type()] for a discussion on information content.
#'
#' @return Motifs See [convert_motifs()] for available output
#'    formats.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.trimmed <- trim_motifs(jaspar)
#'
#' @seealso [create_motif()], [convert_type()]
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
trim_motifs <- function(motifs, min.ic = 0.25) {

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list(min.ic = args$min.ic),
                                1, FALSE, "numeric")
  all_checks <- c(num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, "character")
  else CLASS_IN <- .internal_convert(motifs)

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  mot.names <- vapply(motifs, function(x) x["name"], character(1))

  motifs <- lapply(motifs,
                     function(x) {
                       y <- x["nsites"]
                       if (length(y) == 0) x["nsites"] <- 100
                       x
                     })

  mot.mats <- lapply(motifs, function(x) x["motif"])

  mot.mats.k <- lapply(motifs, function(x) x["multifreq"])

  mot.scores <- lapply(motifs,
                         function(x) {
                          apply(x["motif"], 2, position_icscoreC,
                                bkg = x["bkg"], type = x["type"],
                                pseudocount = x["pseudocount"],
                                nsites = x["nsites"])
                         })

  new.mats <- mapply(function(x, y) trim_motif_internal(x, y, min.ic),
                       mot.mats, mot.scores,
                       SIMPLIFY = FALSE)

  new.mats.k <- mapply(function(x, y) {
                        if (length(x) > 0) {
                          lapply(x, function(z) trim_motif_internal(z, y, min.ic))
                        } else list()
                       }, mot.mats.k, mot.scores,
                       SIMPLIFY = FALSE)

  motifs <- mapply(function(x, y, z) {
                        if (length(x) == 0) return(NULL)
                        z@motif <- x
                        z@multifreq <- y
                        z
                       }, new.mats, new.mats.k, motifs,
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

  motifs <- lapply(motifs,
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
                     })

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs <- .internal_convert(motifs, unique(CLASS_IN))
  motifs

}
