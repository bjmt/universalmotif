#' Trim motifs.
#'
#' Remove edges of motifs with low information content. Currently does
#' not trim `multifreq` representations.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param min.ic `numeric(1)` Minimum allowed information content. See
#'    [convert_type()] for a discussion on information content.
#' @param trim.from `character(1)` Control the direction of trimming.
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
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @export
trim_motifs <- function(motifs, min.ic = 0.25, trim.from = c("both", "left", "right")) {

  # TODO:
  #   - manual trimming?

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list(min.ic = args$min.ic),
                                1, FALSE, TYPE_NUM)
  all_checks <- c(num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  trim.from <- match.arg(trim.from, c("both", "left", "right"))
  trim.from <- switch(trim.from, "both" = 0, "left" = 1, "right" = 2)

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, "character")
  else CLASS_IN <- .internal_convert(motifs)

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)
  mot.names <- vapply(motifs, function(x) x@name, character(1))

  mot.mats <- lapply(motifs, function(x) x@motif)

  mot.mats.k <- lapply(motifs, function(x) x@multifreq)

  mot.scores <- lapply(motifs,
                         function(x) {
                          apply(x@motif, 2, position_icscoreC,
                                bkg = x@bkg[rownames(x@motif)],
                                type = x@type,
                                pseudocount = x@pseudocount,
                                nsites = if (length(x@nsites)) x@nsites else 100)
                         })

  cut.left <- vapply(mot.scores, function(x) {
    if (!trim.from %in% c(0, 1)) return(0)
    kept <- which(x >= min.ic)
    if (length(kept)) kept[1] - 1 else length(x)
  }, numeric(1))

  new.mats <- mapply(function(x, y) trim_motif_internal(x, y, min.ic, trim.from),
    mot.mats, mot.scores, SIMPLIFY = FALSE)

  new.mats.k <- mapply(
    function(x, y) {
      if (length(x) > 0) {
        lapply(x, function(z) trim_motif_internal(z, y, min.ic, trim.from))
      } else list()
     }, mot.mats.k, mot.scores,
     SIMPLIFY = FALSE
  )

  motifs <- mapply(
    function(x, y, z, cut) {
      if (length(x) == 0) return(NULL)
      z@motif <- x
      z@multifreq <- y
      trim_gapinfo(z, cut)
     }, new.mats, new.mats.k, motifs, cut.left,
     SIMPLIFY = FALSE
  )

  dont_keep <- vapply(motifs, is.null, logical(1))
  num_bar <- which(dont_keep)
  if (length(num_bar) > 0) {
    if (length(num_bar) == length(mot.names)) {
      stop("All motifs were completely trimmed")
    }
    message("The following motifs were completely trimmed: ",
            mot.names[num_bar])
  }

  motifs <- motifs[!dont_keep]

  motifs <- lapply(motifs,
                     function(x) {
                       x <- refresh_motif_summary(x)
                       validObject_universalmotif(x)
                       x
                     })

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs <- .internal_convert(motifs, unique(CLASS_IN))
  motifs

}

