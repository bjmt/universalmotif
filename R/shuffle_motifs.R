#' Shuffle motifs by column.
#'
#' Given a set of motifs, shuffle the columns between them.
#'
#' @param motifs Motifs.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
# shuffle_motifs <- function(motifs, BPPARAM = SerialParam()) {
#
  # CLASS_IN <- vapply(motifs, .internal_convert, character(1))
#
  # motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  # motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)
#
  # mot.lens <- vapply(motifs, function(x) ncol(x["motif"]), numeric(1))
  # mot.mats <- lapply(motifs, function(x) x["motif"])
#
  # mot.cols <- do.call(cbind, mot.mats)
  # new.order <- sample(seq_len(ncol(mot.cols)), ncol(mot.cols))
  # mot.cols <- mot.cols[, new.order]
#
  # new.mats <- vector("list", length(mot.mats))
  # for (i in seq_along(mot.mats)) {
    # new.mats[[i]] <- mot.cols[,]
  # }
#
# }
