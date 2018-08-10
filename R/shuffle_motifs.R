#' Shuffle motifs by column.
#'
#' Given a set of motifs, shuffle the columns between them. Currently does not
#' support keeping the 'multifreq' slot. Only the 'bkg', 'nsites', 'strand',
#' and 'bkgsites' slots will be preserved.
#'
#' @param motifs Motifs.
#' @param shuffle.k Numeric.
#' @param shuffle.method Character.
#' @param shuffle.leftovers Character.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
shuffle_motifs <- function(motifs, shuffle.k = 2, shuffle.method = "linear",
                           shuffle.leftovers = "asis", BPPARAM = SerialParam()) {

  CLASS_IN <- vapply(motifs, .internal_convert, character(1))

  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)
  mot.alphs <- vapply(motifs, function(x) x["alphabet"], character(1))
  if (length(unique(mot.alphs)) > 1)
    stop("all motifs must share the same alphabet")

  mot.lens <- vapply(motifs, function(x) ncol(x["motif"]), numeric(1))
  mot.mats <- lapply(motifs, function(x) x["motif"])

  mot.cols <- do.call(cbind, mot.mats)
  col.order <- seq_len(ncol(mot.cols))
  if (shuffle.k == 1) {
    new.order <- sample(col.order, ncol(mot.cols))
  } else {
    if (shuffle.method == "linear") {
      mot.cols2 <- as.character(col.order)
      new.order <- shuffle_linear(mot.cols2, shuffle.k, mode = 2)
      new.order <- as.numeric(new.order)
    } else if (shuffle.method == "random") {
      mot.cols2 <- as.character(col.order)
      new.order <- shuffle_random(mot.cols2, shuffle.k, shuffle.leftovers,
                                  mode = 2)
      new.order <- as.numeric(new.order)
    } else stop("only 'linear' and 'random' are supported")
  }
  mot.cols <- mot.cols[, new.order]

  mot.offsets <- cumsum(c(0, mot.lens[-length(mot.lens)]))

  new.mats <- vector("list", length(mot.mats))
  for (i in seq_along(mot.mats)) {
    new.mats[[i]] <- mot.cols[, seq_len(mot.lens[i]) + mot.offsets[i]]
  }

  new.motifs <- bpmapply(shuffle_new_mot, new.mats, motifs,
                         BPPARAM = BPPARAM, SIMPLIFY = FALSE)

  new.motifs <- .internal_convert(new.motifs, unique(CLASS_IN), BPPARAM = BPPARAM)
  new.motifs

}

shuffle_new_mot <- function(new.mat, motif) {

  mot <- universalmotif_cpp(motif = new.mat, alphabet = motif["alphabet"],
                            bkg = motif["bkg"], bkgsites = motif["bkgsites"],
                            nsites = motif["nsites"], strand = motif["strand"],
                            name = paste(motif["name"], "[shuffled]"))
  msg <- validObject_universalmotif(mot)
  if (length(msg) > 0) stop(msg)
  mot

}
