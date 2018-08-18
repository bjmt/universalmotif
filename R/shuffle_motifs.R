#' Shuffle motifs by column.
#'
#' Given a set of motifs, shuffle the columns between them. Currently does not
#' support keeping the 'multifreq' slot. Only the 'bkg', 'nsites', 'strand',
#' and 'bkgsites' slots will be preserved. Uses the same shuffling methods
#' as \code{\link{shuffle_sequences}}.
#'
#' @param motifs See \code{\link{convert_motifs}} for acceptable formats.
#' @param k \code{numeric(1)} K-let size.
#' @param method \code{character(1)} One of \code{c('linear', 'random')}.
#'    See details.
#' @param leftovers \code{character(1)} For \code{method = 'random'}. One of
#'    \code{c('asis', 'first', 'split', 'discard')}. See details.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Motifs. See \code{\link{convert_motifs}} for available output
#'    formats.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso \code{\link{shuffle_sequences}}
#' @export
shuffle_motifs <- function(motifs, k = 2, method = "linear",
                           leftovers = "asis", BPPARAM = SerialParam()) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(method = args$method,
                                      leftovers = args$leftovers),
                                 numeric(), logical(), "character")
  num_check <- check_fun_params(list(k = args$k), 1, FALSE, "numeric")
  s4_check <- check_fun_params(list(BPPARAM = args$BPPARAM), numeric(), FALSE, "S4")
  all_checks <- c(char_check, num_check, s4_check)
  all_checks <- paste(all_checks, collapse = "\n")
  if (length(all_checks) > 0 && all_checks[1] != "") stop(c("\n", all_checks))
  #---------------------------------------------------------

  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  if (!is.list(motifs)) motifs <- list(motifs)

  motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)
  mot.alphs <- vapply(motifs, function(x) x["alphabet"], character(1))
  if (length(unique(mot.alphs)) > 1)
    stop("all motifs must share the same alphabet")

  mot.lens <- vapply(motifs, function(x) ncol(x["motif"]), numeric(1))
  mot.mats <- lapply(motifs, function(x) x["motif"])

  mot.cols <- do.call(cbind, mot.mats)
  col.order <- seq_len(ncol(mot.cols))
  if (k == 1) {
    new.order <- sample(col.order, ncol(mot.cols))
  } else {
    if (method == "linear") {
      mot.cols2 <- as.character(col.order)
      new.order <- shuffle_linear(mot.cols2, k, mode = 2)
      new.order <- as.numeric(new.order)
    } else if (method == "random") {
      mot.cols2 <- as.character(col.order)
      new.order <- shuffle_random(mot.cols2, k, leftovers,
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
