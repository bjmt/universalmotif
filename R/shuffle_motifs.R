#' Shuffle motifs by column.
#'
#' Given a set of motifs, shuffle the columns between them. Currently does not
#' support keeping the 'multifreq' slot. Only the 'bkg', 'nsites', 'strand',
#' and 'bkgsites' slots will be preserved. Uses the same shuffling methods
#' as [shuffle_sequences()]. When shuffling more than one motif, they are
#' shuffled together.
#'
#' @param motifs See [convert_motifs()] for acceptable formats.
#' @param k `numeric(1)` K-let size.
#' @param method `character(1)` One of `c('linear', 'random')`.
#'    Only relevant if `k > 1`. See details. The `'random'` method will
#'    be removed in the next minor version.
#' @param leftovers `character(1)` For \code{method = 'random'}. One of
#'    `c('asis', 'first', 'split', 'discard')`. See details.
#'
#' @return Motifs. See [convert_motifs()] for available output
#'    formats.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [shuffle_sequences()]
#' @export
shuffle_motifs <- function(motifs, k = 2, method = "linear",
                           leftovers = "asis") {

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  if (!method %in% c("linear", "random")) {
    method_check <- paste0(" * Incorrect 'shuffle.method': expected `linear` ",
                           "or `random`; got `", method, "`")
    method_check <- wmsg2(method_check, 4, 2)
    all_checks <- c(all_checks, method_check)
  }
  if (!leftovers %in% c("asis", "first", "split", "discard")) {
    leftovers_check <- paste0(" * Incorrect 'shuffle.leftovers': expected `asis`,",
                              " `first`, `split` or `discard`; got `",
                              leftovers, "`")
    leftovers_check <- wmsg2(leftovers_check, 4, 2)
    all_checks <- c(all_checks, leftovers_check)
  }
  char_check <- check_fun_params(list(method = args$method,
                                      leftovers = args$leftovers),
                                 numeric(), logical(), "character")
  num_check <- check_fun_params(list(k = args$k), 1, FALSE, "numeric")
  all_checks <- c(all_checks, char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (is(motifs, "universalmotif")) undo.list <- TRUE else undo.list <- FALSE
  if (is.list(motifs)) CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  else CLASS_IN <- .internal_convert(motifs)
  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  motifs <- convert_type_internal(motifs, "PPM")
  mot.alphs <- vapply(motifs, function(x) x@alphabet, character(1))
  if (length(unique(mot.alphs)) > 1)
    stop("all motifs must share the same alphabet")

  mot.lens <- vapply(motifs, function(x) ncol(x@motif), numeric(1))
  mot.mats <- lapply(motifs, function(x) x@motif)

  mot.cols <- do.call(cbind, mot.mats)
  col.order <- seq_len(ncol(mot.cols))
  if (k == 1) {
    new.order <- sample(col.order, ncol(mot.cols))
  } else {
    switch(method,
      "linear" = {
        new.order <- shuffle_linear(col.order, k, mode = 2)
        new.order <- as.numeric(new.order)
      },
      "random" = {
        warning("The 'random' method option will be removed in the next minor version update",
                immediate. = TRUE)
        new.order <- shuffle_random(col.order, k, leftovers, mode = 2)
        new.order <- as.numeric(new.order)
      },
      stop("only 'linear' and 'random' are supported")
    )
  }
  mot.cols <- mot.cols[, new.order]

  mot.offsets <- cumsum(c(0, mot.lens[-length(mot.lens)]))

  # new.mats <- vector("list", length(mot.mats))
  # for (i in seq_along(mot.mats)) {
    # new.mats[[i]] <- mot.cols[, seq_len(mot.lens[i]) + mot.offsets[i]]
  # }
  new.mats <- lapply(seq_along(mot.mats),
                     function(x) mot.cols[, (1 + mot.offsets[x]):(mot.lens[x] + mot.offsets[x])])

  new.motifs <- mapply(shuffle_new_mot, new.mats, motifs,
                         SIMPLIFY = FALSE)

  new.motifs <- .internal_convert(new.motifs, unique(CLASS_IN))
  if (undo.list && is.list(new.motifs)) new.motifs <- new.motifs[[1]]
  new.motifs

}

shuffle_new_mot <- function(new.mat, motif) {

  mot <- universalmotif_cpp(motif = new.mat, alphabet = motif@alphabet,
                            bkg = motif@bkg, bkgsites = motif@bkgsites,
                            nsites = motif@nsites, strand = motif@strand,
                            name = collapse_cpp(c(motif@name, " [shuffled]")))

  validObject_universalmotif(mot)
  mot

}
