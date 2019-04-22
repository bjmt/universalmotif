#' Sequence-related utility functions.
#'
#' @param lets `character` A character vector where each element will be
#'    considered a single unit.
#' @param k `integer(1)` K-let size.
#' @param string `character(1)` A string of characters.
#' @param method `character(1)` Shuffling method. One of `c("euler", "linear",
#'    "markov")`. See [shuffle_sequences()].
#'
#' @return
#'    For [count_klets()]: A `data.frame` with columns `lets` and `counts`.
#'
#'    For [get_klets()]: A `character` vector of k-lets.
#'
#'    For [shuffle_string()]: A single `character` string`.
#'
#' @examples
#' #######################################################################
#' ## count_klets
#' ## Count k-lets for any string of characters
#' count_klets("GCAAATGTACGCAGGGCCGA", k = 2)
#' ## The default 'k' value (1) counts individual letters
#' count_klets("GCAAATGTACGCAGGGCCGA")
#'
#' #######################################################################
#' ## get_klets
#' ## Generate all possible k-lets for a set of characters
#' get_klets(c("A", "C", "G", "T"), 3)
#' ## Note that each element in 'lets' is considered a single unit;
#' ## see:
#' get_klets(c("AA", "B"), k = 2)
#'
#' #######################################################################
#' ## shuffle_string
#' ## Shuffle any string of characters
#' shuffle_string("ASDADASDASDASD", k = 2)
#'
#' @seealso [create_sequences()], [shuffle_sequences()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name utils-sequence
NULL

#' @rdname utils-sequence
#' @export
count_klets <- function(string, k = 1) {

  if (length(string) != 1 || length(k) != 1)
    stop("'string' and 'k' must be length 1")

  k <- as.integer(k)

  if (k < 1) 
    stop("k must be greater than 0")

  string <- safeExplode(string)
  counts <- letter_freqs(string, k, "freqs", FALSE, sort(unique(string)))

  counts$counts

}

# get_klets(lets, k = 1) --> see utils.cpp

#' @rdname utils-sequence
#' @export
shuffle_string <- function(string, k = 1, method = c("euler", "linear", "markov")) {

  method <- match.arg(method, c("euler", "linear", "markov"))
  k <- as.integer(k)
  if (length(string) != 1 || length(k) != 1)
    stop("'string' and 'k' must be length 1")

  if (k == 1) {

    shuffle_k1(string)

  } else if (k > 1) {

    switch(method,
           "euler" = shuffle_euler(string, k),
           "linear" = shuffle_linear(string, k),
           "markov" = shuffle_markov_any(string, k),
           stop("'method' must be one of 'euler', 'linear', 'markov'"))

  } else {

    stop("k must be greater than 0")

  }

}
