#' Sequence-related utility functions.
#'
#' @param lets `character` A character vector of letters.
#' @param k `integer(1)` K-let size.
#' @param string `character(1)` A string of letters.
#' @param method `character(1)` Shuffling method. One of `c("euer", "linear",
#'    "markov")`. See [shuffle_sequences()].
#'
#' @return
#'    For [count_klets()]: A `data.frame` with columns `lets` and `counts`.
#'
#'    For [get_klets()]: A `character` vector of k-lets.
#'
#'    For [shuffle_string()]: A single string of type `character`.
#'
#' @examples
#' #######################################################################
#' ## count_klets
#' ## Count k-lets for any string of characters
#' count_klets("GCAAATGTACGCAGGGCCGA", 2)
#'
#' #######################################################################
#' ## get_klets
#' ## Generate all possible k-lets for a set of letters
#' get_klets(c("A", "C", "G", "T"), 3)
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
count_klets <- function(string, k) {
  if (length(string) > 1)
    stop("please input a single string")
  k <- as.integer(k)
  if (k < 1) 
    stop("k must be greater than 0")
  string <- safeExplode(string)
  counts <- letter_freqs(string, k, "freqs", FALSE, sort(unique(string)))
  counts$counts
}

#' @rdname utils-sequence
#' @export
shuffle_string <- function(string, k = 1, method = c("euler", "linear", "markov")) {
  method <- match.arg(method, c("euler", "linear", "markov"))
  k <- as.integer(k)[1]
  if (length(string) != 1)
    stop("'string' must be a single string")
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
