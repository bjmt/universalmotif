#' Calculate sequence complexity.
#'
#' Calculate sequence complexity using either the Wootton-Federhen, Trifonov,
#' or DUST algorithms.
#'
#' @param seqs \code{\link{XStringSet}} Input sequences.
#' @param window.size `numeric` Window size. If a number between 0 and 1 is
#'    provided, the value is calculated as the number multiplied by the sequence
#'    length. 
#' @param window.overlap `numeric` Overlap between windows. If a number
#'    between 0 and 1 is provided, the value is calculated as the number
#'    multiplied by the sequence length.
#' @param method `character(1)` Choose one of the available methods for calculating
#'    sequence complexity. See details.
#' @param trifonov.max.word.size `numeric(1)` The maximum word size within each window
#'    used to calculate complexity using `method = c("Trifonov", "TrifonovFast")`.
#'    In other words, the Trifonov method involves counting the number of possible
#'    different sub-words in a window at different sizes up to the values provided
#'    by this option. It also involves calculating the product of ever increasing
#'    sequences of numbers and so in order to reduce the computations involed
#'    this can be limited to a specific maximum sub-word size.
#' @param nthreads `numeric(1)` Run [sequence_complexity()] in parallel with `nthreads`
#'    threads. `nthreads = 0` uses all available threads.
#' @param return.granges `logical(1)` Return the results as a `GRanges` object.
#'    Requires the `GenomicRanges` package to be installed.
#'
#' @return `DataFrame`, `GRanges` with each row getting a complexity score for
#'    each window in each input sequence.
#'
#' @details
#' The Wootton-Federhen (Wootton and Federhen, 1993) and Trifonov (Trifonov,
#' 1990) algorithms as well as their faster approximations are well described
#' within Orlov and Potapov (2004). These algorithms score less complex sequences
#' closer to 0, and more complex ones closer to 1. Please note that the
#' 'fast' approximation versions of the two methods are not actually faster
#' within [sequence_complexity()], and so speed should not be a major consideration
#' when choosing which method to use within the `universalmotif` package.
#' The DUST algorithm
#' implementation is described in Morgulis et al. (2006). In this case,
#' less complex sequences score higher, and more complex ones closer
#' to 0.
#'
#' Briefly, the Wootton-Federhen complexity score is a reflection of the
#' numbers of each unique letter found in the window (e.g. for DNA, the
#' more of all four letters can be found in the window the higher the
#' score). An increasing Trifonov score is a relection of the numbers of increasingly
#' larger k-mers (e.g. the count of possible 1-mers, 2-mers, 3-mers, ...,
#' until `trifonov.max.word.size`). Finally, the DUST score approaches 0
#' as the count of unique 3-mers increases. (See the final section in
#' the examples to see how different types of sequence compositions affect
#' the methods.)
#'
#' Please note that the authors of the different methods recommend various
#' window sizes and complexity thresholds. The authors of DUST for example,
#' suggest using a window size of 64 and a threshold of 2 for low complexity.
#' Wootton and Federhen suggest a window size of 40, though show that 10
#' and 20 can be appropriate as well (for amino acid sequences). Keep in
#' mind however that these algorithms were implemented at a time when
#' computers were much slower; perhaps the authors would suggest different
#' window sizes today. One thing to note is that the Wootton-Federhen
#' algorithm has a hard limit due to the need to caculate the product from
#' `1:window.size`. This can end up calculating values which are greater
#' than what a double can hold (e.g. try `prod(1:500)`). Its approximation
#' does not suffer from this though, as it skips calculating the product.
#'
#' In terms of speed, the Wootton-Federhen algorithms are fastest, with DUST
#' being 1-3 times slower and the Trifonov algorithms being several times
#' slower (though the exact amount depends on the max word size).
#'
#' @references
#' Morgulis A, Gertz EM, Schaffer AA, Agarwala R (2006). "A fast and symmetric
#' DUST implementation to mask low-complexity DNA sequences." *Journal of
#' Computational Biology*, **13**, 1028-1040.
#'
#' Orlov YL, Potapov VN (2004). "Complexity: an internet resource for analysis
#' of DNA sequence complexity." *Nucleic Acids Research*, **32**, W628-W633.
#'
#' Trifonov EN (1990). "Making sense of the human genome." In Sarma RH, Sarma
#' MH (Eds), *Structure & Methods* Adenine Press, Albany, **1**, 69-77.
#'
#' Wootton JC, Federhen S (1993). "Statistics of local complexity in amino acid
#' sequences and sequence databases." *Computers & Chemistry*, **17**, 149-163.
#'
#' @examples
#' ## Feel free to play around with different toy sequences to get a feel for
#' ## how the different methods perform
#'
#' library(Biostrings)
#' test.seq <- DNAStringSet(c("AAAAAAAAAAA", "ATGACTGATGC"))
#'
#' sequence_complexity(test.seq, method = "WoottonFederhen")
#' sequence_complexity(test.seq, method = "WoottonFederhenFast")
#' sequence_complexity(test.seq, method = "Trifonov")
#' sequence_complexity(test.seq, method = "TrifonovFast")
#' sequence_complexity(test.seq, method = "DUST")
#'
#' ## You could also use this in conjuction with mask_ranges() to hide
#' ## low complexity regions from scanning, de novo motif discovery, etc
#'
#' if (requireNamespace("GenomicRanges", quiet = TRUE)) {
#' data(ArabidopsisPromoters)
#'
#' # Calculate complexity in 20 bp windows, sliding every 1 bp
#' to.mask <- sequence_complexity(ArabidopsisPromoters, method = "DUST",
#'   window.size = 20, window.overlap = 19, return.granges = TRUE)
#'
#' # Get the ranges with a complexity score greater than 3.5
#' to.mask <- to.mask[to.mask$complexity > 3.5]
#'
#' # See what the low complexity regions look like
#' ArabidopsisPromoters[IRanges::reduce(to.mask)]
#'
#' # Mask them with the default '-' character:
#' mask_ranges(ArabidopsisPromoters, to.mask)
#' }
#'
#' ## To demonstrate how the methods work, consider:
#' ## (These examples use the calc_complexity() utility which utilizes
#' ## the same algorithms and works on character vectors, but lacks
#' ## the ability to use sliding windows.)
#' a <- "ACGT"
#'
#' # For Wootton-Federhen, it can be easily shown it is only dependent
#' # on the counts of individual letters (though do note that the
#' # original paper discusses this method in the context of amino acid
#" # sequences and not DNA):
#' calc_complexity("AAACCCGGGTTT", alph = a)  # 0.7707
#' calc_complexity("AACCGGTTACGT", alph = a)  # 0.7707
#' calc_complexity("ACGTACGTACGT", alph = a)  # 0.7707
#'
#' # As letters start to see drops in counts, the scores go down too:
#' calc_complexity("AAAACCCCGGGG", alph = a)  # 0.6284
#' calc_complexity("AAAAAACCCCCC", alph = a)  # 0.4105
#' calc_complexity("AAAAAAAAAACC", alph = a)  # 0.2518
#'
#' # Trifonov on the other hand is greatly affected by the number
#' # of higher order combinations:
#' calc_complexity("AAACCCGGGTTT", c = "Trifonov", alph = a)  # 0.6364
#' calc_complexity("AACCGGTTACGT", c = "Trifonov", alph = a)  # 0.7273
#'
#' # This next one may seem surprising, but it indeed scores very low.
#' # This is because although it has many of each individual letter,
#' # the number of higher order letter combinations in fact is quite
#' # low due to this particular repeating pattern!
#' calc_complexity("ACGTACGTACGT", c = "Trifonov", alph = a)  # 0.01231
#'
#' # By extension, this means it scores sequences with fewer
#' # counts of individual letters lower too.
#' calc_complexity("AAAACCCCGGGG", c = "Trifonov", alph = a)  # 0.2386
#' calc_complexity("AAAAAACCCCCC", c = "Trifonov", alph = a)  # 0.0227
#' calc_complexity("AAAAAAAAAACC", c = "Trifonov", alph = a)  # 0.0011
#'
#' # As for DUST, it considers the number of 3-mers in the sequence.
#' # The higher the numbers of 3-mers, the lower the score.
#' # (0 = the max possible number of DNA 3-mers for the window size)
#' calc_complexity("AAACCCGGGTTT", c = "DUST", alph = a)  # 0
#' calc_complexity("AACCGGTTACGT", c = "DUST", alph = a)  # 0
#' calc_complexity("ACGTACGTACGT", c = "DUST", alph = a)  # 0.8889
#' calc_complexity("AAAACCCCGGGG", c = "DUST", alph = a)  # 0.333
#' calc_complexity("ACGACGACGACG", c = "DUST", alph = a)  # 1.333
#' calc_complexity("AAAAAACCCCCC", c = "DUST", alph = a)  # 1.333
#' # Similarly to Trifonov, the next one also scores as less complex
#' # compared to the previous one:
#' calc_complexity("ACACACACACAC", c = "DUST", alph = a)  # 2.222
#' calc_complexity("AAAAAAAAAACC", c = "DUST", alph = a)  # 3.111
#' calc_complexity("AAAAAAAAAAAC", c = "DUST", alph = a)  # 4
#' calc_complexity("AAAAAAAAAAAA", c = "DUST", alph = a)  # 5
#'
#' # Just to show once more why the seemingly more complex sequences
#' # such as "ACACACACACAC" score as less complex than "AAAAAACCCCCC"
#' # for the Trifonov and DUST methods:
#' count_klets("ACACACACACAC", k = 3)  # Only 2 possible 3-mers
#' count_klets("AAAAAACCCCCC", k = 3)  # Now 4 possible 3-mers!
#' 
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @seealso [calc_complexity()], [count_klets()], [get_bkg()], [mask_ranges()],
#'  [mask_seqs()]
#' @export
sequence_complexity <- function(seqs, window.size = 20,
  window.overlap = round(window.size / 2),
  method = c("WoottonFederhen", "WoottonFederhenFast", "Trifonov", "TrifonovFast", "DUST"),
  trifonov.max.word.size = 7, nthreads = 1, return.granges = FALSE) {

  method <- match.arg(method)

  if (!is(seqs, "XStringSet")) {
    stop(wmsg("`seqs` should be an `XStringSet` object"), call. = FALSE)
  }

  if (method == "DUST") {
    if (!seqtype(seqs) %in% c("DNA", "RNA")) {
      stop(wmsg("If `method = \"DUST\"`, then `seqs` must be DNA/RNA"),
        call. = FALSE)
    }
    seqs <- DNAStringSet(seqs)
  }

  # param check --------------------------------------------
  args <- as.list(environment())
  all_checks <- character(0)
  num_check <- check_fun_params(list(window.size = args$window.size,
                                     window.overlap = args$window.overlap,
                                     trifonov.max.word.size = args$trifonov.max.word.size,
                                     nthreads = args$nthreads),
                                c(0, 0, 1, 1), FALSE, TYPE_NUM)
  logi_check <- check_fun_params(list(return.granges = args$return.granges),
                                 1, FALSE, TYPE_LOGI)
  all_checks <- c(all_checks, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  seq.names <- names(seqs)
  if (is.null(seq.names)) seq.names <- as.character(seq_len(length(seqs)))

  window.size <- rep_len(window.size, length(seqs))
  window.overlap <- rep_len(window.overlap, length(seqs))

  if (any(window.size > width(seqs))) {
    # warning(wmsg("Found sequences smaller than `window.size`"), call. = FALSE,
    #   immediate. = FALSE)
    window.overlap[window.size > width(seqs)] <- 0
    window.size[window.size > width(seqs)] <- width(seqs)[window.size > width(seqs)]
  }

  if (any(window.size == 0))
    stop("`window.size` must be greater than 0")

  seqlens <- width(seqs)

  window.size[window.size < 1] <- as.integer(window.size[window.size < 1] *
    seqlens[window.size < 1])
  window.overlap[window.overlap < 1 & window.overlap > 0] <- as.integer(
    window.overlap[window.overlap < 1 & window.overlap > 0] *
    seqlens[window.overlap < 1 & window.overlap > 0])

  if (any(window.size <= window.overlap)) {
    stop("`window.overlap` cannot be larger than or equal to `window.size`")
  }

  wins <- mapply(calc_wins, seqlens, window.size, window.overlap, 
    SIMPLIFY = FALSE)
  starts <- lapply(wins, function(x) x$starts)
  stops <- lapply(wins, function(x) x$stops)

  seqs.c <- as.character(seqs)

  alph <- vapply(seqs.c, get_alphabet_cpp, character(1))
  alph <- get_alphabet_cpp(collapse_cpp(alph))

  alph <- switch(seqtype(seqs),
    DNA = get_alphabet_cpp(paste0("ACGT", alph)),
    RNA = get_alphabet_cpp(paste0("ACGU", alph)),
    AA  = get_alphabet_cpp(paste0(AA_STANDARD2, alph, collapse = "")),
    alph
  )

  if (nchar(alph) <= 1) {
    stop(wmsg("Can't calculate complexity when alphabet is a single letter [",
        alph, "]"), call. = FALSE)
  }

  complx <- mapply(sliding_complexity_cpp, seqs.c, window.size, window.overlap,
    MoreArgs = list(metric = method, alph = alph, maxWordSize = trifonov.max.word.size,
      nthreads = nthreads), SIMPLIFY = FALSE)

  res <- mapply(build_complx_res, complx, starts, stops, seq.names,
    SIMPLIFY = FALSE)

  res <- as(do.call(rbind, res), "DataFrame")
  rownames(res) <- NULL

  if (return.granges) {
    colnames(res)[1] <- "seqname"
    colnames(res)[3] <- "end"
    res <- granges_fun(GenomicRanges::GRanges(res,
        seqlengths = structure(width(seqs), names = seq.names)))
  }

  res

}

build_complx_res <- function(complexities, starts, stops, seqname) {
  data.frame(
    sequence = seqname, start = as.integer(starts),
    stop = as.integer(stops), complexity = complexities
  )
}
