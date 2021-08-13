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
#' Please note that the authors of the different methods recommend various
#' window sizes and complexity thresholds. The authors of DUST for example,
#' suggest using a window size of 64 and a threshold of 2 for low complexity.
#' Wootton and Federhen suggest a window size of 40, though show that 10
#' and 20 can be appropriate as well.
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
#' # Mask them with the '-' character:
#' mask_ranges(ArabidopsisPromoters, to.mask)
#' }
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
