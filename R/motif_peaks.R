#' Look for overrepresented motif position peaks in a set of sequences.
#'
#' Using the motif position data from [scan_sequences()] (or elsewhere),
#' test whether certain positions in the sequences have significantly higher
#' motif density.
#'
#' @param hits `numeric` A vector of sequence positions indicating motif sites.
#' @param seq.length `numeric(1)` Length of sequences. Only one number is
#'    allowed, as all sequences must be of identical length. If missing,
#'    then the largets number from `hits` is used.
#' @param seq.count `numeric(1)` Number of sequences with motif sites. If
#'    missing, then the number of unique `hits` values is used.
#' @param bandwidth `numeric(1)` Peak smoothing parameter. Smaller numbers
#'    will result in skinnier peaks, larger numbers will result in wider
#'    peaks. Leaving this empty will cause [motif_peaks()] to generate one
#'    by itself.
#' @param max.p `numeric(1)` Maximum P-value allowed for finding significant
#'    motif site peaks.
#' @param peak.width `numeric(1)` Minimum peak width. A peak is defined as
#'    as the highest point within the value set by `peak.width`.
#' @param nrand `numeric(1)` Number of random permutations for generating a
#'    null distribution. In order to calculate P-values, a set of random
#'    motif site positions are generated `nrand` times.
#' @param plot `logical(1)` Will create a `ggplot2` object displaying motif
#'    peaks.
#' @param BP `logical(1)` Allows for the use of \pkg{BiocParallel} within
#'    [motif_peaks()]. See [BiocParallel::register()] to change the
#'    default backend. Setting `BP = TRUE` is only recommended for
#'    exceptionally large jobs. Keep in mind that this function will not
#'    attempt to limit its memory usage.
#'
#' @details
#'    Kernel smoothing is used to calculate motif position density. The
#'    implementation for this process is based on code from the
#'    \pkg{KernSmooth} R package. These density estimates are used to
#'    determine peak locations and heights. To calculate the P-values of
#'    these peaks, a null distribution is calculated from peak heights of
#'    randomly generated motif positions.
#'
#'    If the `bandwidth` option is not supplied, then the following code is used
#'    (from \pkg{KernSmooth}):
#'
#'    `del0 <- (1 / (4 * pi))^(1 / 10)`
#'
#'    `bandwidth <- del0 * (243 / (35 * length(hits)))^(1 / 5) * sqrt(var(hits))`
#'
#' @return A `data.frame` with peak positions and P-values. If `plot = TRUE`,
#'    then a list is returned with the `data.frame` as the first item and
#'    the `ggplot2` object as the second item.
#'
#' @examples
#' data(ArabidopsisMotif)
#' data(ArabidopsisPromoters)
#' hits <- scan_sequences(ArabidopsisMotif, ArabidopsisPromoters, RC = FALSE,
#'                        verbose = 0, progress = FALSE, threshold = 0,
#'                        threshold.type = "logodds")
#' res <- motif_peaks(hits$start, 1000, 50)
#' # Open plot:
#' res$Plot
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [scan_sequences()]
#' @export
motif_peaks <- function(hits, seq.length, seq.count, bandwidth, max.p = 10^-6,
                        peak.width = 3, nrand = 1000, plot = TRUE, BP = FALSE) {

# TODO: add tests and vignette section

# Plans:
#  - Compare peak locations to a set of bkg peaks
#  - Look for peak co-occurence between multiple motifs

  # param check --------------------------------------------
  args <- as.list(environment())
  num_check <- check_fun_params(list(hits = args$hits,
                                     seq.length = args$seq.length,
                                     seq.count = args$seq.count,
                                     bandwidth = args$bandwidth,
                                     max.p = args$max.p,
                                     peak.width = args$peak.width,
                                     nrand = args$nrand),
                                c(0, rep(1, 6)),
                                c(FALSE, TRUE, TRUE, TRUE,
                                  FALSE, FALSE, FALSE), "numeric")
  logi_check <- check_fun_params(list(plot = args$plot, BP = args$BP),
                                 numeric(), logical(), "logical")
  all_checks <- c(num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  if (missing(seq.length)) seq.length <- max(hits)
  if (missing(seq.count)) seq.count <- length(unique(hits))
  if (missing(bandwidth)) {
    del0 <- (1 / (4 * pi))^(1 / 10)
    bandwidth <- del0 * (243 / (35 * length(hits)))^(1 / 5) * sqrt(var(hits))
  }

  hit.count.perseq <- length(hits) / seq.count

  rand.hits <- lapply(seq_len(nrand),
                      function(x) sample(seq_len(seq.length), length(hits),
                                         replace = TRUE))
  # this is the slowest step
  rand.kern <- lapply_(rand.hits, function(x) kern_fun(x, bandwidth, seq.length,
                                                     seq.length),
                       BP = BP, PB = FALSE)
  # second slowest
  rand.peaks <- lapply_(rand.kern, function(x) x$y[peakfinder_cpp(x$y, peak.width)],
                        BP = BP, PB = FALSE)
  rand.peaks <- do.call(c, rand.peaks)

  data.kern <- kern_fun(hits, bandwidth, seq.length, seq.length)
  data.loc <- peakfinder_cpp(data.kern$y, peak.width)
  data.peaks <- data.kern$y[data.loc]

  peak.pvals <- 1 - ecdf(rand.peaks)(data.peaks)

  if (plot) {
    pval.lim <- quantile(rand.peaks, 1 - max.p)
    kern.df <- data.frame(x = data.kern$x, y = data.kern$y)
    p <- ggplot(kern.df, aes(x, y)) +
           geom_line() +
           geom_point(data = data.frame(x = data.kern$x[data.loc],
                                        y = data.kern$y[data.loc]),
                      colour = "red") +
           xlab("Sequence location") +
           ylab("Kernel density") +
           geom_hline(yintercept = pval.lim, colour = "blue") +
           geom_text(data = data.frame(x = 0, y = pval.lim), aes(x, y),
                     hjust = 0, vjust = -1,
                     label = paste("Cutoff for P-value =<", max.p)) +
           theme_bw()
  }

  out <- data.frame(Peak = data.kern$x[data.loc], Pval = peak.pvals)
  out <- out[order(out$Pval), ]
  out <- out[out$Pval <= max.p, ]
  if (nrow(out) == 0) {
    message("No significant peaks found.")
    if (plot) return(p) else return(invisible(NULL))
  }
  rownames(out) <- NULL
  if (plot) return(list(Peaks = out, Plot = p)) else return(out)

}

kern_fun <- function(x, bandwidth, gridsize, range.x) {

  # modified from KernSmooth::bkde

  n <- length(x)
  M <- gridsize
  tau <- 4
  h <- bandwidth
  a <- 1
  b <- range.x

  gpoints <- seq(a, b, length = M)
  gcounts <- linbin_cpp(x, gpoints)

  delta <- (b - a) / (h * (M - 1L))
  L <- min(floor(tau / delta), M)

  lvec <- seq(0L, L, 1)
  kappa <- dnorm(lvec * delta) / (n * h)

  P <- 2^(ceiling(log(M + L + 1L) / log(2)))
  kappa <- c(kappa, rep(0, P - 2L * L - 1L), rev(kappa[-1L]))
  tot <- sum(kappa) * (b - a) / (M - 1L) * n
  gcounts <- c(gcounts, rep(0, P - M))
  kappa <- fft(kappa / tot)
  gcounts <- fft(gcounts)

  list(x = gpoints, y = (Re(fft(kappa * gcounts, TRUE)) / P)[1L:M])

}
