#' Look for overrepresented motif position peaks in a set of sequences.
#'
#' Using the motif position data from [scan_sequences()] (or elsewhere),
#' test whether certain positions in the sequences have significantly higher
#' motif density.
#'
#' @param hits `numeric` A vector of sequence positions indicating motif sites.
#' @param seq.length `numeric(1)` Length of sequences. Only one number is
#'    allowed, as all sequences must be of identical length.
#' @param seq.count `numeric(1)` Number of sequences with motif sites.
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
                                c(FALSE, FALSE, FALSE, TRUE,
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
  rand.kern <- lapply_(rand.hits, function(x) my_kern(x, bandwidth, seq.length,
                                                     c(1, seq.length)),
                       BP = BP, PB = FALSE)
  # second slowest
  rand.peaks <- lapply_(rand.kern, function(x) x$y[my_peakfinder(x$y, peak.width)],
                        BP = BP, PB = FALSE)
  rand.peaks <- do.call(c, rand.peaks)

  data.kern <- my_kern(hits, bandwidth, seq.length, c(1, seq.length))
  data.loc <- my_peakfinder(data.kern$y, peak.width)
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

my_kern <- function(x, bandwidth, gridsize, range.x) {

  # modified from KernSmooth::bkde

  n <- length(x)
  M <- gridsize

  tau <- 4
  h <- bandwidth

  if (missing(range.x)) range.x <- c(min(x) - tau * h, max(x) + tau * h)

  a <- range.x[1L]
  b <- range.x[2L]

  gpoints <- seq(a, b, length = M)
  gcounts <- my_linbin(x, gpoints)

  delta <- (b - a) / (h * (M - 1L))
  L <- min(floor(tau / delta), M)

  lvec <- 0L:L
  kappa <- dnorm(lvec * delta) / (n * h)

  P <- 2^(ceiling(log(M + L + 1L) / log(2)))
  kappa <- c(kappa, rep(0, P - 2L * L - 1L), rev(kappa[-1L]))
  tot <- sum(kappa) * (b - a) / (M - 1L) * n
  gcounts <- c(gcounts, rep(0, P - M))
  kappa <- fft(kappa / tot)
  gcounts <- fft(gcounts)

  list(x = gpoints, y = (Re(fft(kappa * gcounts, TRUE)) / P)[1L:M])

}

my_peakfinder <- function(x, m = 3) {

  # modified from ggpmisc:::find_peaks

  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- vapply(which(shape < 0), function(i) peak_finder_single(i, x, m),
                numeric(1))
  pks[!is.na(pks)]

}

peak_finder_single <- function(i, x, m) {

  z <- i - m + 1
  z <- ifelse(z > 0, z, 1)
  w <- i + m + 1
  w <- ifelse(w < length(x), w, length(x))

  if (all(x[c(z:i, (i + 2):w)] <= x[i + 1]))
    return(i + 1)
  else
    return(as.numeric(NA))

}

my_linbin <- function(x, gpoints) {

  # modified from KernSmooth:::linbin

  M <- length(gpoints)
  gcnts <- rep(0, M)

  delta <- (gpoints[M] - gpoints[1]) / (M - 1)

  for (i in seq_along(x)) {

    lxi <- ((x[i] - gpoints[1]) / delta) + 1
    li <- as.integer(lxi)
    rem <- lxi - li

    if (li >= 1 && li < M) {
      gcnts[li] <- gcnts[li] + (1 - rem)
      gcnts[li + 1] <- gcnts[li + 1] + rem
    }

  }

  gcnts

}
