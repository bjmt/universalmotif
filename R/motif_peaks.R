# Plans:
#  - Compare peak locations to a set of bkg peaks
#  - Look for peak co-occurence between multiple motifs

motif_peaks <- function(hits, seq.length, seq.count, bandwidth, max.p = 10^-6,
                        peak.width = 3, nrand = 1000, plot = TRUE, BP = FALSE) {

  if (missing(bandwidth)) {
    del0 <- (1 / (4 * pi))^(1 / 10)
    bandwidth <- del0 * (243 / (35 * length(hits)))^(1 / 5) * sqrt(var(hits))
  }

  hit.count.perseq <- length(hits) / seq.count

  rand.hits <- lapply(seq_len(nrand),
                      function(x) sample(seq_len(seq.length), length(hits),
                                         replace = TRUE))
  # this is slowest step
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
    # plot(data.kern, type = "l")
    # points(x = data.kern$x[data.loc], y = data.kern$y[data.loc],
           # col = "red", pch = 21, bg = "red")
    # abline(h = pval.lim, col = "blue")
    kern.df <- data.frame(x = data.kern$x, y = data.kern$y)
    p <- ggplot(kern.df, aes(x, y)) +
           geom_line() +
           # stat_peaks(span = peak.width, colour = "red") +  # ggpmisc::stat_peaks
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
    print(p)
  }

  out <- data.frame(Peak = data.kern$x[data.loc], Pval = peak.pvals)
  out <- out[order(out$Pval), ]
  out <- out[out$Pval <= max.p, ]
  if (nrow(out) == 0) {
    message("No significant peaks found.")
    return(invisible(NULL))
  }
  rownames(out) <- NULL
  out

}

my_kern <- function(x, bandwidth, gridsize, range.x) {

  # source: KernSmooth::bkde

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

  # source: ggpmisc:::find_peaks

  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i) {
                                    z <- i - m + 1
                                    z <- ifelse(z > 0, z, 1)
                                    w <- i + m + 1
                                    w <- ifelse(w < length(x), w, length(x))
                                    if (all(x[c(z:i, (i + 2):w)] <= x[i + 1]))
                                      return(i + 1)
                                    else
                                      return(numeric(0))
                                  })
  pks <- unlist(pks)
  pks

}

my_linbin <- function(x, gpoints) {

  # source: KernSmooth:::linbin

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
