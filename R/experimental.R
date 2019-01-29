# UNIVERSALMOTIF CLUSTERS

# Idea:
# Have an object with multiple motifs, and some additional meta info. This meta
# info includes total max/min width, max/min inter-motif width, motif order.
# The idea is to allow for scanning and enriching for clusters of motifs.

universalcluster <- setClass("universalcluster",
                             slots = list(motifs = "list",
                                          meta = "character"))

setValidity("universalcluster",
            function(object) {

              msg <- vector()
              valid <- TRUE

              motifs <- object@motifs
              motif.types <- vapply(motifs, class, character(1))
              if (!all(motif.types == "universalcluster")) {
                msg <- c(msg, "All motifs must be 'universalmotif' objects")
                valid <- FALSE
              }

            })

setMethod("initialize", signature = "universalcluster",
          definition = function(.Object, motifs, meta) {

            .Object@motifs <- motifs
            .Object@meta <- meta

            .Object

          })

setMethod("show", signature = "universalcluster",
          definition = function(object) {

            cat("universalcluster object with", length(object@motifs),
                "motifs\n")
            cat("Meta info:\n")
            print(object@meta)

          })

create_cluster <- function(motifs, meta) {

  motifs <- convert_motifs(motifs)
  motifs <- convert_type(motifs, "PWM")

  cluster <- universalcluster(motifs, meta)

  cluster

}

# Finding binding hotspots in sequences:
data(ArabidopsisMotif)
data(ArabidopsisPromoters)
res <- scan_sequences(ArabidopsisMotif, ArabidopsisPromoters)
# get kernel density (from sequences of length 1000)
kdens <- KernSmooth::bkde(res$start, range.x = c(0, 1000), gridsize = 1000)
# get local maxima
lmax <- zoo::rollapply(kdens$y, 3, function(x) which.max(x) == 2,
                       align = "center")
# plot
plot(kdens, col = as.integer(lmax) + 1)

plot_peaks <- function(res, seq.length, gridsize = seq.length,
                       bandwidth = "auto", maxima.window = 3,
                       use.ggplot = TRUE) {

  if (bandwidth == "auto") {
    kdens <- KernSmooth::bkde(res$start, range.x = c(0, seq.length),
                              gridsize = gridsize)
  } else {
    kdens <- KernSmooth::bkde(res$start, range.x = c(0, seq.length),
                              gridsize = gridsize, bandwidth = bandwidth)
  }

  if (!use.ggplot) {

    lmax <- zoo::rollapply(kdens$y, maxima.window,
                           function(x) which.max(x) == (maxima.window + 1) / 2,
                           align = "center")
    plot(kdens, col = as.integer(lmax) + 1, type = "l")
    points(x = kdens$x[lmax], y = kdens$y[lmax], col = "red", pch = 21, bg = "red")

  } else {

    kdens.df <- data.frame(x = kdens$x, y = kdens$y)
    ggplot2::ggplot(kdens.df, aes(x, y)) +
      ggplot2::geom_line() +
      ggpmisc::stat_peaks(colour = "red", span = maxima.window) +
      theme_bw()

  }

}

# Calculating peak significance: need to get probability of having a peak at
# a certain height. Get null distribution from peak heights  in random data to
# find peak height above 95th percentile. Any peaks in real data above that
# height are then significant. Might need to keep bandwidth param constant.
# Solution: use bkg to calculate expected number of hits in sequences,
# then generate random hit locations for getting random peaks.

# example: AAAAA (probA = 0.3) in 1000 length sequences
# bkg: A=0.3477, C=0.1616, G=0.1517, T=0.3390)
# 0.3477^5 * 1000 = 5.082
# Then repeat sample(1000, sample(5:6, 1, prob = c(.57, .43)) for the number
# of sequences in set. Get distribution of peak heights.

bkgA.freqs <- Biostrings::oligonucleotideFrequency(ArabidopsisPromoters, 1)[, 1] / 1000
hit.counts <- bkgA.freqs^5 * 1000
random.hits <- do.call(c, sapply(hit.counts, function(x) sample(1:1000, x)))
random.dens <- KernSmooth::bkde(random.hits, range.x = c(1, 1000), gridsize = 1000)
peak.heights <- random.dens$y[find_peaks(random.dens$y)]

# Now repeat this like 1000 times or so to get a null distribution.

get_heights <- function(res.hits, seq.len, bandwidth = NULL) {

  if (is.null(bandwidth)) 
    x <- KernSmooth::bkde(res.hits, range.x = c(1, seq.len), gridsize = seq.len)
  else
    x <- KernSmooth::bkde(res.hits, range.x = c(1, seq.len), gridsize = seq.len,
                          bandwidth = bandwidth)
  x$y[ggpmisc:::find_peaks(x$y)]

}

get_heights_ran <- function(hit.counts, seq.len, repeats = 1000, bandwidth = NULL) {

  heights <- sapply(seq_len(repeats),
                    function(x) {

                      x1 <- do.call(c, sapply(hit.counts,
                                              function(x) sample(1:seq.len, x)))
                      if (is.null(bandwidth))
                        x2 <- KernSmooth::bkde(x1, range.x = c(1, seq.len),
                                               gridsize = seq.len)
                      else
                        x2 <- KernSmooth::bkde(x1, range.x = c(1, seq.len),
                                               gridsize = seq.len,
                                               bandwidth = bandwidth)
                      x2$y[ggpmisc:::find_peaks(x2$y)]

                    })

  do.call(c, heights)

}

######
# Complete workthrough

bandwidth.all <- 20

random.heights <- get_heights_ran(hit.counts, 1000, bandwidth = bandwidth.all)
AtProm.heights <- get_heights(res2$start, 1000, bandwidth = bandwidth.all)

1 - ecdf(random.heights)(AtProm.heights)

percentiles <- quantile(random.heights, c(0.95, 0.99, 0.999, 0.999999))

plot_peaks(res2, 1000, bandwidth = bandwidth.all, use.ggplot = F)
abline(h = percentiles[4], col = "blue")
