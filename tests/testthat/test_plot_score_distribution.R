context("plot_score_distribution()")

small_motif <- function() {
  create_motif(Biostrings::DNAStringSet(c("TCAG", "TCAG", "TGAG", "TCTG")),
    pseudocount = 1, name = "small")
}

## Index of the first layer drawn with the given geom.
geom_layer <- function(p, geom) {
  which(vapply(p$layers, function(l) class(l$geom)[1], character(1)) == geom)
}

## The bin edges and probabilities of the filled distribution in one panel.
dist_bars <- function(p, panel = 1) {
  d <- ggplot2::layer_data(p, geom_layer(p, "GeomRect"))
  d <- d[d$PANEL == panel, ]
  list(edges = c(d$xmin, d$xmax[nrow(d)]), mass = d$ymax)
}

## The probability of each bin from scoring every word, weighted by its
## background probability, with the PWM put on the same 0.001 grid as
## motif_pvalue().
brute_force <- function(motif, edges) {
  motif <- convert_type(motif, "PWM")
  pwm <- trunc(motif@motif * 1000)
  bkg <- motif@bkg[rownames(pwm)]
  idx <- as.matrix(expand.grid(rep(list(seq_len(nrow(pwm))), ncol(pwm))))
  scores <- rowSums(matrix(pwm[cbind(as.vector(idx), rep(seq_len(ncol(pwm)),
    each = nrow(idx)))], nrow = nrow(idx)))
  probs <- apply(matrix(bkg[idx], nrow = nrow(idx)), 1, prod)
  grid <- round(edges * 1000)
  bin <- findInterval(scores, grid)
  expect_true(min(bin) >= 1 && max(bin) < length(grid))
  vapply(seq_len(length(grid) - 1), function(i) sum(probs[bin == i]),
    numeric(1))
}

test_that("probabilities sum to one and match enumeration", {
  m <- small_motif()
  for (bins in c(7, 30, 1e5)) {
    p <- plot_score_distribution(m, pvalues = NULL, bins = bins)
    expect_s3_class(p, "ggplot")
    d <- dist_bars(p)
    expect_true(length(d$mass) <= bins)
    expect_equal(sum(d$mass), 1, tolerance = 1e-12)
    expect_equal(d$mass, brute_force(m, d$edges), tolerance = 1e-12)
    r <- motif_range(m)
    expect_equal(d$edges[1], unname(r["min"]))
    expect_equal(d$edges[length(d$edges)], unname(r["max"]) + 0.001)
  }
})

test_that("one bin per score gives the enumerated distribution", {
  m <- small_motif()
  d <- dist_bars(plot_score_distribution(m, pvalues = NULL, bins = 1e5))
  scores <- round(get_scores(m) * 1000)
  want <- table(scores) / length(scores)
  grid <- round(d$edges[-length(d$edges)] * 1000)
  got <- d$mass[d$mass > 0]
  expect_equal(grid[d$mass > 0], as.numeric(names(want)))
  expect_equal(got, as.numeric(want), tolerance = 1e-12)
})

test_that("an uneven background matches enumeration", {
  m <- small_motif()
  m["bkg"] <- c(A = 0.1, C = 0.4, G = 0.3, T = 0.2)
  d <- dist_bars(plot_score_distribution(m, pvalues = NULL, bins = 25))
  expect_equal(sum(d$mass), 1, tolerance = 1e-12)
  expect_equal(d$mass, brute_force(m, d$edges), tolerance = 1e-12)
})

test_that("bkg.probs scores against that background", {
  data(ArabidopsisMotif)
  m <- ArabidopsisMotif
  bkg <- c(0.3, 0.2, 0.2, 0.3)
  p <- plot_score_distribution(m, bkg.probs = bkg, pvalues = 1e-4, bins = 60)
  d <- dist_bars(p)
  expect_equal(sum(d$mass), 1, tolerance = 1e-12)
  tail <- motif_pvalue(m, score = d$edges[-c(1, length(d$edges))],
    bkg.probs = bkg)
  expect_equal(d$mass, -diff(c(1, tail, 0)), tolerance = 1e-15)
  expect_true(d$mass[length(d$mass)] > 0)
  cut <- ggplot2::layer_data(p, geom_layer(p, "GeomVline"))
  expect_equal(cut$xintercept, motif_pvalue(m, pvalue = 1e-4, bkg.probs = bkg))
  ## A named background, and one per motif, give the same picture.
  named <- c(A = 0.3, C = 0.2, G = 0.2, T = 0.3)
  p2 <- plot_score_distribution(list(m, m), bkg.probs = list(named, bkg),
    pvalues = 1e-4, bins = 60)
  expect_equal(dist_bars(p2, 1)$mass, d$mass)
  expect_equal(dist_bars(p2, 2)$mass, d$mass)
})

test_that("a cutoff line sits at the threshold and bounds the tail", {
  data(ArabidopsisMotif)
  m <- ArabidopsisMotif
  p <- plot_score_distribution(m, pvalues = 1e-4, bins = 1e5)
  threshold <- motif_pvalue(m, pvalue = 1e-4)
  cut <- ggplot2::layer_data(p, geom_layer(p, "GeomVline"))
  expect_equal(cut$xintercept, threshold)
  d <- dist_bars(p)
  grid <- round(d$edges[-length(d$edges)] * 1000)
  at <- match(round(threshold * 1000), grid)
  expect_false(is.na(at))
  n <- length(d$mass)
  expect_true(sum(d$mass[at:n]) <= 1e-4)
  expect_equal(sum(d$mass[at:n]), motif_pvalue(m, score = threshold),
    tolerance = 1e-9)
  ## The threshold is the lowest such score.
  lower <- max(which(d$mass[seq_len(at - 1)] > 0))
  expect_true(sum(d$mass[lower:n]) > 1e-4)
  lab <- ggplot2::layer_data(p, geom_layer(p, "GeomLabel"))
  hi <- motif_range(m)["max"]
  expect_equal(lab$label, sprintf("P = 1e-4, score %.2f (%.0f%%)",
    threshold, 100 * threshold / hi))
})

test_that("several cutoffs are drawn", {
  data(ArabidopsisMotif)
  ps <- c(1e-3, 1e-5, 1e-4)
  p <- plot_score_distribution(ArabidopsisMotif, pvalues = ps)
  cut <- ggplot2::layer_data(p, geom_layer(p, "GeomVline"))
  expect_equal(sort(cut$xintercept),
    sort(motif_pvalue(ArabidopsisMotif, pvalue = ps)))
  lab <- ggplot2::layer_data(p, geom_layer(p, "GeomLabel"))
  expect_equal(sort(sub(",.*", "", lab$label)),
    c("P = 0.001", "P = 1e-4", "P = 1e-5"))
  p <- plot_score_distribution(ArabidopsisMotif, pvalues = NULL)
  expect_length(geom_layer(p, "GeomVline"), 0)
})

test_that("an unreachable cutoff warns and is skipped", {
  m <- create_motif("TGACT", pseudocount = 1, name = "short")
  expect_warning(p <- plot_score_distribution(m, pvalues = c(1e-2, 1e-9)),
    "cannot reach a P-value of 1e-9")
  cut <- ggplot2::layer_data(p, geom_layer(p, "GeomVline"))
  expect_equal(cut$xintercept, motif_pvalue(m, pvalue = 1e-2))
})

test_that("hits come from scan tables, lists and vectors", {
  data(ArabidopsisMotif)
  data(ArabidopsisPromoters)
  m1 <- ArabidopsisMotif
  m2 <- create_motif("TATAWAWR", name = "TATA", pseudocount = 1)
  hits <- scan_sequences(list(m1, m2), ArabidopsisPromoters,
    threshold = 1e-3, threshold.type = "pvalue", nthreads = 1)
  p <- plot_score_distribution(list(m1, m2), hits = hits, pvalues = 1e-3)
  rug <- ggplot2::layer_data(p, geom_layer(p, "GeomRug"))
  key <- ggplot2::layer_data(p, geom_layer(p, "GeomText"))
  for (i in 1:2) {
    want <- sort(hits$score[hits$motif == c(m1@name, m2@name)[i]])
    expect_true(length(want) > 0)
    expect_equal(sort(rug$x[rug$PANEL == i]), want)
    expect_equal(key$label[key$PANEL == i],
      paste0("hits (n = ", length(want), ")"))
  }

  lite <- scan_sequences_lite(m1, ArabidopsisPromoters, pvalue = 1e-3,
    nthreads = 1)
  p <- plot_score_distribution(m1, hits = lite)
  rug <- ggplot2::layer_data(p, geom_layer(p, "GeomRug"))
  expect_equal(sort(rug$x), sort(lite$score))

  scores <- c(2.5, 9, 12.25)
  p <- plot_score_distribution(m1, hits = scores)
  expect_equal(ggplot2::layer_data(p, geom_layer(p, "GeomRug"))$x, scores)
  p <- plot_score_distribution(list(m1, m2), hits = list(TATA = scores),
    pvalues = 1e-3)
  rug <- ggplot2::layer_data(p, geom_layer(p, "GeomRug"))
  expect_equal(as.integer(rug$PANEL), rep(2L, 3))
  expect_equal(rug$x, scores)
})

test_that("log.y draws the tail on a log scale", {
  data(ArabidopsisMotif)
  m <- ArabidopsisMotif
  p <- plot_score_distribution(m, log.y = TRUE)
  expect_equal(p$scales$get_scales("y")$trans$name, "log-10")
  expect_equal(p$labels$y, "Probability (log10)")
  lin <- dist_bars(plot_score_distribution(m))$mass
  d <- ggplot2::layer_data(p, geom_layer(p, "GeomRect"))
  pos <- lin > 0
  expect_equal(10^d$ymax[pos], lin[pos])
  ## The floor sits below the smallest probability, and the top bin, holding
  ## the best words, is drawn.
  expect_true(all(10^d$ymin < min(lin[pos])))
  expect_true(lin[length(lin)] >= motif_pvalue(m, score = motif_range(m)["max"]))
  expect_equal(p$labels$y, "Probability (log10)")
})

test_that("panels, titles and the percentage axis", {
  data(ArabidopsisMotif)
  m <- ArabidopsisMotif
  motifs <- list(m, create_motif("CACGTG", name = "G-box", pseudocount = 1),
    create_motif("TATAWAWR", name = "TATA", pseudocount = 1))
  p <- plot_score_distribution(motifs, ncol = 2, pvalues = 1e-3)
  b <- ggplot2::ggplot_build(p)
  expect_equal(levels(b$layout$layout$motif), c(m@name, "G-box", "TATA"))
  expect_equal(max(b$layout$layout$COL), 2)
  expect_equal(max(b$layout$layout$ROW), 2)
  expect_false(inherits(p$scales$get_scales("x")$secondary.axis,
    "AxisSecondary"))

  p <- plot_score_distribution(m)
  sec <- p$scales$get_scales("x")$secondary.axis
  expect_true(inherits(sec, "AxisSecondary"))
  expect_equal(sec$name, "% of maximum score")
  hi <- unname(motif_range(m)["max"])
  expect_equal(sec$trans(c(-10, hi)), c(-1000 / hi, 100))
  expect_equal(p$labels$x, "Score")
  expect_equal(p$labels$y, "Probability")

  ## Duplicated names get their own panels.
  p <- plot_score_distribution(list(m, m))
  expect_equal(levels(ggplot2::ggplot_build(p)$layout$layout$motif),
    c(m@name, paste0(m@name, ".1")))
})

test_that("bad input is refused", {
  data(ArabidopsisMotif)
  m <- ArabidopsisMotif
  expect_error(plot_score_distribution(list()), "no motifs")
  for (bad in list(0, -1, 1.5, TRUE, c(10, 20)))
    expect_error(plot_score_distribution(m, bins = bad), "bins")
  for (bad in list(0, 1.5, c(1e-4, -1), "a"))
    expect_error(plot_score_distribution(m, pvalues = bad), "pvalues")
  expect_error(plot_score_distribution(list(m, m, m),
    bkg.probs = list(rep(0.25, 4), rep(0.25, 4))), "one per motif")
  expect_error(plot_score_distribution(m,
    hits = data.frame(motif = m@name)), "score")
  expect_error(plot_score_distribution(list(m, m),
    hits = data.frame(score = 1)), "motif")
  expect_error(plot_score_distribution(list(m, m), hits = c(1, 2)),
    "several motifs")
  expect_error(plot_score_distribution(m, hits = list(`G-box` = 1)),
    "not plotted: G-box")
  expect_error(plot_score_distribution(m, log.y = NA), "log.y")
})

test_that("P-values are formatted for labels", {
  expect_equal(format_pvalue_label(c(0.05, 1e-3, 1e-4, 2.5e-6, 9.999e-5)),
    c("0.05", "0.001", "1e-4", "2.5e-6", "1e-4"))
})
