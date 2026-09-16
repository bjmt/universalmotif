context("background rebasing and integer-model null probabilities")

test_that("explicit backgrounds rebase scores and null weights in scalar and batch calls", {
  ppm <- matrix(rep(c(.5, .25, .125, .125), 3), 4,
                dimnames = list(c("A", "C", "G", "T"), NULL))
  m <- create_motif(ppm, type = "PPM", pseudocount = 0)
  old_pwm <- convert_type(m, "PWM")
  backgrounds <- list(c(A = .1, C = .4, G = .4, T = .1),
                      c(A = .4, C = .1, G = .1, T = .4))
  words <- as.matrix(expand.grid(rep(list(1:4), 3)))
  queries <- c(-1.2345, 0, .5, 2.7)
  expected <- lapply(backgrounds, function(bg) {
    pwm <- log2(ppm / bg)
    scores <- rowSums(vapply(1:3, function(j)
      trunc(pwm[words[, j], j] * 1000), numeric(nrow(words))))
    weights <- apply(words, 1, function(w) prod(bg[w]))
    vapply(queries, function(q) sum(weights[scores >= ceiling(q * 1000)]), numeric(1))
  })
  for (threads in c(1, 2)) {
    for (source in list(m, old_pwm)) {
      actual <- motif_pvalue(list(source, source), score = rep(list(queries), 2),
                            bkg.probs = backgrounds, nthreads = threads)
      expect_equal(actual, expected, tolerance = 1e-12)
      for (i in seq_along(backgrounds))
        expect_equal(motif_pvalue(source, score = queries, bkg.probs = backgrounds[[i]],
                                 nthreads = threads), expected[[i]], tolerance = 1e-12)
    }
  }
})

test_that("both scanners use the rebased motif background rather than sequence composition", {
  bg <- c(A = .1, C = .4, G = .4, T = .1)
  ppm <- matrix(rep(c(.5, .25, .125, .125), 3), 4,
                dimnames = list(names(bg), NULL))
  m <- create_motif(ppm, type = "PPM", bkg = bg, pseudocount = 0)
  pwm <- log2(ppm / bg)
  words <- as.matrix(expand.grid(rep(list(1:4), 3)))
  scores <- rowSums(vapply(1:3, function(j)
    trunc(pwm[words[, j], j] * 1000), numeric(nrow(words))))
  weights <- apply(words, 1, function(w) prod(bg[w]))
  tails <- vapply(scores, function(s) sum(weights[scores >= s]), numeric(1))
  text <- apply(words, 1, function(w) paste0(names(bg)[w], collapse = ""))
  # Every word occurs once: sequence composition is uniform, unlike bg.
  seqs <- Biostrings::DNAStringSet(text)
  cutoff <- .03
  wanted <- which(tails <= cutoff)
  for (threads in c(1, 2)) {
    results <- list(
      suppressMessages(scan_sequences(m, seqs, threshold = cutoff,
        RC = FALSE, calc.pvals = TRUE, calc.qvals = FALSE, verbose = 0,
        nthreads = threads)),
      scan_sequences_lite(m, seqs, pvalue = cutoff, RC = FALSE,
                          return.granges = FALSE, nthreads = threads))
    for (hits in results) {
      expect_setequal(hits$sequence.i, wanted)
      expect_equal(unname(hits$score), unname(scores[hits$sequence.i]) / 1000)
      expect_equal(unname(hits$pvalue), unname(tails[hits$sequence.i]), tolerance = 1e-12)
      expect_true(all(hits$pvalue <= cutoff))
    }
  }
})
