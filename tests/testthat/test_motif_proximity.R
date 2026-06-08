context("motif_proximity() / plot_motif_proximity()")

PROX_COLS <- c("motif", "motif.i", "method", "count", "orientation",
               "nhits", "best.dist", "hits.near", "hits.far", "frac.near",
               "expected.near", "enrichment", "log2.enrichment", "pvalue",
               "qvalue", "ks.D", "dist.to.nearest")

## A 3-"chromosome" synthetic genome of unequal lengths.
make_genome <- function() {
  c1 <- create_sequences(seqnum = 1, seqlen = 6000, rng.seed = 11)
  c2 <- create_sequences(seqnum = 1, seqlen = 4000, rng.seed = 12)
  c3 <- create_sequences(seqnum = 1, seqlen = 3000, rng.seed = 13)
  g <- c(c1, c2, c3)
  names(g) <- c("chr1", "chr2", "chr3")
  g
}

## Plant `motif_str` at `offset` from each position on chromosome `chr`.
plant_at <- function(genome, chr, positions, motif_str = "TTGACATA",
                     offset = 5L) {
  ds <- Biostrings::DNAString(motif_str)
  w  <- nchar(motif_str)
  for (p in positions) {
    s <- p + offset
    s <- max(1L, min(length(genome[[chr]]) - w + 1L, s))
    Biostrings::subseq(genome[[chr]], start = s, width = w) <- ds
  }
  genome
}

skip_no_granges <- function() {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomeInfoDb")
}

test_that("output shape and column types are stable", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  r <- suppressMessages(motif_proximity(hits, anchors, qvalue = 1))
  expect_true(is.data.frame(r))
  expect_equal(names(r), PROX_COLS)
  expect_type(r$motif, "character")
  expect_type(r$nhits, "integer")
  expect_type(r$best.dist, "integer")
  expect_type(r$pvalue, "double")
  expect_type(r$qvalue, "double")
  expect_true(is.list(r$dist.to.nearest))
})

test_that("binomial detects a motif planted near anchors", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  r <- suppressMessages(motif_proximity(hits, anchors, qvalue = 1))
  expect_equal(nrow(r), 1L)
  expect_lt(r$pvalue, 1e-4)
  expect_gt(r$enrichment, 2)
})

test_that("ks detects, and ks is rejected for anchors / directional", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  r <- suppressMessages(motif_proximity(hits, anchors, method = "ks",
                                        qvalue = 1))
  expect_equal(nrow(r), 1L)
  expect_lt(r$pvalue, 1e-3)
  expect_gt(r$ks.D, 0)

  expect_error(motif_proximity(hits, anchors, method = "ks", count = "anchors"),
               regexp = "anchors")
  expect_error(motif_proximity(hits, anchors, method = "ks",
                               orientation = "downstream"),
               regexp = "directional")
})

test_that("permutation is reproducible after set.seed() (serial)", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  set.seed(1)
  r1 <- suppressMessages(motif_proximity(hits, anchors, method = "permutation",
                                         n.perm = 99, qvalue = 1))
  set.seed(1)
  r2 <- suppressMessages(motif_proximity(hits, anchors, method = "permutation",
                                         n.perm = 99, qvalue = 1))
  expect_equal(r1, r2)
  expect_lte(r1$pvalue, 0.05)
})

test_that("binomial/ks leave the global RNG untouched; permutation uses it", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  ## A non-random method must not advance the stream.
  set.seed(123); before <- runif(1)
  set.seed(123)
  invisible(suppressMessages(motif_proximity(hits, anchors, qvalue = 1)))
  after <- runif(1)
  expect_equal(before, after)

  ## Permutation does advance it (it draws from the global stream).
  set.seed(123)
  invisible(suppressMessages(
    motif_proximity(hits, anchors, method = "permutation", n.perm = 20,
                    qvalue = 1)))
  after_perm <- runif(1)
  expect_false(isTRUE(all.equal(before, after_perm)))
})

test_that("negative control: random hits are not enriched near anchors", {
  skip_no_granges()
  g <- make_genome()
  ## Plant the motif at random positions unrelated to the anchors.
  set.seed(99)
  rand.pos <- sample.int(5900, 25)
  g <- plant_at(g, "chr1", rand.pos, offset = 0L)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(seq(300, 5700, by = 300), width = 1),
                strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  r <- suppressMessages(motif_proximity(hits, anchors, qvalue = 0.1))
  expect_equal(nrow(r), 0L)
  expect_equal(names(r), PROX_COLS)
})

test_that("count = 'anchors' runs and counts anchors", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  set.seed(1)
  r <- suppressMessages(motif_proximity(hits, anchors, count = "anchors",
                                        method = "permutation", n.perm = 99,
                                        qvalue = 1))
  expect_equal(nrow(r), 1L)
  expect_equal(r$nhits, length(anchors))   # anchors mode counts anchors
  expect_lte(r$hits.near, length(anchors))
})

test_that("internal scan path equals the precomputed-hits path", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  r1 <- suppressMessages(motif_proximity(hits, anchors, qvalue = 1))
  r2 <- suppressMessages(motif_proximity(motifs = m, sequences = g,
                                         anchors = anchors, scan.pvalue = 1e-3,
                                         qvalue = 1))
  expect_equal(r1$pvalue, r2$pvalue)
  expect_equal(r1$enrichment, r2$enrichment)

  expect_error(motif_proximity(hits, anchors, motifs = m, sequences = g),
               regexp = "either")
  expect_error(motif_proximity(anchors = anchors), regexp = "either")
})

test_that("data.frame input requires an explicit universe", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  df <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = FALSE)

  expect_error(suppressMessages(motif_proximity(df, anchors)),
               regexp = "universe")
  U <- GenomicRanges::GRanges(
    c("chr1", "chr2", "chr3"),
    IRanges::IRanges(1, c(6000, 4000, 3000)))
  r <- suppressMessages(motif_proximity(df, anchors, universe = U, qvalue = 1))
  expect_equal(nrow(r), 1L)
})

test_that("U* restriction drops hits on anchorless chromosomes", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  ## also plant on chr2 (which has no anchors) -- those hits must be dropped
  g <- plant_at(g, "chr2", seq(400, 3600, by = 200))
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  expect_message(motif_proximity(hits, anchors, qvalue = 1),
                 regexp = "dropped")
  r <- suppressMessages(motif_proximity(hits, anchors, qvalue = 1))
  ## nhits counts only chr1 (anchor-bearing) hits, fewer than all hits.
  expect_lt(r$nhits, length(hits))
})

test_that("mismatched seqnames raise an informative error", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("Chr1",   # capital C: does not match chr1
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  expect_error(suppressMessages(motif_proximity(hits, anchors)),
               regexp = "seqnames|seqlevels")
})

test_that("orientation: downstream-planted motif is found downstream only", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  ## Plant 5 bp downstream of '+' anchors.
  g <- plant_at(g, "chr1", pos, offset = 5L)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)

  down <- suppressMessages(motif_proximity(hits, anchors,
                            orientation = "downstream", qvalue = 1))
  up   <- suppressMessages(motif_proximity(hits, anchors,
                            orientation = "upstream", qvalue = 1))
  expect_lt(down$pvalue, 1e-4)
  expect_gt(down$enrichment, up$enrichment)
  expect_gt(up$pvalue, down$pvalue)

  ## directional requires stranded anchors
  unstr <- anchors
  GenomicRanges::strand(unstr) <- "*"
  expect_error(motif_proximity(hits, unstr, orientation = "downstream"),
               regexp = "stranded")
})

test_that("dist.to.nearest is signed when anchors are stranded", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos, offset = 5L)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)
  r <- suppressMessages(motif_proximity(hits, anchors, qvalue = 1))
  d <- r$dist.to.nearest[[1]]
  ## Motif sits downstream of '+' anchors, so signed distances skew positive.
  expect_true(any(d > 0))
  expect_gt(median(d, na.rm = TRUE), 0)
})

test_that("collapse_hits_per_cluster keeps one best hit per overlap cluster", {
  skip_no_granges()
  ## 3 overlapping hits (one cluster) + 1 separate hit, same motif.
  gr <- GenomicRanges::GRanges(
    "chr1",
    IRanges::IRanges(start = c(100, 103, 106, 500), width = 8),
    strand = "+")
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(
    motif = rep("m", 4), score = c(1, 5, 2, 9))
  out <- universalmotif:::collapse_hits_per_cluster(gr, ignore.strand = TRUE)
  expect_equal(length(out), 2L)
  ## best score of the overlapping cluster (5) is kept, not 1 or 2.
  expect_true(5 %in% S4Vectors::mcols(out)$score)
  expect_true(9 %in% S4Vectors::mcols(out)$score)
  expect_false(1 %in% S4Vectors::mcols(out)$score)
})

test_that("plot_motif_proximity returns a ggplot and guards bad input", {
  skip_no_granges()
  g <- make_genome()
  pos <- seq(500, 5500, by = 250)
  g <- plant_at(g, "chr1", pos)
  anchors <- GenomicRanges::GRanges("chr1",
                IRanges::IRanges(pos, width = 1), strand = "+")
  m <- create_motif("TTGACATA", name = "x")
  hits <- scan_sequences2(m, g, pvalue = 1e-3, return.granges = TRUE)
  r <- suppressMessages(motif_proximity(hits, anchors, qvalue = 1))

  gg <- plot_motif_proximity(r)
  expect_s3_class(gg, "ggplot")
  geoms <- vapply(gg$layers, function(l) class(l$geom)[1], character(1))
  expect_true(any(grepl("GeomBar|GeomHistogram", geoms)))

  expect_error(plot_motif_proximity(universalmotif:::empty_proximity_result()),
               regexp = "0 rows")
  expect_error(plot_motif_proximity(data.frame(motif = "x")),
               regexp = "dist.to.nearest")
})
