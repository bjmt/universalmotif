context("plot_hits()")

## A small hit table in the shape of scan_sequences_lite(return.granges =
## FALSE): s1 has two overlapping hits and an unstranded one, s2 one hit.
hits_df <- function() {
  data.frame(
    motif = c("A", "B", "A", "B"),
    motif.i = c(1L, 2L, 1L, 2L),
    sequence = c("s1", "s1", "s1", "s2"),
    sequence.i = c(1L, 1L, 1L, 2L),
    start = c(11, 16, 41, 6),
    end = c(20, 25, 50, 15),
    strand = c("+", "-", "*", "-"),
    score = c(9, 8, 7, 8),
    pvalue = c(1e-5, 1e-4, 1e-3, 1e-4),
    stringsAsFactors = FALSE
  )
}
LENGTHS <- c(s1 = 60, s2 = 40, s3 = 50)

layer_classes <- function(p) {
  unname(vapply(p$layers, function(l) class(l$geom)[1], character(1)))
}

## Built data of the first layer with the given geom.
built_layer <- function(p, geom) {
  b <- ggplot2::ggplot_build(p)
  b$data[[which(layer_classes(p) == geom)[1]]]
}

## Each hit's left, right, bottom and top, centre line and the x of its tip
## (NA for blocks), from the built polygon data.
hit_boxes <- function(p) {
  d <- built_layer(p, "GeomPolygon")
  do.call(rbind, lapply(split(d, d$group), function(v) {
    tip <- if (v$x[3] > v$x[2]) v$x[3] else if (v$x[6] < v$x[1]) v$x[6]
      else NA
    data.frame(left = min(v$x), right = max(v$x), bottom = min(v$y),
      top = max(v$y), y = v$y[3], tip = tip, fill = v$fill[1],
      alpha = v$alpha[1])
  }))
}

test_that("hits are drawn where they are, as arrows along their strands", {
  p <- plot_hits(hits_df(), LENGTHS, labels = FALSE)
  expect_s3_class(p, "ggplot")
  expect_equal(layer_classes(p)[1:2], c("GeomSegment", "GeomPolygon"))
  expect_equal(p$scales$get_scales("y")$labels, c("s1", "s2", "s3"))

  ## Each track's line runs the length of its sequence.
  seg <- built_layer(p, "GeomSegment")
  expect_equal(seg$x, rep(0.5, 3))
  expect_equal(seg$xend, unname(LENGTHS) + 0.5)

  box <- hit_boxes(p)
  expect_equal(nrow(box), 4)
  expect_equal(box$left, hits_df()$start - 0.5)
  expect_equal(box$right, hits_df()$end + 0.5)
  expect_equal(box$tip[1], 20.5)
  expect_equal(box$tip[2], 15.5)
  expect_true(is.na(box$tip[3]))
  ## The hit on s2 is centred on the s2 line.
  expect_equal(box$y[4], seg$y[2])
  ## Overlapping hits on s1 go in separate lanes.
  expect_false(box$y[1] == box$y[2])
  expect_equal(box$fill, c("#F8766D", "#00BFC4", "#F8766D", "#00BFC4"))
})

test_that("scan_sequences() and scan_sequences_lite() hits draw the same", {
  m <- create_motif("TTGACC", name = "W", nsites = 20)
  s <- Biostrings::DNAStringSet(c(a = "AAAAGGTCAAAATTGACCAAA",
    b = "AGGTCAACCCCCCCCC"))
  h1 <- scan_sequences(m, s, RC = TRUE, nthreads = 1, threshold = 0.9,
    threshold.type = "logodds", calc.pvals = FALSE)
  h2 <- scan_sequences_lite(m, s, pvalue = 1e-2, return.granges = FALSE)
  h2 <- h2[h2$score >= min(h1$score), ]
  ## scan_sequences() reports reverse strand hits with start after stop.
  expect_true(any(h1$start > h1$stop))
  b1 <- hit_boxes(plot_hits(h1, s))
  b2 <- hit_boxes(plot_hits(h2, s))
  cols <- c("left", "right", "tip")
  expect_equal(unname(as.matrix(b1[order(b1$left, b1$y), cols])),
    unname(as.matrix(b2[order(b2$left, b2$y), cols])))
  expect_equal(sort(b1$left), sort(c(5, 13, 2) - 0.5))
})

test_that("GRanges hits are accepted", {
  skip_if_not_installed("GenomicRanges")
  m <- create_motif("TTGACC", name = "W", nsites = 20)
  s <- Biostrings::DNAStringSet(c(a = "AAAAGGTCAAAATTGACCAAA",
    b = "AGGTCAACCCCCCCCC"))
  gr <- scan_sequences_lite(m, s, pvalue = 1e-2, return.granges = TRUE)
  df <- scan_sequences_lite(m, s, pvalue = 1e-2, return.granges = FALSE)
  expect_equal(hit_boxes(plot_hits(gr, s)), hit_boxes(plot_hits(df, s)))
})

test_that("lanes never put hits on top of one another", {
  set.seed(1)
  starts <- sample(1:400, 300, replace = TRUE)
  h <- data.frame(motif = ifelse(seq_along(starts) %% 2 == 0, "A", "B"),
    sequence = paste0("s", seq_along(starts) %% 3), start = starts,
    end = starts + 11, strand = ifelse(seq_along(starts) %% 4 == 0, "-", "+"))
  box <- hit_boxes(plot_hits(h))
  expect_equal(nrow(box), 300)
  expect_gt(length(unique(box$y)), 3)
  apart <- outer(box$right, box$left, "<=") | outer(box$left, box$right, ">=") |
    outer(box$top, box$bottom, "<=") | outer(box$bottom, box$top, ">=")
  diag(apart) <- TRUE
  expect_true(all(apart))
  ln <- plot_hits_lanes(c(0, 5, 12, 30, 31), c(10, 15, 20, 40, 33), 0)
  expect_equal(ln$lane, c(0L, 1L, 0L, 0L, 1L))
  expect_equal(ln$n, 2L)
  expect_equal(plot_hits_lanes(c(0, 10.5), c(10, 20), 1)$lane, c(0L, 1L))
})

test_that("tracks follow the sequences given, or the hits", {
  s <- Biostrings::DNAStringSet(c(s1 = strrep("A", 60), s2 = strrep("C", 40)))
  p <- plot_hits(hits_df(), s)
  expect_equal(built_layer(p, "GeomSegment")$xend, c(60.5, 40.5))
  ## Unnamed sequences are named by their position, as scan_sequences()
  ## names them.
  h <- hits_df()
  h$sequence <- c("1", "1", "1", "2")
  p <- plot_hits(h, unname(s))
  expect_equal(p$scales$get_scales("y")$labels, c("1", "2"))
  ## Without sequences, tracks run to their last hit in sequence.i order.
  h <- hits_df()
  h$sequence.i <- c(2L, 2L, 2L, 1L)
  p <- plot_hits(h)
  expect_equal(p$scales$get_scales("y")$labels, c("s2", "s1"))
  expect_equal(built_layer(p, "GeomSegment")$xend, c(15.5, 50.5))
  ## Tracks without hits are drawn when sequences are given.
  p <- plot_hits(hits_df()[0, ], LENGTHS)
  expect_equal(p$scales$get_scales("y")$labels, names(LENGTHS))
  expect_equal(nrow(built_layer(p, "GeomPolygon")), 0)
})

test_that("names and max.tracks pick the tracks", {
  p <- plot_hits(hits_df(), LENGTHS, names = c("s3", "s1"))
  expect_equal(p$scales$get_scales("y")$labels, c("s3", "s1"))
  expect_equal(nrow(hit_boxes(p)), 3)
  p <- plot_hits(hits_df(), LENGTHS, max.tracks = 2)
  expect_equal(p$scales$get_scales("y")$labels, c("s1", "s2"))
  expect_match(p$labels$subtitle, "first 2 of 3 sequences")
  expect_null(plot_hits(hits_df(), LENGTHS)$labels$subtitle)
  expect_error(plot_hits(hits_df(), LENGTHS, names = "nope"),
    "no sequences or ranges named: nope")
  expect_error(plot_hits(hits_df(), max.tracks = 0), "max.tracks")
})

test_that("motifs set the order, the colours and the logos", {
  m <- list(create_motif("CACGTG", name = "B"), create_motif("TTGAC",
    name = "A"))
  p <- plot_hits(hits_df(), LENGTHS, motifs = m)
  expect_equal(p$scales$get_scales("fill")$get_limits(), c("B", "A"))
  box <- hit_boxes(p)
  expect_equal(box$fill, c("#00BFC4", "#F8766D", "#00BFC4", "#F8766D"))
  ## Logos come from a layer of invisible points keyed by motif name.
  expect_true(!is.null(p$scales$get_scales("colour")))
  ## Motifs given by name filter the hits and have no logos.
  p <- plot_hits(hits_df(), LENGTHS, motifs = "A")
  expect_equal(nrow(hit_boxes(p)), 2)
  expect_null(p$scales$get_scales("colour"))
  ## Colours by motif follow motif.i, whichever tracks are drawn.
  h <- hits_df()
  h$motif.i <- c(3L, 12L, 3L, 12L)
  box <- hit_boxes(plot_hits(h, LENGTHS))
  expect_equal(unique(box$fill), c("#B79F00", "#FF64B0"))
  box <- hit_boxes(plot_hits(hits_df(), colours = c(A = "black",
    B = "red")))
  expect_equal(unique(box$fill), c("black", "red"))
  expect_error(plot_hits(hits_df(), colours = c(A = "black")),
    "no colours for: B")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_s3_class(ggplot2::ggplotGrob(plot_hits(hits_df(), LENGTHS,
    motifs = m, shade.by = "pvalue", size.by = "score")), "gtable")
})

test_that("shade.by and size.by scale the hits", {
  h <- hits_df()
  p <- plot_hits(h, LENGTHS, shade.by = "pvalue")
  box <- hit_boxes(p)
  v <- -log10(h$pvalue)
  frac <- (v - min(v)) / (max(v) - min(v))
  expect_equal(box$alpha, 0.25 + 0.75 * frac)
  expect_equal(p$scales$get_scales("alpha")$name, "-log10(pvalue)")
  p <- plot_hits(h, LENGTHS, size.by = "score")
  box <- hit_boxes(p)
  frac <- (h$score - 7) / 2
  expect_equal(box$top - box$bottom, 0.8 * (0.35 + 0.65 * frac))
  expect_equal(p$scales$get_scales("size")$name, "score")
  ## The same column for both gives a single legend.
  p <- plot_hits(h, LENGTHS, shade.by = "score", size.by = "score")
  g <- ggplot2::get_guide_data(p, "size")
  expect_true(all(c("alpha", "size") %in% colnames(g)))
  expect_error(plot_hits(h, shade.by = "nope"), "no `nope` column")
  expect_error(plot_hits(h, size.by = "motif"), "numeric column")
})

test_that("labels go on hits wide enough to hold them", {
  h <- data.frame(motif = c("GATA", "W-box"), sequence = "a",
    start = c(3, 40), end = c(30, 42), strand = c("-", "+"))
  p <- plot_hits(h, c(a = 60))
  expect_true("GeomText" %in% layer_classes(p))
  txt <- built_layer(p, "GeomText")
  expect_equal(txt$label, "GATA")
  ## Black text on the light red of the first motif.
  expect_equal(txt$colour, "black")
  expect_false("GeomText" %in% layer_classes(plot_hits(h, c(a = 60),
    labels = FALSE)))
})

test_that("legend can be turned off", {
  expect_equal(plot_hits(hits_df(), legend = FALSE)$theme$legend.position,
    "none")
})

test_that("ranges draw the hits that overlap them", {
  h <- data.frame(motif = c("A", "B", "A", "B"),
    sequence = c("chr1", "chr1", "chr1", "chr2"),
    start = c(1001, 1095, 5001, 101), end = c(1010, 1104, 5010, 110),
    strand = "+", pvalue = 1e-4)
  r <- data.frame(seqnames = c("chr1", "chr2"), start = c(1001, 51),
    end = c(1100, 200), name = c("peak1", "peak2"))
  p <- plot_hits(h, ranges = r)
  expect_equal(p$scales$get_scales("y")$labels, c("peak1", "peak2"))
  expect_equal(p$labels$x, "Position from the start of each range")
  box <- hit_boxes(p)
  expect_equal(nrow(box), 3)
  ## Positions count from the start of each range, and the hit running off
  ## the end of peak1 is cut there.
  expect_equal(box$left, c(0.5, 94.5, 50.5))
  expect_equal(box$right, c(10.5, 100.5, 60.5))
  ## A single range keeps genomic positions.
  p <- plot_hits(h, ranges = r[1, c("seqnames", "start", "end")])
  expect_equal(p$scales$get_scales("y")$labels, "chr1:1001-1100")
  expect_equal(p$labels$x, "Position on chr1")
  expect_equal(hit_boxes(p)$left, c(1000.5, 1094.5))
  expect_error(plot_hits(h, ranges = r[0, ]), "no rows")
  expect_error(plot_hits(h, ranges = data.frame(start = 1, end = 2)),
    "seqnames")
})

test_that("GRanges ranges are accepted", {
  skip_if_not_installed("GenomicRanges")
  h <- data.frame(motif = "A", sequence = "chr1", start = 1001, end = 1010,
    strand = "-")
  r <- GenomicRanges::GRanges("chr1", IRanges::IRanges(990, 1100))
  names(r) <- "peak1"
  p <- plot_hits(h, ranges = r)
  expect_equal(p$scales$get_scales("y")$labels, "peak1")
  expect_equal(hit_boxes(p)$left, 1000.5)
})

test_that("bad input is refused", {
  expect_error(plot_hits(list()), "must be a DataFrame")
  expect_error(plot_hits(hits_df()[, -1]), "`motif` column")
  expect_error(plot_hits(hits_df()[, setdiff(names(hits_df()), "end")]),
    "`end` column")
  expect_error(plot_hits(hits_df(), c(s1 = 15)), "run past the end of s1")
  expect_error(plot_hits(hits_df(), c(x = 15)), "none of the hits")
  expect_error(plot_hits(hits_df(), 1:3), "named numeric")
  expect_error(plot_hits(hits_df()[0, ]), "no rows to plot")
  expect_error(plot_hits(hits_df(), text.size = -1), "text.size")
  h <- hits_df()
  h$start[1] <- NA
  expect_error(plot_hits(h), "NA")
})
