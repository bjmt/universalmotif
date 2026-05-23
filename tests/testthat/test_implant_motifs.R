context("implant_motifs()")

## Helpers ----------------------------------------------------------------

## Deterministic 6-bp motif: one-hot CACGTG ("E-box" style).
make_fixed_motif <- function(consensus = "CACGTG") {
  alph <- c("A", "C", "G", "T")
  mat  <- matrix(0, nrow = 4, ncol = nchar(consensus),
                 dimnames = list(alph, NULL))
  for (j in seq_len(nchar(consensus))) {
    let <- substr(consensus, j, j)
    mat[let, j] <- 1
  }
  create_motif(mat, alphabet = "DNA", name = "fixed", type = "PPM")
}

## A simple random DNA motif for shape tests.
make_random_motif <- function(width = 8, seed = 1) {
  set.seed(seed)
  create_motif(width, alphabet = "DNA")
}

make_seqs <- function(n, len, gc = 0.4, seed = 1, alph = "DNA") {
  set.seed(seed)
  letters_ <- if (alph == "DNA") c("A", "C", "G", "T") else c("A", "C", "G", "U")
  prob <- c((1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2)
  out <- vapply(seq_len(n), function(i)
    paste0(sample(letters_, len, replace = TRUE, prob = prob), collapse = ""),
    character(1))
  if (alph == "DNA") Biostrings::DNAStringSet(out)
  else               Biostrings::RNAStringSet(out)
}

## Tests ------------------------------------------------------------------

test_that("output shape and class preserved", {
  m <- make_random_motif(8)
  s <- make_seqs(10, 200, seed = 1)
  out <- implant_motifs(m, s, n.per.seq = 2L, rng.seed = 1)
  expect_s4_class(out, "DNAStringSet")
  expect_equal(length(out), length(s))
  expect_equal(as.integer(width(out)), as.integer(width(s)))
  expect_equal(names(out), names(s))
})

test_that("return.indices = TRUE produces the documented columns", {
  m <- make_random_motif(8)
  s <- make_seqs(5, 200, seed = 2)
  df <- implant_motifs(m, s, n.per.seq = 3L, return.indices = TRUE,
                       rng.seed = 1)
  expect_s3_class(df, "data.frame")
  expect_equal(names(df),
               c("sequence.i", "motif.i", "start", "width", "strand", "planted"))
  expect_equal(nrow(df), 5L * 3L)
  expect_true(all(df$width == 8L))
  expect_true(all(df$motif.i == 1L))
  expect_true(all(df$strand %in% c("+", "-")))
})

test_that("fixed-motif round trip: scan recovers every implant (+ strand)", {
  m <- make_fixed_motif("CACGTG")
  s <- make_seqs(20, 200, gc = 0.4, seed = 3)
  out <- implant_motifs(m, s, n.per.seq = 1L, strand = "+", rng.seed = 1)
  df  <- implant_motifs(m, s, n.per.seq = 1L, strand = "+",
                        return.indices = TRUE, rng.seed = 1)
  ## Every planted instance should literally appear as "CACGTG" at df$start.
  for (k in seq_len(nrow(df))) {
    seq_i <- df$sequence.i[k]
    st    <- df$start[k]
    sub   <- as.character(Biostrings::subseq(out[seq_i], start = st, width = 6L))
    expect_equal(sub, "CACGTG")
  }
})

test_that("n.per.seq mode produces N * n.per.seq rows", {
  m <- make_random_motif(8)
  s <- make_seqs(10, 200, seed = 4)
  df <- implant_motifs(m, s, n.per.seq = 3L, return.indices = TRUE,
                       rng.seed = 1)
  expect_equal(nrow(df), 30L)
  expect_setequal(sort(unique(df$sequence.i)), seq_len(10))
  expect_equal(sort(table(df$sequence.i)), as.table(rep(3L, 10),
               which = NULL), check.attributes = FALSE)
})

test_that("rate mode yields roughly Poisson(rate * width) implants total", {
  m <- make_random_motif(8)
  s <- make_seqs(500, 200, seed = 5)
  df <- implant_motifs(m, s, rate = 5e-3,
                       return.indices = TRUE, rng.seed = 1)
  expected <- 5e-3 * 200 * 500
  expect_gt(nrow(df), 0.5 * expected)
  expect_lt(nrow(df), 1.5 * expected)
})

test_that("positions mode honours user starts exactly", {
  m <- make_random_motif(8)
  s <- make_seqs(4, 200, seed = 6)
  pos <- list(c(10L, 50L), c(20L), c(100L, 150L, 180L), c(1L))
  df  <- implant_motifs(m, s, positions = pos, return.indices = TRUE,
                        rng.seed = 1)
  expect_equal(df$start[df$sequence.i == 1L], pos[[1]])
  expect_equal(df$start[df$sequence.i == 2L], pos[[2]])
  expect_equal(df$start[df$sequence.i == 3L], pos[[3]])
  expect_equal(df$start[df$sequence.i == 4L], pos[[4]])
})

test_that("strand = 'both' produces both + and -", {
  m <- make_random_motif(8)
  s <- make_seqs(50, 200, seed = 7)
  df <- implant_motifs(m, s, n.per.seq = 2L, strand = "both",
                       return.indices = TRUE, rng.seed = 1)
  expect_true(any(df$strand == "+"))
  expect_true(any(df$strand == "-"))
})

test_that("strand = '+' produces only + strand", {
  m <- make_random_motif(8)
  s <- make_seqs(50, 200, seed = 8)
  df <- implant_motifs(m, s, n.per.seq = 2L, strand = "+",
                       return.indices = TRUE, rng.seed = 1)
  expect_true(all(df$strand == "+"))
})

test_that("centre.bias clusters implants near sequence centres", {
  m <- make_random_motif(8)
  s <- make_seqs(500, 1000, seed = 9)
  df.unif <- implant_motifs(m, s, n.per.seq = 1L, centre.bias = 1L,
                            return.indices = TRUE, rng.seed = 1)
  df.cb   <- implant_motifs(m, s, n.per.seq = 1L, centre.bias = 8L,
                            return.indices = TRUE, rng.seed = 1)
  ## Distance from sequence midpoint should be smaller for centre-biased.
  midpt <- 500
  d.unif <- abs(df.unif$start - midpt)
  d.cb   <- abs(df.cb$start   - midpt)
  expect_lt(mean(d.cb), mean(d.unif))
})

test_that("min.spacing prevents overlaps and respects the gap", {
  m <- make_random_motif(8)
  ## One long sequence + many implants
  s <- make_seqs(1, 5000, seed = 10)
  df <- implant_motifs(m, s, n.per.seq = 20L, min.spacing = 50L,
                       return.indices = TRUE, rng.seed = 1)
  starts <- sort(df$start)
  ends   <- starts + df$width[order(df$start)] - 1L
  if (length(starts) > 1L) {
    gaps <- starts[-1] - ends[-length(ends)] - 1L
    expect_true(all(gaps >= 50L))
  }
})

test_that("determinism: same rng.seed -> identical output and indices", {
  m <- make_random_motif(8)
  s <- make_seqs(10, 200, seed = 11)
  o1 <- implant_motifs(m, s, n.per.seq = 2L, rng.seed = 42L)
  o2 <- implant_motifs(m, s, n.per.seq = 2L, rng.seed = 42L)
  ## sanity: the implant should actually change the sequences
  expect_false(identical(as.character(o1), as.character(s)))
  expect_equal(as.character(o1), as.character(o2))
  d1 <- implant_motifs(m, s, n.per.seq = 2L, return.indices = TRUE,
                       rng.seed = 42L)
  d2 <- implant_motifs(m, s, n.per.seq = 2L, return.indices = TRUE,
                       rng.seed = 42L)
  expect_equal(d1, d2)
})

test_that("AA / non-DNA-RNA motif rejection", {
  m <- create_motif(5, alphabet = "AA")
  s <- make_seqs(2, 100, seed = 12)
  expect_error(implant_motifs(m, s), regexp = "DNA/RNA")
})

test_that("DNA motif + RNA sequences mismatch rejected", {
  m <- make_random_motif(8)
  s <- make_seqs(2, 100, seed = 13, alph = "RNA")
  expect_error(implant_motifs(m, s), regexp = "alphabet")
})

test_that("capacity overflow emits a warning", {
  m <- make_random_motif(8)
  ## Very short sequence, demand many non-overlapping wide-spaced implants
  s <- make_seqs(1, 50, seed = 14)
  expect_warning(
    implant_motifs(m, s, n.per.seq = 20L, min.spacing = 30L,
                   max.retries = 5L, rng.seed = 1),
    regexp = "could not be placed"
  )
})

test_that("multiple motifs: motif.i covers all input motifs", {
  motifs <- list(make_random_motif(8, seed = 1),
                 make_random_motif(10, seed = 2),
                 make_random_motif(6, seed = 3))
  s <- make_seqs(100, 200, seed = 15)
  df <- implant_motifs(motifs, s, n.per.seq = 1L, return.indices = TRUE,
                       rng.seed = 1)
  expect_setequal(sort(unique(df$motif.i)), 1:3)
  ## width column should reflect the actual motif chosen
  widths_by_id <- c(8L, 10L, 6L)
  expect_true(all(df$width == widths_by_id[df$motif.i]))
})
