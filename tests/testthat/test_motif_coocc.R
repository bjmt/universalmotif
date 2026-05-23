context("motif_coocc()")

## Helpers ----------------------------------------------------------------

## Fixed one-hot motif. Default width 12 keeps random matches at
## scan p ~ 1e-3 negligible (random rate ~ 4^-12 = 6e-8 per bp).
make_fixed <- function(consensus, name = consensus) {
  alph <- c("A", "C", "G", "T")
  mat  <- matrix(0, nrow = 4, ncol = nchar(consensus),
                 dimnames = list(alph, NULL))
  for (j in seq_len(nchar(consensus)))
    mat[substr(consensus, j, j), j] <- 1
  create_motif(mat, alphabet = "DNA", name = name, type = "PPM")
}

make_seqs <- function(n, len, gc = 0.4, seed = 1, alph = "DNA") {
  set.seed(seed)
  letters_ <- if (alph == "DNA") c("A","C","G","T") else c("A","C","G","U")
  prob <- c((1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2)
  out <- vapply(seq_len(n), function(i)
    paste0(sample(letters_, len, replace = TRUE, prob = prob),
           collapse = ""), character(1))
  if (alph == "DNA") Biostrings::DNAStringSet(out)
  else               Biostrings::RNAStringSet(out)
}

## Build a synthetic hit table from per-motif sequence-presence lists.
## Each motif gets a list of (sequence.i, start) tuples.
fake_hits <- function(presence_list, starts_list = NULL) {
  rows <- list()
  for (i in seq_along(presence_list)) {
    si <- presence_list[[i]]
    if (length(si) == 0L) next
    st <- if (!is.null(starts_list)) starts_list[[i]] else rep(50L, length(si))
    rows[[length(rows) + 1L]] <- data.frame(
      motif.i    = rep(as.integer(i), length(si)),
      sequence.i = as.integer(si),
      start      = as.integer(st),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

## Tests ------------------------------------------------------------------

test_that("output shape and column types", {
  motifs <- list(make_fixed("CACGTGCACGTG", "M1"),
                 make_fixed("AAATTTAAATTT", "M2"),
                 make_fixed("GGGCCCGGGCCC", "M3"))
  s <- make_seqs(50, 300, seed = 1)
  co <- motif_coocc(motifs, s, pvalue = 1e-3)
  expect_s3_class(co, "data.frame")
  expect_equal(names(co),
               c("motif_a", "motif_b", "a_only", "b_only", "both",
                 "neither", "odds_ratio", "pvalue", "qvalue"))
  expect_equal(nrow(co), choose(3L, 2L))
  expect_true(all(co$a_only + co$b_only + co$both + co$neither ==
                  length(s)))
})

test_that("strong-signal pair detected (synthetic hit table)", {
  ## 200 sequences. M1 in seqs 1..50 + 51..60 (60 total). M2 in seqs
  ## 1..50 + 61..70 (60 total). Both = 50, a_only = 10, b_only = 10,
  ## neither = 130. Under independence the (both) cell would average
  ## ~18; observed = 50, so strong over-co-occurrence.
  hits <- fake_hits(list(c(1:50, 51:60), c(1:50, 61:70)))
  motifs <- list(make_fixed("AAAAAAAAAAAA", "M1"),
                 make_fixed("CCCCCCCCCCCC", "M2"))
  co <- motif_coocc(motifs, hits = hits, n.sequences = 200L)
  expect_equal(nrow(co), 1L)
  expect_equal(co$both,     50L)
  expect_equal(co$a_only,   10L)
  expect_equal(co$b_only,   10L)
  expect_equal(co$neither, 130L)
  expect_gt(co$odds_ratio, 1)
  expect_lt(co$pvalue, 1e-3)
  expect_lt(co$qvalue, 1e-3)
})

test_that("no-signal pair has qvalue not significant (synthetic)", {
  ## Independent draws -- M1 in 40 random seqs, M2 in 40 random seqs.
  set.seed(31)
  S1 <- sample.int(200, 40)
  S2 <- sample.int(200, 40)
  hits <- fake_hits(list(S1, S2))
  motifs <- list(make_fixed("AAAAAAAAAAAA", "M1"),
                 make_fixed("CCCCCCCCCCCC", "M2"))
  co <- motif_coocc(motifs, hits = hits, n.sequences = 200L)
  expect_gt(co$qvalue, 0.05)
})

test_that("hit-table path agrees with internal-scan path", {
  motifs <- list(make_fixed("CACGTGCACGTG", "M1"),
                 make_fixed("AAATTTAAATTT", "M2"),
                 make_fixed("GGGCCCGGGCCC", "M3"))
  s <- make_seqs(60, 300, seed = 4)
  hits <- scan_sequences2(motifs, s, pvalue = 1e-3, return.granges = FALSE)
  co.internal <- motif_coocc(motifs, s, pvalue = 1e-3)
  co.hits     <- motif_coocc(motifs, hits = hits, n.sequences = length(s))
  ord <- function(x) x[order(x$motif_a, x$motif_b), ]
  expect_equal(ord(co.internal)[, c("motif_a","motif_b","a_only","b_only",
                                    "both","neither")],
               ord(co.hits)[, c("motif_a","motif_b","a_only","b_only",
                                "both","neither")],
               check.attributes = FALSE)
})

test_that("spatial mode adds descriptive columns without changing Fisher", {
  ## 40 sequences. In all 40, M1 has a hit at start=100 and M2 has a hit
  ## at start=120 (dist=20).
  hits <- rbind(
    data.frame(motif.i = 1L, sequence.i = 1:40,  start = 100L),
    data.frame(motif.i = 2L, sequence.i = 1:40,  start = 120L)
  )
  motifs <- list(make_fixed("AAAAAAAAAAAA", "M1"),
                 make_fixed("CCCCCCCCCCCC", "M2"))
  co_all   <- motif_coocc(motifs, hits = hits, n.sequences = 40L)
  co_close <- motif_coocc(motifs, hits = hits, n.sequences = 40L,
                          max.distance = 50L)
  co_tight <- motif_coocc(motifs, hits = hits, n.sequences = 40L,
                          max.distance = 10L)

  ## `both` is the same across all three -- the Fisher test doesn't
  ## see the spatial filter.
  expect_equal(co_all$both,   40L)
  expect_equal(co_close$both, 40L)
  expect_equal(co_tight$both, 40L)
  ## p-value is identical across all three.
  expect_equal(co_all$pvalue, co_close$pvalue)
  expect_equal(co_all$pvalue, co_tight$pvalue)

  ## Spatial mode adds the two descriptive columns.
  expect_true(all(c("both.clustered", "median.distance") %in% names(co_close)))
  expect_false("both.clustered" %in% names(co_all))

  ## With max.distance = 50 (>= 20), all 40 sequences qualify.
  expect_equal(co_close$both.clustered, 40L)
  expect_equal(co_close$median.distance, 20)

  ## With max.distance = 10, none qualify -- both.clustered = 0,
  ## median.distance = NA.
  expect_equal(co_tight$both.clustered, 0L)
  expect_true(is.na(co_tight$median.distance))
})

test_that("min.coocc filter yields NA pvalue/qvalue (synthetic)", {
  ## M1 in 20 seqs, M2 in DISJOINT 20 seqs -> both == 0.
  hits <- fake_hits(list(1:20, 21:40))
  motifs <- list(make_fixed("AAAAAAAAAAAA", "M1"),
                 make_fixed("CCCCCCCCCCCC", "M2"))
  co <- motif_coocc(motifs, hits = hits, n.sequences = 40L,
                    min.coocc = 1L)
  expect_equal(co$both, 0L)
  expect_true(is.na(co$pvalue))
  expect_true(is.na(co$qvalue))
})

test_that("internal-scan path rejects AA motifs", {
  m <- create_motif(5, alphabet = "AA")
  s <- make_seqs(2, 100, seed = 7)
  expect_error(motif_coocc(list(m, m), s),
               regexp = "DNA/RNA")
})

test_that("hit-table path accepts AA motifs (no alphabet rejection)", {
  m1 <- create_motif(5, alphabet = "AA", name = "AA1")
  m2 <- create_motif(5, alphabet = "AA", name = "AA2")
  ## 30 seqs, M1 in seqs {1..15, 20..25}, M2 in seqs {10..25}.
  hits <- data.frame(
    motif.i    = c(rep(1L, 15 + 6), rep(2L, 16)),
    sequence.i = c(1:15, 20:25, 10:25),
    start      = c(rep(50L, 21), rep(100L, 16)),
    stringsAsFactors = FALSE
  )
  co <- motif_coocc(list(m1, m2), hits = hits, n.sequences = 30L)
  expect_s3_class(co, "data.frame")
  expect_equal(nrow(co), 1L)
  expect_true(is.finite(co$pvalue))
  ## Sanity on counts: M1 in 21 seqs, M2 in 16 seqs, intersect = seqs 10:15
  ## union 20:25 = 12 seqs.
  expect_equal(co$both, 12L)
})

test_that("DNA motif + RNA sequences mismatch rejected (internal-scan)", {
  m <- make_fixed("CACGTGCACGTG", "M1")
  s <- make_seqs(5, 100, seed = 8, alph = "RNA")
  expect_error(motif_coocc(list(m, m), s), regexp = "alphabet")
})

test_that("self.pairs flag flips (i,i) inclusion", {
  motifs <- list(make_fixed("AAAAAAAAAAAA", "M1"),
                 make_fixed("CCCCCCCCCCCC", "M2"))
  hits <- fake_hits(list(1:10, 11:20))
  co.no  <- motif_coocc(motifs, hits = hits, n.sequences = 20L,
                        self.pairs = FALSE)
  co.yes <- motif_coocc(motifs, hits = hits, n.sequences = 20L,
                        self.pairs = TRUE)
  expect_equal(nrow(co.no),  1L)   # choose(2,2) = 1
  expect_equal(nrow(co.yes), 3L)   # 1 + 2 self-pairs
  expect_true(any(co.yes$motif_a == co.yes$motif_b))
})

test_that("BH ordering: q-values >= p-values; sorted ascending", {
  ## Build a hit table with random per-motif presence sets.
  motifs <- list(make_fixed("AAAAAAAAAAAA", "M1"),
                 make_fixed("CCCCCCCCCCCC", "M2"),
                 make_fixed("GGGGGGGGGGGG", "M3"),
                 make_fixed("TTTTTTTTTTTT", "M4"))
  n.seq <- 80L
  set.seed(31)
  presence <- lapply(seq_along(motifs),
                     function(i) sample.int(n.seq, sample(20:50, 1)))
  hits <- fake_hits(presence)
  co <- motif_coocc(motifs, hits = hits, n.sequences = n.seq)
  tested <- !is.na(co$pvalue)
  expect_true(all(co$qvalue[tested] >= co$pvalue[tested] -
                  sqrt(.Machine$double.eps)))
  ## Within tested rows, q must be sorted ascending in the returned order.
  expect_equal(co$qvalue[tested], sort(co$qvalue[tested]))
})

test_that("hits as GRanges round-trips correctly", {
  motifs <- list(make_fixed("CACGTGCACGTG", "M1"),
                 make_fixed("AAATTTAAATTT", "M2"))
  s <- make_seqs(40, 300, seed = 41)
  hits_df <- scan_sequences2(motifs, s, pvalue = 1e-3, return.granges = FALSE)
  hits_gr <- scan_sequences2(motifs, s, pvalue = 1e-3, return.granges = TRUE)
  co_df <- motif_coocc(motifs, hits = hits_df, n.sequences = length(s))
  co_gr <- motif_coocc(motifs, hits = hits_gr, n.sequences = length(s))
  expect_equal(co_df[, c("both","a_only","b_only","neither")],
               co_gr[, c("both","a_only","b_only","neither")])
})

