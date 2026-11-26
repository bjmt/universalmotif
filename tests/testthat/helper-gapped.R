## Shared fixtures for the gapped scanning tests in test_scan_sequences.R and
## test_scan_sequences_lite.R.

gapped_fixture <- function() {
  seqs <- Biostrings::DNAStringSet(c(
    s1 = "GCTAAAGACATGACGAGTCAATTACATAACATCACGTGACACGTCAGCTGACGTCAACGAAACTTG",
    s2 = "TTGGCCCAAGGTTTCACGGAGTGTGAATCGCTTAATCCTGAACCTGGGTTANNAGTAAGTGTG",
    s3 = "ATGCATACGCCTTGACTGTCATTACTTGCTGTGACGGGTCATGTCCACCCCATCG"))
  pcm <- function(cons, name) {
    l <- strsplit(cons, "")[[1]]
    m <- vapply(l, function(x) 1 + 3 * (c("A", "C", "G", "T") == x), numeric(4))
    rownames(m) <- c("A", "C", "G", "T")
    create_motif(unname(m), type = "PCM", alphabet = "DNA", pseudocount = 1,
                 name = name)
  }
  ## An ungapped motif, then one whose gap can have size zero, then one with
  ## two gaps and four combinations of gap sizes.
  motifs <- list(
    pcm("CACGTG", "u"),
    add_gap(pcm("TGACGTCA", "g"), gaploc = 4, mingap = 0, maxgap = 2),
    add_gap(pcm("AGGTCAGGA", "g2"), gaploc = c(3, 6), mingap = c(1, 0),
            maxgap = c(2, 1)))
  list(seqs = seqs, motifs = motifs)
}

read_hits <- function(text) {
  x <- read.table(text = text, stringsAsFactors = FALSE,
    col.names = c("motif", "sequence", "start", "end", "strand", "score",
                  "match", "pvalue"))
  x$start <- as.integer(x$start)
  x$end <- as.integer(x$end)
  x
}
