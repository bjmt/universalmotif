#' Find motif binding sites in a set of sequences.
#'
#' DNA/RNA only.
#'
#' @param motifs List of motifs or a single motif.
#' @param sequences XStringSet object. List of sequences to scan
#' @param threshold Numeric. Logodds threshold.
#' @param threshold.type Character. One of 'logodds' and 'pvalue'.
#' @param RC Logical. Check reverse strand.
#' @param use.freq Numeric.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Site search results as a data.frame object.
#'
#' @details
#'    Benchmarking: scanning a 1 billion bp sequence with HMMorder = 1 and 
#'    RC = FALSE took about 4.5 minutes. Scanning 1000 sequences 1 million bp
#'    each took about 6 minutes. As for HMMorder = 0; for shorter sequences
#'    it is amazingly fast.. but when I tried with a 1 billion bp long
#'    sequence it took about 5 min. Not sure why it slows down so much.
#'
#'    Careful with large k. A motif with freqs 10 is 80 MB.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.6,
                           threshold.type = "logodds", RC = TRUE,
                           use.freq = 1, BPPARAM = bpparam()) {

  if (is.list(motifs)) {
    results <- lapply(motifs, scan_sequences, threshold = threshold,
                      threshold.type = threshold.type,
                      sequences = sequences, RC = RC, 
                      use.freq = use.freq, BPPARAM = BPPARAM)
    results <- do.call(rbind, results)
    rownames(results) <- seq_len(nrow(results))
    return(results)
  }

  if (missing(motifs) || missing(sequences)) {
    stop("missing 'motifs' and/or 'sequences'")
  }

  if (motifs["alphabet"] == "RNA") motifs <- switch_alph(motifs)
  if (sequences@elementType == "RNAString") {
    sequences <- DNAStringSet(sequences)
    RNA <- TRUE
  } else RNA <- FALSE


  mot.name <- motifs["name"]
  mot.mat <- convert_type(motifs, "PWM")["motif"]
  mot.pfm <- convert_type(motifs, "PPM")["motif"]
  mot.len <- ncol(mot.mat)
  max.score <- sum(apply(mot.mat, 2, max))

  bkg <- motifs["bkg"]
  names(bkg) <- DNA_BASES

  if (threshold.type == "pvalue") {
    threshold <- TFMpv2sc(mot.pfm, threshold, bkg, type = "PFM")
    threshold <- threshold / sum(apply(mot.mat, 2, max))
  }
  if (threshold < 0) stop("cannot have negative threshold")

  min.score <- threshold * 100
  min.score <- paste0(min.score, "%")

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))


  if (use.freq == 1) {

    motif <- motifs
    motif <- convert_type(motif, "PWM")
    motif.rc <- suppressWarnings(motif_rc(motif))
    motif <- motif["motif"]
    motif.rc <- motif.rc["motif"]

    parse_hits <- function(x, y, strand = "+") {
      if (length(x) == 0) {
        data.frame(motif = NULL, sequence = NULL, start = NULL, stop = NULL,
                   max.score = NULL, score.pct = NULL, match = NULL,
                   strand = NULL)
      } else {
        motif <- rep(mot.name, length(x))
        sequence = rep(y, length(x))
        start <- x@ranges@start
        stop <- start + (mot.len - 1)
        if (strand == "-") {
          tmp <- start
          start <- stop
          stop <- tmp
        }
        score <- x@elementMetadata$score 
        max.score <- rep(max.score, length(x))
        score.pct <- (score / max.score) * 100
        match <- as.character(x)

        data.frame(motif = motif, sequence = sequence, start = start, stop = stop,
                   score = score, max.score = max.score, score.pct = score.pct,
                   match = match, strand = rep(strand, length(x)))
      }
    }

    sequence.hits <- bplapply(seq_len(length(sequences)),
                              function(x) matchPWM(motif, sequences[[x]],
                                                   min.score = min.score,
                                                   with.score = TRUE),
                              BPPARAM = BPPARAM)

    sequence.hits <- bpmapply(parse_hits, sequence.hits, seq.names,
                              BPPARAM = BPPARAM, SIMPLIFY = FALSE)

    sequence.hits <- do.call(rbind, sequence.hits)
    if (nrow(sequence.hits) > 0) sequence.hits$strand <- "+"

    if (RC) {
      sequence.hits.rc <- bplapply(seq_len(length(sequences)),
                                   function(x) matchPWM(motif.rc, sequences[[x]],
                                                        min.score = min.score,
                                                        with.score = TRUE),
                                   BPPARAM = BPPARAM)

      sequence.hits.rc <- bpmapply(parse_hits, sequence.hits.rc, seq.names,
                                   "-", BPPARAM = BPPARAM, SIMPLIFY = FALSE)

      sequence.hits.rc <- do.call(rbind, sequence.hits.rc)
      
      if (nrow(sequence.hits.rc) > 0) sequence.hits.rc$strand <- "-"
      match <- sequence.hits.rc$match
      match <- as.character(reverseComplement(DNAStringSet(match)))
      sequence.hits.rc$match <- match

      sequence.hits <- rbind(sequence.hits, sequence.hits.rc)

    }

  } else {

    if (use.freq > 6) {
      stop("'use.freq' > 6 is not supported as this crashes R")
      # warning("with 'use.freq' > 6 R may crash; parallel cores disabled")
      # BPPARAM <- SerialParam()
    }
    results <- bpmapply(function(x, y) get_res_k(motifs, x, y, "+", use.freq,
                                                 threshold),
                        sequences, seq.names, SIMPLIFY = FALSE, BPPARAM = BPPARAM)

    if (RC) {
      results.rc <- bpmapply(function(x, y) get_res_k(motifs, x, y, "-", use.freq,
                                                      threshold),
                             reverseComplement(sequences), seq.names,
                             SIMPLIFY = FALSE, BPPARAM = BPPARAM)

      results <- c(results, results.rc)
    }

    sequence.hits <- do.call(rbind, results)

  } 

  if (is.null(sequence.hits)) {
    message("no matches found using current threshold for motif ", mot.name)
    return(NULL)
  } else if (nrow(sequence.hits) == 0) {
    message("no matches found using current threshold for motif ", mot.name)
    return(NULL)
  }

  sequence.hits <- sequence.hits[order(sequence.hits$score.pct,
                                       decreasing = TRUE), ]
  rownames(sequence.hits) <- seq_len(nrow(sequence.hits))

  if (RNA) {
    sequence.hits$match <- gsub("T", "U", sequence.hits$match)
  }

  sequence.hits$p.value <- vapply(sequence.hits$score,
                                  function(x) TFMsc2pv(mot.pfm, x, bkg, "PFM"),
                                  numeric(1))

  sequence.hits[, c(1:7, 10, 8:9)]

}

get_res_k <- function(motif, seq, seq.name, seqstrand, k, threshold) {

  sequence <- as.character(seq)
  sequence <- strsplit(sequence, "")[[1]]

  if (motif["pseudocount"] == 0) motif["pseudocount"] <- 0.0001
  if (length(motif["nsites"]) == 0) motif["nsites"] <- 100

  # warning: r can run crash here with k > 6
  # skip.bigmem <- TRUE
  # if (k > 6) {
    # if (!requireNamespace("bigmemory", quietly = TRUE)) {
      # skip.bigmem <- TRUE
    # } else {
      # skip.bigmem <- FALSE
      # scores <- bigmemory::as.big.matrix(motif@multifreq[[as.character(k)]])
      # max.score <- vector("numeric", ncol(scores))
      # for (i in colnames(scores)) {
        # scores[, i] <- ppm_to_pwm(scores[, i], motif["bkg"], motif["pseudocount"],
                                  # motif["nsites"])
        # max.score[i] <- sum(scores[, i])
      # # }
      # max.score <- sum(max.score)
    # }
  # } else if (k <= 6) {
    scores <- apply(motif@multifreq[[as.character(k)]], 2, ppm_to_pwmC,
                    nsites = motif["nsites"],
                    pseudocount = motif["pseudocount"])
    max.score <- sum(apply(scores, 2, max))
  # }
  min.score <- max.score * threshold

  max.len <- length(sequence)
  seq.mat <- matrix(ncol = max.len - k + 1, nrow = k)

  for (i in seq_len(k)) {
    to.remove <- i - 1
    if (to.remove == 0) {
      seq.mat[i, ] <- sequence[seq_len(ncol(seq.mat))]
      next
    }
    sequence.i <- sequence[-seq_len(to.remove)]
    seq.mat[i, ] <- sequence.i[seq_len(ncol(seq.mat))]
  }

  seqs <- DNA_to_int_k(as.character(seq.mat), k)
  # warning: r can crash here with k > 6
  # if (!skip.bigmem) {
    # to_keep <- scan_seq_internal_bigmem(seqs, scores@address, min.score)
  # } else {
    to_keep <- scan_seq_internal(seqs, scores, min.score)
  # }
  to_keep <- which(as.logical(to_keep))

  mot_len <- ncol(motif["motif"]) - k + 1
  name <- motif["name"]

  res <- lapply(to_keep, function(x) parse_k_res(x, sequence, seqs,
                mot_len, min.score, max.score, name, seq.name, scores,
                seqstrand, max.len, k))

  do.call(rbind, res)

}

parse_k_res <- function(to_keep, sequence, seqs, mot_len, min.score,
                        max.score, mot.name, seq.name, score.mat,
                        strand, seq.length, k) {

  if (strand == "+") {
    hit <- seqs[to_keep:(to_keep + mot_len - 1)]
    score <- score_seq(hit, score.mat)
    data.frame(motif = mot.name, sequence = seq.name, start = to_keep,
               stop = to_keep + mot_len + k - 2, score = score,
               max.score = max.score, score.pct = score / max.score * 100,
               match = paste(sequence[to_keep:(to_keep + mot_len + k - 2)],
                             collapse = ""), strand = strand)
  } else if (strand == "-") {
    hit <- seqs[to_keep:(to_keep + mot_len - 1)]
    score <- score_seq(hit, score.mat)
    start <- seq.length - to_keep + 1
    stop <- to_keep + mot_len
    stop <- seq.length - stop - k + 3
    match <- paste(sequence[to_keep:(to_keep + mot_len + k - 2)],
                   collapse = "")
    data.frame(motif = mot.name, sequence = seq.name, start = start,
               stop = stop, score = score, max.score = max.score,
               score.pct = score / max.score * 100,
               match = as.character(match), strand = strand)
  }
  
}
