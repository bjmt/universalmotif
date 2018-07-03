#' Find motif binding sites in a set of sequences.
#'
#' @param motifs List of motifs or a single motif.
#' @param sequences XStringSet object. List of sequences to scan
#' @param threshold Numeric. Logodds threshold.
#' @param RC Logical. Check reverse strand.
#' @param HMMorder Numeric.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Site search results as a data.frame object.
#'
#' @details
#'    Benchmarking: scanning a 1 billion bp sequence with HMMorder = 1 and 
#'    RC = FALSE took about 4.5 minutes. Scanning 1000 sequences 1 million bp
#'    each took about 6 minutes.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
scan_sequences <- function(motifs, sequences, threshold = 0.6,
                           RC = FALSE, HMMorder = 0, BPPARAM = bpparam()) {

  if (is.list(motifs)) {
    results <- lapply(motifs, scan_sequences, threshold = threshold,
                      sequences = sequences, RC = RC, HMMorder = 0,
                      BPPARAM = BPPARAM)
    results <- do.call(rbind, results)
    rownames(results) <- seq_len(nrow(results))
    return(results)
  }

  if (missing(motifs) || missing(sequences)) stop("missing 'motifs' and/or 'sequences'")

  min.score <- threshold * 100
  min.score <- paste0(min.score, "%")

  mot.name <- motifs["name"]
  mot.mat <- convert_type(motifs, "PWM")["motif"]
  mot.len <- ncol(mot.mat)
  max.score <- sum(apply(mot.mat, 2, max))

  seq.names <- names(sequences)
  if (is.null(seq.names)) seq.names <- seq_len(length(sequences))


  if (HMMorder == 0) {

    motif <- motifs
    motif <- convert_type(motif,"PWM")
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

      sequence.hits <- rbind(sequence.hits, sequence.hits.rc)

    }

  } else if (HMMorder == 1) {

    if (length(motifs["hmmfirst"]) <= 1) stop("no 1st-order HMM in motif ",
                                              mot.name)

    parse_1st_res <- function(to_keep, sequence, seqs, mot_len, min.score,
                              max.score, mot.name, seq.name, score.mat, strand,
                              seq.length) {
      if (strand == "+") {
        hit <- seqs[to_keep:(to_keep + mot_len - 1)]
        score <- score_seq(hit, score.mat)
        data.frame(motif = mot.name, sequence = seq.name, start = to_keep,
                   stop = to_keep + mot_len, score = score, max.score = max.score,
                   score.pct = score / max.score * 100,
                   match = paste(sequence[to_keep:(to_keep + mot_len)],
                                 collapse = ""),
                   strand = strand)
      } else if (strand == "-") {
        hit <- seqs[to_keep:(to_keep + mot_len - 1)]
        score <- score_seq(hit, score.mat)
        start <- seq.length - to_keep + 1
        stop <- to_keep + mot_len
        stop <- seq.length - stop + 1
        match <- paste(sequence[to_keep:(to_keep + mot_len)], collapse = "")
        match <- reverseComplement(DNAString(match))
        data.frame(motif = mot.name, sequence = seq.name, start = start,
                   stop = stop, score = score, max.score = max.score,
                   score.pct = score / max.score * 100,
                   match = as.character(match), strand = strand)
      }
    }

    get_res <- function(motif, seq, seq.name, seqstrand) {

      sequence <- as.character(seq)
      sequence <- strsplit(sequence, "")[[1]]

      if (motif["pseudocount"] == 0) motif["pseudocount"] <- 0.8
      scores <- apply(motif["hmmfirst"], 2, ppm_to_pwm,
                      nsites = motif["nsites"],
                      pseudocount = motif["pseudocount"])
      max.score <- sum(apply(scores, 2, max))
      min.score <- max.score * threshold

      seq.mat <- matrix(ncol = length(sequence) - 1, nrow = 2)
      seq.mat[1, ] <- sequence[-length(sequence)]
      seq.mat[2, ] <- sequence[-1]

      lookup <- as.integer(0:15)
      names(lookup) <- DNA_DI

      seqs <- DNA_to_int_di(as.character(seq.mat))

      to_keep <- scan_1st_order(seqs, scores, min.score)
      to_keep <- which(as.logical(to_keep))

      mot_len <- ncol(motif["motif"]) - 1
      name <- motif["name"]
      res <- lapply(to_keep,
                    function(x) parse_1st_res(x, sequence, seqs, mot_len,
                                              min.score, max.score, name,
                                              seq.name, scores, seqstrand,
                                              length(sequence)))
      do.call(rbind, res)

    }

    results <- bpmapply(function(x, y) get_res(motifs, x, y, "+"),
                        sequences, seq.names, SIMPLIFY = FALSE, BPPARAM = BPPARAM)

    if (RC) {
      results.rc <- bpmapply(function(x, y) get_res(motifs, x, y, "-"),
                             reverseComplement(sequences), seq.names,
                             SIMPLIFY = FALSE, BPPARAM = BPPARAM)
      results <- c(results, results.rc)
    }

    sequence.hits <- do.call(rbind, results)

  }  else if (HMMorder == 2) {
  
    if (length(motifs["hmmsecond"]) <= 1) stop("no 2nd-order HMM in motif ",
                                               mot.name)

    parse_2nd_res <- function(to_keep, sequence, seqs, mot_len, min.score,
                              max.score, mot.name, seq.name, score.mat,
                              strand, seq.length) {
      if (strand == "+") {
        hit <- seqs[to_keep:(to_keep + mot_len - 1)]
        score <- score_seq(hit, score.mat)
        data.frame(motif = mot.name, sequence = seq.name, start = to_keep,
                   stop = to_keep + mot_len + 1, score = score,
                   max.score = max.score,
                   score.pct = score / max.score * 100,
                   match = paste(sequence[to_keep:(to_keep + mot_len + 1)],
                                 collapse = ""), strand = strand)
      } else if (strand == "-") {
        hit <- seqs[to_keep:(to_keep + mot_len - 1)]
        score <- score_seq(hit, score.mat)
        start <- seq.length - to_keep + 1
        stop <- to_keep + mot_len
        stop <- seq.length - stop
        match <- paste(sequence[to_keep:(to_keep + mot_len + 1)],
                       collapse = "")
        match <- reverseComplement(DNAString(match))
        data.frame(motif = mot.name, sequence = seq.name, start = start,
                   stop = stop, score = score, max.score = max.score,
                   score.pct = score / max.score * 100,
                   match = as.character(match), strand = strand)
      }
    }

    get_res_2nd <- function(motif, seq, seq.name, seqstrand) {
    
      sequence <- as.character(seq)
      sequence <- strsplit(sequence, "")[[1]]

      if (motif["pseudocount"] == 0) motif["pseudocount"] <- 0.8
      scores <- apply(motif["hmmsecond"], 2, ppm_to_pwm,
                      nsites = motif["nsites"],
                      pseudocount = motif["pseudocount"])
      max.score <- sum(apply(scores, 2, max))
      min.score <- max.score * threshold

      seq.mat <- matrix(ncol = length(sequence) - 2, nrow = 3)
      seq.mat[1, ] <- sequence[-c(length(sequence) - 1, length(sequence))]
      seq.mat[2, ] <- sequence[-c(1, length(sequence))]
      seq.mat[3, ] <- sequence[-c(1, 2)]

      lookup <- as.integer(0:63)
      names(lookup) <- DNA_TRI

      seqs <- DNA_to_int_tri(as.character(seq.mat))

      to_keep <- scan_2nd_order(seqs, scores, min.score)
      to_keep <- which(as.logical(to_keep))

      mot_len <- ncol(motif["motif"]) - 2
      name <- motif["name"]
      res <- lapply(to_keep,
                    function(x) parse_2nd_res(x, sequence, seqs, mot_len,
                                              min.score, max.score, name,
                                              seq.name, scores, seqstrand,
                                              length(sequence)))
      do.call(rbind, res)
    
    }

    results <- bpmapply(function(x, y) get_res_2nd(motifs, x, y, "+"),
                        sequences, seq.names, SIMPLIFY = FALSE, BPPARAM = BPPARAM)

    if (RC) {
      results.rc <- bpmapply(function(x, y) get_res_2nd(motifs, x, y, "-"),
                             reverseComplement(sequences), seq.names,
                             SIMPLIFY = FALSE, BPPARAM = BPPARAM)
      results <- c(results, results.rc)
    }

    sequence.hits <- do.call(rbind, results)

  } else stop("incorrect 'HMMorder'")

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

  sequence.hits

}
