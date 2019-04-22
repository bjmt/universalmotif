#' Utility functions.
#'
#' Utility functions have been split into two categories: those related to
#' motifs ?`utils-motif`, and those related to sequences ?`utils-sequence`.
#'
#' @seealso [utils-motif], [utils-sequence]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name utilities
NULL

# PUBLIC UTILS LIST
#
#   motif-related:
#
#     get_matches()
#     score_match()
#     motif_score()
#     summarise_motifs()
#     make_DBscores()
#     ppm_to_icm()
#     icm_to_ppm()
#     pcm_to_ppm()
#     ppm_to_pcm()
#     ppm_to_pwm()
#     pwm_to_ppm()
#     position_icscore()
#     get_consensus()
#     consensus_to_ppm()
#     consensus_to_ppmAA()
#     get_consensusAA()
#
#   sequence-related:
#     
#     get_klets()
#     count_klets()
#     shuffle_string(x, k = 1, method = "euler")

# INTERNAL CONSTANTS -----------------------------------------------------------

DNA_DI <- c("AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT")

DNA_TRI <- c("AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG",
             "AGT","ATA","ATC","ATG","ATT","CAA","CAC","CAG","CAT","CCA","CCC",
             "CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT","GAA",
             "GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT",
             "GTA","GTC","GTG","GTT","TAA","TAC","TAG","TAT","TCA","TCC","TCG",
             "TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT")

# TYPES:
#
#        NULL   0
#      symbol   1
#     closure   3
# environment   4
#    promises   5
#    logical   10
#    integer   13
#       real   14
#    complex   15
#     string   16
#        ...   17
#        ANY   18
#    generic   19
# expression   20
#  byte code   21
#    ext ptr   22
#        raw   24
#         S4   25

TYPE_LOGI <- 10L
TYPE_NUM  <- 14L
TYPE_CHAR <- 16L
TYPE_S4   <- 25L

UNIVERSALMOTIF_SLOTS <- c(

  "name",
  "altname",
  "family",
  "organism",
  "motif",
  "alphabet",
  "type",
  "icscore",
  "nsites",
  "pseudocount",
  "bkg",
  "bkgsites",
  "consensus",
  "strand",
  "pval",
  "qval",
  "eval",
  "multifreq",
  "extrainfo"

)

COMPARE_METRICS <- c("PCC", "MPCC", "EUCL", "MEUCL", "SW", "MSW", "KL", "MKL")

# INTERNAL UTILITIES ----------------------------------------------------------- 

.internal_convert <- function(motifs, class = NULL) {

  if (is.null(class)) {

    CLASS <- class(motifs)
    CLASS_PKG <- attributes(CLASS)$package
    CLASS_IN <- collapse_cpp(c(CLASS_PKG, "-", CLASS))

    CLASS_IN

  } else {

    if (length(class) == 1 && class[1] != "universalmotif-universalmotif") {

      tryCatch(motifs <- convert_motifs(motifs, class),
               error = function(e) message("motifs converted to class 'universalmotif'"))

    } else if (length(class) > 1)
      message("motifs converted to class 'universalmotif'")

    motifs

  }

}

# for a motif of length 4, the transition matrix is something like this:
#       bkg pos1 pos2 pos3 pos4
#  bkg    0    1    0    0    0
# pos1    0    0    1    0    0
# pos2    0    0    0    1    0
# pos3    0    0    0    0    1
# pos4    1    0    0    0    0

# Inspired from S4Vectors::wmsg
wmsg2 <- function(..., exdent = 0, indent = 0)
  paste0(strwrap(paste0(..., collapse = ""), exdent = exdent, indent = indent),
         collapse = "\n")

lapply_ <- function(X, FUN, ..., BP = FALSE, PB = FALSE) {

  FUN <- match.fun(FUN)

  if (!BP) {

    if (!PB) {

      out <- lapply(X, FUN, ...)

    } else {

      out <- vector("list", length(X))
      max <- length(X)
      print_pb(0)
      if (is.list(X)) {
        for (i in seq_along(X)) {
          out[[i]] <- do.call(FUN, list(X[[i]], ...))
          update_pb(i, max)
        }
      } else {
        for (i in seq_along(X)) {
          out[[i]] <- do.call(FUN, list(X[i], ...))
          update_pb(i, max)
        }
      }

    }

  } else {

    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      out <- BiocParallel::bplapply(X, FUN, ...)
    } else {
      stop("'BiocParallel' is not installed")
    }
    # BPPARAM <- BiocParallel::bpparam()
    # if (PB) BPPARAM$progressbar <- TRUE
    # out <- BiocParallel::bplapply(X, FUN, ..., BPPARAM = BPPARAM)

  }

  out

}

mapply_ <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE,
                    USE.NAMES = TRUE, BP = FALSE, PB = FALSE) {

  FUN <- match.fun(FUN)

  if (!BP) {

    if (!PB) {

      out <- mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY,
                    USE.NAMES = USE.NAMES)

    } else {

      # not sure how to implement USE.NAMES here, get error sometimes
      dots <- list(...)
      dots.len <- vapply(dots, length, numeric(1))
      dots.len.max <- max(dots.len)
      dots <- lapply(dots, rep, length.out = dots.len.max)
      out <- vector("list", dots.len.max)

      print_pb(0)
      for (i in seq_len(dots.len.max)) {
        dots.i <- mapply(function(dots, i) {
                           if (is.list(dots)) dots[[i]]
                           else dots[i]
                    }, dots, i, SIMPLIFY = FALSE)
        out[[i]] <- do.call(FUN, c(dots.i, MoreArgs))
        update_pb(i, dots.len.max)
      }

      if (SIMPLIFY && length(dots))
        out <- simplify2array(out, higher = (SIMPLIFY == "array"))

    }

  } else {

    if (requireNamespace("BiocParallel", quietly = TRUE)) {
      BPPARAM <- BiocParallel::bpparam()
      if (PB) BPPARAM$progressbar <- TRUE
      out <- BiocParallel::bpmapply(FUN, ..., MoreArgs = MoreArgs,
                                    SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES,
                                    BPPARAM = BPPARAM)
    } else {
      stop("'BiocParallel' is not installed")
    }

  }

  out

}
