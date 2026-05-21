#' Utility functions.
#'
#' Utility functions have been split into two categories: those related to
#' motifs ?'utils-motif', and those related to sequences ?'utils-sequence'.
#'
#' @seealso [utils-motif], [utils-sequence]
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @name utilities
NULL

# INTERNAL CONSTANTS -----------------------------------------------------------

DNA_DI <- c("AA", "AC", "AG", "AT",
            "CA", "CC", "CG", "CT",
            "GA", "GC", "GG", "GT",
            "TA", "TC", "TG", "TT")

# AA_STANDARD2 <- sort(AA_STANDARD)
AA_STANDARD2 <- safeExplode("ACDEFGHIKLMNPQRSTVWY")

# TYPE_NULL <- 0L
# TYPE_SYM  <- 1L
# TYPE_ENV  <- 4L
TYPE_LOGI <- 10L
# TYPE_INT  <- 13L
TYPE_NUM  <- 14L
# TYPE_COMP <- 15L
TYPE_CHAR <- 16L
# TYPE_DOT  <- 17L
# TYPE_ANY  <- 18L
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

COMPARE_METRICS <- c("PCC", "EUCL", "SW", "KL", "WEUCL",
                     "ALLR", "BHAT", "HELL", "WPCC",
                     "SEUCL",  "MAN", "ALLR_LL")

# Credit to https://github.com/omarwagih/ggseqlogo/blob/master/R/col_schemes.r
# for the colours.
DNA_COLOURS <- c(A = "#109648", C = "#255C99", G = "#F7B32B", T = "#D62839")
RNA_COLOURS <- c(A = "#109648", C = "#255C99", G = "#F7B32B", U = "#D62839")
AA_COLOURS <- c(G = "#058644", S = "#058644", T = "#058644", Y = "#058644",
  C = "#058644", Q = "#720091", N = "#720091", K = "#0046C5", R = "#0046C5",
  H = "#0046C5", D = "#C5003E", E = "#C5003E", A = "#2E2E2E", V = "#2E2E2E",
  L = "#2E2E2E", I = "#2E2E2E", P = "#2E2E2E", W = "#2E2E2E", F = "#2E2E2E",
  M = "#2E2E2E")

# INTERNAL UTILITIES ----------------------------------------------------------- 

shrink_string <- function(name, maxLen = 5, suffix = "..") {
  if (nchar(name) > maxLen) {
    name <- paste0(substr(name, 1, maxLen), suffix)
  }
  name
}

# Resolve the user-supplied `nthreads` value before handing it to C++.
#
# The documented contract for every user-facing function in this package is
# that `nthreads = 0` means "use all available threads". However,
# RcppThread::parallelFor(..., nthreads = 0) actually sets the global pool to
# zero workers, and then ThreadPool::parallelFor() clamps the worker count
# back up to 1 (RcppThread/ThreadPool.hpp). The work completes but runs
# serially -- silently violating the documented behaviour.
#
# This helper translates `0` into the detected hardware concurrency before
# the value ever reaches C++, so the documented behaviour matches reality
# for every nthreads-aware function in the package.
#
# - `nthreads = 0`     -> parallel::detectCores() (with a floor of 1)
# - `nthreads = NA`    -> 1 (defensive)
# - `nthreads < 0`     -> error (caller is asking for nonsense)
# - otherwise          -> as.integer(nthreads)
resolve_nthreads <- function(nthreads) {
  if (length(nthreads) != 1L)
    stop("'nthreads' must be a single non-negative integer", call. = FALSE)
  if (is.na(nthreads)) return(1L)
  nthreads <- as.integer(nthreads)
  if (is.na(nthreads))
    stop("'nthreads' must be a single non-negative integer", call. = FALSE)
  if (nthreads < 0L)
    stop("'nthreads' cannot be less than 0", call. = FALSE)
  if (nthreads == 0L) {
    n <- parallel::detectCores(logical = TRUE)
    if (is.na(n) || n < 1L) n <- 1L
    return(as.integer(n))
  }
  nthreads
}

# Emit a one-time-per-call hint if a scan_sequences() invocation could be
# served by the faster, multi-threaded scan_sequences2() without losing any
# functionality. Only fires when every argument the user passed maps cleanly
# onto scan_sequences2()'s feature set; advanced features (multifreq,
# gapped motifs, q-values, exhaustive p-values, respect.strand,
# allow.nonfinite, intra-motif overlap removal via the connected-components
# algorithm, non-pvalue threshold types, amino-acid alphabets) suppress it.
#
# Opt out: options(universalmotif.suggest.scan_sequences2 = FALSE).
suggest_scan_sequences2 <- function(threshold.type, use.freq, use.gaps,
                                    allow.nonfinite, no.overlaps,
                                    calc.qvals, respect.strand,
                                    motif_pvalue.method, alphabet,
                                    mot.hasgap) {
  if (!isTRUE(getOption("universalmotif.suggest.scan_sequences2"))) return(invisible())
  if (threshold.type != "pvalue")        return(invisible())
  if (use.freq != 1)                      return(invisible())
  if (any(mot.hasgap) && isTRUE(use.gaps)) return(invisible())
  if (isTRUE(allow.nonfinite))            return(invisible())
  if (isTRUE(no.overlaps))                return(invisible())
  if (isTRUE(calc.qvals))                 return(invisible())
  if (isTRUE(respect.strand))             return(invisible())
  if (motif_pvalue.method != "dynamic")   return(invisible())
  if (!alphabet %in% c("DNA", "RNA"))     return(invisible())

  message(wmsg(
    "Tip: this scan_sequences() call uses only arguments supported by ",
    "scan_sequences2(), a leaner counterpart that parallelises better. ",
    "See ?scan_sequences2. ",
    "Silence with `options(universalmotif.suggest.scan_sequences2 = FALSE)`."
  ))
  invisible()
}

# Emit a one-time-per-call hint if a compare_motifs() invocation could be
# served by compare_motifs2() without losing any functionality. Same
# pattern as suggest_scan_sequences2(). Fires when every argument maps
# cleanly onto compare_motifs2()'s feature set (PCC + sum, default IC
# filters, no multifreq, no report, no normalisation, no db.scores
# lookup, DNA/RNA alphabet).
#
# Opt out: options(universalmotif.suggest.compare_motifs2 = FALSE).
suggest_compare_motifs2 <- function(method, use.freq, use.type,
                                    min.mean.ic, min.position.ic,
                                    relative_entropy, normalise.scores,
                                    score.strat, has.db.scores,
                                    has.output.report, alphabet) {
  if (!isTRUE(getOption("universalmotif.suggest.compare_motifs2")))
    return(invisible())
  if (method != "PCC")                                return(invisible())
  if (!score.strat %in% c("sum", "a.mean"))           return(invisible())
  if (use.freq != 1)                                  return(invisible())
  if (!use.type %in% c("PPM"))                        return(invisible())
  if (!isTRUE(all.equal(min.mean.ic, 0.25)) &&
      !isTRUE(all.equal(min.mean.ic, 0)))             return(invisible())
  if (!isTRUE(all.equal(min.position.ic, 0)))         return(invisible())
  if (isTRUE(relative_entropy))                       return(invisible())
  if (isTRUE(normalise.scores))                       return(invisible())
  if (isTRUE(has.db.scores))                          return(invisible())
  if (isTRUE(has.output.report))                      return(invisible())
  if (!alphabet %in% c("DNA", "RNA"))                 return(invisible())

  message(wmsg(
    "Tip: this compare_motifs() call uses only arguments supported by ",
    "compare_motifs2(), a leaner counterpart that computes empirical-null",
    "p-values and parallelises better. ",
    "See ?compare_motifs2. ",
    "Silence with `options(universalmotif.suggest.compare_motifs2 = FALSE)`."
  ))
  invisible()
}

warn_pseudo <- function(v = 1) {
  # Let's calm down on the warnings maybe...
  if (isTRUE(getOption("pseudocount.warning"))) {
    message(wmsg("Note: Added a pseudocount."))
    # if (v == 1) {
    #   message(wmsg("Note: found -Inf values in motif PWM, adding a pseudocount. ",
    #     "(To turn off this message: `options(pseudocount.warning=FALSE)`.) ",
    #     "Set `allow.nonfinite = TRUE` to prevent this behaviour."))
    # } else if (v == 2) {
    #   message(wmsg("Note: found -Inf values in motif PWM, adding a pseudocount. ",
    #     "(To turn off this message: `options(pseudocount.warning=FALSE)`.)"))
    # } else {
    #   message(wmsg("Note: found -Inf values in motif PWM, adding a pseudocount. ",
    #     "(To turn off this message: `options(pseudocount.warning=FALSE)`.) ", v))
    # }
  }
  invisible()
}

get_nsites <- function(motifs) {
  out <- numeric(length(motifs))
  for (i in seq_along(out)) {
    n <- motifs[[i]]@nsites
    out[i] <- ifelse(length(n) == 1 && n > 1, n, 100)
  }
  out
}

get_bkgs <- function(motifs, use.freq = 1) {

  if (use.freq == 1) {

    out <- lapply(motifs, function(x) x@bkg[seq_len(nrow(x@motif))])

  } else {

    out <- vector("list", length(motifs))
    for (i in seq_along(out)) {
      alph <- rownames(motifs[[i]]@motif)
      alph <- get_klets(alph, use.freq)
      bkg <- motifs[[i]]@bkg[alph]
      if (length(bkg) != nrow(motifs[[i]]@multifreq[[as.character(use.freq)]]))
        stop("Missing higher order background in motif: ", motifs[[i]]@name)
      out[[i]] <- bkg
    }

  }

  out

}

.internal_convert <- function(motifs, class = NULL) {

  if (is.null(class)) {

    CLASS <- class(motifs)
    CLASS_PKG <- attributes(CLASS)$package
    CLASS_IN <- collapse_cpp(c(CLASS_PKG, "-", CLASS))

    CLASS_IN

  } else {

    if (length(class) == 1 && class[1] != "universalmotif-universalmotif" &&
        class[1] != "MotifDb-MotifList") {

      tryCatch(motifs <- convert_motifs(motifs, class),
               error = function(e) message("motifs converted to class 'universalmotif'"))

    } else if (length(class) > 1 || class[1] == "MotifDb-MotifList")
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
