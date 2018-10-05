#' Run MEME from within R.
#'
#' De novo motif discovery via MEME. For a detailed description of the command,
#' see \url{http://meme-suite.org/doc/meme.html}. For a brief description of
#' the command parameters, run `run_meme()`. Parameters in [run_meme()]
#' which are directly taken from the MEME program are tagged with \[MEME\].
#'
#' @param target.sequences `XStringSet` List of sequences to get motifs from.
#' @param output `character(1)` Name of the output folder. If `NULL`, MEME
#'    output will be deleted.
#' @param overwrite.dir `logical(1)` If `output` is set but already exists,
#'    allow overwritting.
#' @param control.sequences `XStringSet` List of negative sequences. Only
#'    used if `objfun = c("de", "se")`.
#' @param weights `numeric` Vector of numbers between 0 and 1, representing
#'    sequence weights.
#' @param text `logical(1)` \[MEME\]
#' @param brief `numeric(1)` \[MEME\]
#' @param objfun `character(1)` \[MEME\]
#' @param test `character(1)` \[MEME\]
#' @param use_llr `logical(1)` \[MEME\]
#' @param shuf `numeric(1)` \[MEME\]
#' @param hsfrac `numeric(1)` \[MEME\]
#' @param cefrac `numeric(1)` \[MEME\]
#' @param searchsize `numeric(1)` \[MEME\]
#' @param norand `logical(1)` \[MEME\]
#' @param csites `numeric(1)` \[MEME\]
#' @param seed `numeric(1)` \[MEME\]
#' @param alph `character(1)` \[MEME\]
#' @param revcomp `logical(1)` \[MEME\]
#' @param pal `logical(1)` \[MEME\]
#' @param mod `character(1)` \[MEME\]
#' @param nmotifs `numeric(1)` \[MEME\]
#' @param evt `numeric(1)` \[MEME\]
#' @param nsites `numeric(1)` \[MEME\]
#' @param minsites `numeric(1)` \[MEME\]
#' @param maxsites `numeric(1)` \[MEME\]
#' @param wnsites `numeric(1)` \[MEME\]
#' @param w `numeric(1)` \[MEME\]
#' @param minw `numeric(1)` \[MEME\]
#' @param maxw `numeric(1)` \[MEME\]
#' @param allw `numeric(1)` \[MEME\]
#' @param nomatrim `logical(1)` \[MEME\]
#' @param wg `numeric(1)` \[MEME\]
#' @param ws `numeric(1)` \[MEME\]
#' @param noendgaps `logical(1)` \[MEME\]
#' @param bfile `character(1)` \[MEME\]
#' @param markov_order `numeric(1)` \[MEME\]
#' @param psp `character(1)` \[MEME\]
#' @param maxiter `numeric(1)` \[MEME\]
#' @param distance `numeric(1)` \[MEME\]
#' @param prior `character(1)` \[MEME\]
#' @param b `numeric(1)` \[MEME\]
#' @param plib `character(1)` \[MEME\]
#' @param spfuzz `numeric(1)` \[MEME\]
#' @param spmap `character(1)` \[MEME\]
#' @param cons `character(1)` \[MEME\]
#' @param p `numeric(1)` \[MEME\]
#' @param maxsize `numeric(1)` \[MEME\]
#' @param maxtime `numeric(1)` \[MEME\]
#' @param wd `character(1)` Working directory to run MEME in.
#' @param logfile `character(1)` File to dump MEME stderr. If `NULL`, no logs
#'    will be saved.
#' @param readsites `logical(1)` Read sites from MEME output (from [read_meme()]).
#' @param echo `logical(1)` Dump MEME output to console.
#' @param verbose `numeric(1)` Set `verbose = 0` to quiet [run_meme()].
#' @param timeout `numeric(1)` Stop MEME program past `timeout` (seconds). See
#'    [processx::run()].
#' @param bin `character(1)` Location of MEME binary. Alternatively, set this
#'    via `options(meme.bin = '/path/to/meme/bin')`.
#'
#' @return `list` The output file is read with [read_meme()].
#'
#' @examples
#' \dontrun{
#' ## To check that you are properly linking to the binary:
#' run_meme()
#' }
#'
#' @references
#'    \insertRef{meme3}{universalmotif}
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @seealso [read_meme()], [create_sequences()], [shuffle_sequences()],
#'    [processx::run()]
#' @export
run_meme <- function(target.sequences, output = NULL,
                     overwrite.dir = FALSE, control.sequences = NULL,
                     weights = NULL, text = FALSE, brief = 1000, objfun = "classic",
                     test = NULL, use_llr = FALSE, shuf = 2, hsfrac = NULL,
                     cefrac = NULL, searchsize = NULL, norand = FALSE,
                     csites = 1000, seed = 0, alph = NULL, revcomp = FALSE,
                     pal = FALSE, mod = "zoops", nmotifs = 3, evt = NULL,
                     nsites = NULL, minsites = NULL,
                     maxsites = NULL, wnsites = 0.8, w = NULL, minw = 8,
                     maxw = 50, allw = NULL, nomatrim = FALSE, wg = 11,
                     ws = 1, noendgaps = FALSE, bfile = NULL,
                     markov_order = 0, psp = NULL, maxiter = 50,
                     distance = 0.001, prior = NULL, b = NULL, plib = NULL,
                     spfuzz = NULL, spmap = NULL, cons = NULL, p = NULL,
                     maxsize = NULL, maxtime = NULL,
                     wd = getwd(), logfile = paste0(wd, "/memerun.log"),
                     readsites = TRUE, echo = FALSE,
                     verbose = 1, timeout = Inf, bin = getOption("meme.bin")) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(output = args$output,
                                      objfun = args$objfun,
                                      test = args$test,
                                      alph = args$alph,
                                      mod = args$mod,
                                      bfile = args$bfile,
                                      psp = args$psp,
                                      plib = args$lib,
                                      spmap = args$spmap,
                                      cons = args$cons,
                                      prior = args$prior,
                                      wd = args$wd,
                                      logfile = args$logfile,
                                      bin = args$bin),
                                 numeric(), c(TRUE, FALSE, rep(TRUE, 9),
                                              FALSE, TRUE, TRUE),
                                 "character")
  num_check <- check_fun_params(list(brief = args$brief,
                                     shuf = args$shuf,
                                     hsfrac = args$hsfrac,
                                     cefrac = args$cefrac,
                                     searchsize = args$searchsize,
                                     csites = args$csites,
                                     seed = args$seed,
                                     nmotifs = args$nmotifs,
                                     evt = args$evt,
                                     nsites = args$nsites,
                                     minsites = args$minsites,
                                     maxsites = args$maxsites,
                                     wnsites = args$wnsites,
                                     w = args$w,
                                     minw = args$minw,
                                     maxw = args$maxw,
                                     allw = args$allw,
                                     wg = args$wg,
                                     ws = args$ws,
                                     markov_order = args$markov_order,
                                     maxiter = args$maxiter,
                                     distance = args$distance,
                                     b = args$b,
                                     spfuzz = args$spfuzz,
                                     p = args$p,
                                     maxsize = args$maxsize,
                                     maxtime = args$maxtime,
                                     verbose = args$verbose,
                                     timeout = args$timeout),
                                numeric(), c(rep(TRUE, 27), FALSE, FALSE),
                                             "numeric")
  logi_check <- check_fun_params(list(overwrite.dir = args$overwrite.dir,
                                      text = args$text,
                                      use_llr = args$use_llr,
                                      norand = args$norand,
                                      revcomp = args$revcomp,
                                      pal = args$pal,
                                      nomatrim = args$nomatrim,
                                      noendgaps = args$noendgaps,
                                      echo = args$echo),
                                 numeric(), logical(), "logical")
  all_checks <- c(char_check, num_check, logi_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  args <- as.list(environment())
  v <- verbose
  os.check <- .Platform$OS.type
  if (os.check == "windows") stop("MEME is not available on windows")

  if (is.null(bin)) stop("please specify the location of the MEME binary")
  meme.version <- tryCatch(processx::run(bin, "-version", error_on_status = FALSE),
                           error = function(e) stop("could not find the MEME binary"))
  meme.version <- sub("\n", "", meme.version$stdout, fixed = TRUE)
  meme.major <- as.numeric(strsplit(meme.version, ".", fixed = TRUE)[[1]][1])
  if (!meme.major %in% 4:5)
    warning("'run_meme' has been optimized for MEME versions 4-5",
            immediate. = TRUE)

  if (missing(target.sequences)) {
    help <- processx::run(bin, "-h", error_on_status = FALSE)$stderr
    cat("MEME version ", meme.version, "\n", bin, "\n\n", sep = "")
    cat(help)
    return(invisible(NULL))
  }

  if (!is(target.sequences, "XStringSet"))
    stop("'sequences' must be an 'XStringSet' object")

  if (v>0) cat(paste0("Using MEME version ", meme.version, "\n"))
  if (is.null(output)) {
    if (v>0) cat("No output folder specified, output will be deleted.\n")
    delete.ouput <- TRUE
    output <- tempdir()
  } else {
    if (v>0) cat(paste0("Output folder: ", output, "\n"))
    delete.ouput <- FALSE
    if (dir.exists(output)) {
      if (!overwrite.dir) {
        stop("Output folder exists but 'overwrite.dir' is set to FALSE")
      } else {
        if (v>0) cat("NOTE: output folder already exists, will be overwritten\n")
      }
    }
  }

  if (v>0) {

    if (objfun == "classic") cat("Search mode: Classic\n")
    else if (objfun == "de") cat("Search mode: Differential Enrichment\n")
    else if (objfun == "se") cat("Search mode: Selective Enrichment\n")
    else if (objfun == "cd") cat("Search mode: Central Distance\n")
    else if (objfun == "ce") cat("Search mode: Central Enrichment\n")
    else if (objfun == "nc") cat("Search mode: Numerically Correct\n")

    if (objfun %in% c("de", "se")) {
      if (test == "mhg" || is.null(test)) cat("Test: Multiple Hypergeometric\n")
      else if (test == "mbn") cat("Test: Multiple Binomial\n")
      else if (test == "mrs") cat("Test: Multiple Rank-Sum\n")
    }

    if (mod == "zoops" || is.null(mod))
      cat("Model: Zero or One Occurrence Per Sequence\n")
    else if (mod == "oops") cat("Model: One Occurrence Per Sequence\n")
    else if (mod == "anr") cat("Model: Any Number of Repetitions\n")

    if (is.null(nmotifs) || nmotifs == 1)
      cat("Looking for 1 motif\n")
    else
      cat("Looking for", nmotifs, "motifs\n")

  }

  if (is.null(alph)) {
    if (is(target.sequences, "DNAStringSet")) alph.arg <- "-dna"
    else if (is(target.sequences, "RNAStringSet")) alph.arg <- "-rna"
    else if (is(target.sequences, "AAStringSet")) alph.arg <- "-protein"
    else if (is(target.sequences, "BStringSet"))
      stop("for custom alphabets, please pass an alphabet file to 'alph'\n",
           "  (http://meme-suite.org/doc/alphabet-format.html)")
  } else alph.arg <- "-file"

  if (is.null(names(target.sequences)))
    names(target.sequences) <- as.character(seq_len(length(target.sequences)))
  dataset <- paste0(wd, "/target.temp.fasta")
  if (is.null(weights))
    writeXStringSet(target.sequences, dataset)
  else {
    if (length(target.sequences != length(weights)))
      stop("length of 'weights' must match length of 'target.sequences'")
    weights <- paste(weights, collapse = " ")
    weights <- paste(">WEIGHTS", weights, "\n")
    cat(weights, file = dataset)
    writeXStringSet(target.sequences, dataset, append = TRUE)
  }

  if (delete.ouput) to.delete <- c(dataset, output) else to.delete <- dataset

  if (!is.null(control.sequences)) {
    if (!is(target.sequences, class(control.sequences)))
      stop(paste0("'target.sequences' has class ", class(target.sequences),
                  " whereas 'control.sequences' has class ",
                  class(control.sequences)))
    bkg <- paste0(wd, "/control.temp.fasta")
    writeXStringSet(control.sequences, bkg)
    to.delete <- c(to.delete, bkg)
    has.bkg <- TRUE
  } else has.bkg <- FALSE

  on.exit(unlink(to.delete, TRUE))

  meme.args <- dataset
  if (overwrite.dir) output.arg <- "-oc" else output.arg <- "-o"
  meme.args <- c(meme.args, output.arg, output)

  if (alph.arg != "file") meme.args <- c(meme.args, alph.arg)
  else meme.args <- c(meme.args, "-alph", alph)

  if (has.bkg) meme.args <- c(meme.args, "-neg", bkg)

  if (is.null(control.sequences)) meme.args <- c(meme.args, "-shuf", shuf)

  meme.args <- c(meme.args, "-objfun", objfun)

  if (text) meme.args <- c(meme.args, "-text")
  if (use_llr) meme.args <- c(meme.args, "-use_llr")
  if (norand) meme.args <- c(meme.args, "-norand")
  if (revcomp) meme.args <- c(meme.args, "-revcomp")
  if (pal) meme.args <- c(meme.args, "-pal")
  if (nomatrim) meme.args <- c(meme.args, "-nomatrim")
  if (noendgaps) meme.args <- c(meme.args, "-noendgaps")

  if (!is.null(mod)) meme.args <- c(meme.args, "-mod", mod)
  if (!is.null(minw)) meme.args <- c(meme.args, "-minw", minw)
  if (!is.null(maxw)) meme.args <- c(meme.args, "-maxw", maxw)
  if (!is.null(wg)) meme.args <- c(meme.args, "-wg", wg)
  if (!is.null(ws)) meme.args <- c(meme.args, "-ws", ws)
  if (!is.null(markov_order)) meme.args <- c(meme.args, "-markov_order", markov_order)
  if (!is.null(maxiter)) meme.args <- c(meme.args, "-maxiter", maxiter)
  if (!is.null(distance)) meme.args <- c(meme.args, "-distance", distance)
  if (!is.null(brief)) meme.args <- c(meme.args, "-brief", brief)
  if (!is.null(shuf)) meme.args <- c(meme.args, "-shuf", shuf)
  if (!is.null(csites)) meme.args <- c(meme.args, "-csites", csites)
  if (!is.null(nmotifs)) meme.args <- c(meme.args, "-nmotifs", nmotifs)
  if (!is.null(seed)) meme.args <- c(meme.args, "-seed", seed)
  if (!is.null(test)) meme.args <- c(meme.args, "-test", test)
  if (!is.null(hsfrac)) meme.args <- c(meme.args, "-hsfrac", hsfrac)
  if (!is.null(cefrac)) meme.args <- c(meme.args, "-cefrac", cefrac)
  if (!is.null(searchsize)) meme.args <- c(meme.args, "-searchsize", searchsize)
  if (!is.null(evt)) meme.args <- c(meme.args, "-evt", evt)
  if (!is.null(nsites)) meme.args <- c(meme.args, "-nsites", nsites)
  if (!is.null(minsites)) meme.args <- c(meme.args, "-minsites", minsites)
  if (!is.null(wnsites)) meme.args <- c(meme.args, "-wnsites", wnsites)
  if (!is.null(w)) meme.args <- c(meme.args, "-w", w)
  if (!is.null(allw)) meme.args <- c(meme.args, "-allw", allw)
  if (!is.null(bfile)) meme.args <- c(meme.args, "-bfile", bfile)
  if (!is.null(psp)) meme.args <- c(meme.args, "-psp", psp)
  if (!is.null(prior)) meme.args <- c(meme.args, "-prior", prior)
  if (!is.null(b)) meme.args <- c(meme.args, "-b", b)
  if (!is.null(plib)) meme.args <- c(meme.args, "-plib", plib)
  if (!is.null(spfuzz)) meme.args <- c(meme.args, "-spfuzz", spfuzz)
  if (!is.null(spmap)) meme.args <- c(meme.args, "-spmap", spmap)
  if (!is.null(cons)) meme.args <- c(meme.args, "-cons", cons)
  if (!is.null(p)) meme.args <- c(meme.args, "-p", p)
  if (!is.null(maxsize)) meme.args <- c(meme.args, "-maxsize", maxsize)

  if (v>0) cat(paste0("\n *** Starting MEME ***\n\n", pdate(), "\n"))

  t.start <- Sys.time()
  run.res <- processx::run(bin, meme.args, error_on_status = FALSE,
                           wd = wd, timeout = timeout,
                           echo_cmd = ifelse(v>0, TRUE, FALSE),
                           echo = ifelse(echo, TRUE, FALSE),
                           stderr_line_callback = ifelse(echo||v==0, NULL, meme_cb))
  t.stop <- Sys.time()

  if (!is.null(logfile)) cat(run.res$stderr, file = logfile)
  if (run.res$timeout) stop("MEME process timed out")
  if (run.res$status != 0) {
    if (!is.null(logfile)) {
      logs <- readLines(con <- file(logfile)); close(con)
      stop("MEME had a non-zero exit, dumping last line of logfile:\n  ",
           logs[length(logs)])
    } else
      stop("MEME had a non-zero exit, dumping stderr\n", run.res$stderr)
  }

  motifs <- read_meme(paste0(output, "/meme.txt"), readsites = readsites)

  t.diff <- format(difftime(t.stop, t.start))
  if (v>0) cat(paste(pdate(), "\n\n *** Run over ***\n\nTotal runtime:", t.diff, "\n"))

  motifs

}

pdate <- function() paste("[", date(), "]")

meme_cb <- function(line, proc) {
  if (grepl("motif=", line)) {
    motif.num <- strsplit(line, "=")[[1]][2]
    cat(pdate(), "\n")
    cat("Generating motif", motif.num, "\n")
  }
}