## ---------------------------------------------------------------------------
## Shared HTML report generator for compare_motifs() and compare_motifs2().
##
## Both functions normalise their results into a common `matches` data.frame
## (query / query.i / target / target.i / offset / strand + statistic columns)
## and hand it to compare_report(), which writes a temporary R Markdown file
## (a self-contained html_document with a sortable paged results table and a
## per-match query-vs-target logo alignment) and renders it with rmarkdown.
## All functions here are internal.
## ---------------------------------------------------------------------------

## Orchestrator. `view.fun` is the name of the logo function emitted in the
## figure chunks ("view_motifs" or "view_motifs2"); `view.args` is a named
## list of scalar arguments splatted into that call; `summary` is a named
## list of label -> value bullets for the header.
compare_report <- function(call, matches, motifs, view.fun, view.args,
                           summary, output, max.print) {

  output <- .cr_check_output(output)

  if (!is.data.frame(matches) || nrow(matches) == 0L) {
    warning("No matches to report", call. = FALSE)
    return(FALSE)
  }
  max.print <- as.integer(max.print)
  if (max.print < nrow(matches))
    matches <- matches[seq_len(max.print), , drop = FALSE]

  ## Only serialise the motifs actually shown (a query-vs-DB run may hold
  ## thousands), remapping the pair indices into the subset.
  needed <- sort(unique(c(matches$query.i, matches$target.i)))
  remap  <- setNames(seq_along(needed), as.character(needed))
  matches$query.i.local  <- unname(remap[as.character(matches$query.i)])
  matches$target.i.local <- unname(remap[as.character(matches$target.i)])
  motifs.sub <- motifs[needed]

  table_df <- .cr_build_table(matches)

  ## Render in a private directory so html_document's self_contained step
  ## base64-embeds the figures (it does not when output_dir differs from the
  ## Rmd's location), then copy the single self-contained file to `output`.
  tmpd <- tempfile("ureport")
  dir.create(tmpd)
  on.exit(unlink(tmpd, recursive = TRUE), add = TRUE)

  rds <- file.path(tmpd, "payload.rds")
  saveRDS(list(table = table_df, motifs = motifs.sub), rds)

  rmd <- .cr_build_rmd(call, summary, matches, view.fun, view.args, rds)
  f <- file.path(tmpd, "report.Rmd")
  writeLines(rmd, f)

  rmarkdown::render(f, output_file = "report.html", quiet = TRUE)
  if (!file.copy(file.path(tmpd, "report.html"), output, overwrite = TRUE))
    stop("could not write the report to '", output, "'", call. = FALSE)
  TRUE
}

## Validate the output path: single non-empty .html/.htm filename in an
## existing, writable directory.
.cr_check_output <- function(output) {
  if (!is.character(output) || length(output) != 1L || is.na(output) ||
      !nzchar(output))
    stop("`output.report` must be a single non-empty filename", call. = FALSE)
  ext <- tolower(tools::file_ext(output))
  if (!ext %in% c("html", "htm"))
    stop("`output.report` must end in .html or .htm (got '", output, "')",
         call. = FALSE)
  dir <- dirname(output)
  if (!dir.exists(dir))
    stop("`output.report` directory does not exist: ", dir, call. = FALSE)
  if (file.access(dir, mode = 2L) != 0L)
    stop("`output.report` directory is not writable: ", dir, call. = FALSE)
  output
}

## Display table: drop internal index columns and round statistic columns.
.cr_build_table <- function(matches) {
  drop <- c("query.i", "target.i", "query.i.local", "target.i.local")
  tab <- matches[, setdiff(names(matches), drop), drop = FALSE]
  for (nm in intersect(c("score", "logPval", "Pval", "Eval", "pvalue",
                         "qvalue"), names(tab)))
    tab[[nm]] <- signif(as.numeric(tab[[nm]]), 3L)
  rownames(tab) <- NULL
  tab
}

## Escape markdown/HTML metacharacters in a string destined for a header.
.cr_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  gsub("([][`*_{}()#+.!<>|\\-])", "\\\\\\1", x, perl = TRUE)
}

## Assemble the R Markdown lines. `rds` is the path the setup chunk reads.
.cr_build_rmd <- function(call, summary, matches, view.fun, view.args, rds) {

  arg_str <- if (length(view.args))
    paste(names(view.args),
          vapply(view.args, function(v) deparse(v, width.cutoff = 500L),
                 character(1)),
          sep = " = ", collapse = ", ")
  else ""

  ## Rendered as inline code, so it must not be markdown-escaped. trimws()
  ## strips the indentation deparse() adds to wrapped continuation lines,
  ## which otherwise shows as a gap in the joined one-line call.
  call_str <- paste(trimws(deparse(call)), collapse = " ")

  hdr <- c(
    "---",
    'title: "universalmotif comparison report"',
    paste0('date: "', format(Sys.time()), '"'),
    "output:",
    "  html_document:",
    "    df_print: paged",
    "    toc: true",
    "    toc_float: true",
    "    self_contained: true",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "suppressPackageStartupMessages(library(universalmotif))",
    paste0('.report <- readRDS("', rds, '")'),
    "```",
    "",
    "## Summary",
    "",
    paste0("- **Call:** `", call_str, "`"),
    unlist(lapply(names(summary),
                  function(nm) paste0("- **", nm, ":** ", summary[[nm]]))),
    "",
    "## Results",
    "",
    "```{r results-table}",
    ".report$table",
    "```",
    "",
    "## Match alignments",
    ""
  )

  body <- unlist(lapply(seq_len(nrow(matches)), function(k) {
    q   <- .cr_escape(matches$query[k])
    tg  <- .cr_escape(matches$target[k])
    rcf <- if (!is.na(matches$strand[k]) && matches$strand[k] == "-")
      " [RC]" else ""
    qi  <- matches$query.i.local[k]
    ti  <- matches$target.i.local[k]
    fig <- paste0(view.fun, "(.report$motifs[c(", qi, ", ", ti, ")]",
                  if (nzchar(arg_str)) paste0(", ", arg_str) else "", ")")
    c(
      paste0("### ", k, ". ", q, " vs ", tg, rcf),
      "",
      paste0("*score: ", signif(as.numeric(matches$score[k]), 3L),
             ", offset: ", matches$offset[k],
             ", strand: ", matches$strand[k], "*"),
      "",
      paste0("```{r logo-", k, ", error=TRUE, fig.width=7, fig.height=2.6, ",
             "dpi=100, fig.alt='motif alignment'}"),
      fig,
      "```",
      ""
    )
  }))

  c(hdr, body)
}

## v1 helper: compare_motifs() does not return alignment offset/strand, so
## recompute them per top pair with the same C++ aligner view_motifs() uses
## (view_motifs_prep). Passing the pair query-first makes the query the
## anchor (the IC re-sort in view_motifs() happens in R, before this call),
## so the query-relative offset is target-pad minus query-pad. This is a
## fresh alignment, independent of the C++ path that produced `score`, and
## can rarely disagree on ties; that is acceptable for a display report.
.cr_offsets_via_prep <- function(mot.mats, mot.bkgs, mot.nsites, alph, pairs,
                                 prep.args) {
  n <- nrow(pairs)
  offs <- integer(n)
  strs <- character(n)
  for (k in seq_len(n)) {
    qi <- pairs$qi[k]; ti <- pairs$ti[k]
    res <- tryCatch(
      view_motifs_prep(mot.mats[c(qi, ti)], prep.args$method, prep.args$tryRC,
                       prep.args$min.overlap, prep.args$min.mean.ic,
                       prep.args$min.position.ic, mot.bkgs[c(qi, ti)],
                       prep.args$relative_entropy, prep.args$normalise.scores,
                       alph, mot.nsites[c(qi, ti)], prep.args$score.strat),
      error = function(e) NULL)
    if (is.null(res) || length(res$offsets) < 2L) {
      offs[k] <- NA_integer_
      strs[k] <- NA_character_
    } else {
      offs[k] <- as.integer(res$offsets[2] - res$offsets[1])
      strs[k] <- if (isTRUE(res$motIsRC[1])) "-" else "+"
    }
  }
  data.frame(offset = offs, strand = strs, stringsAsFactors = FALSE)
}

## Summary bullets: query/target motif counts and their mean widths.
## `query.i` indexes the query motifs within `motifs`; the targets are the
## remaining motifs searched against (queries excluded).
.cr_motif_stats <- function(motifs, query.i) {
  w  <- vapply(motifs, function(m) ncol(m@motif), integer(1))
  qi <- unique(query.i)
  ti <- setdiff(seq_along(w), qi)
  meanw <- function(i) if (length(i)) sprintf("%.1f", mean(w[i])) else "NA"
  list(
    "Query motifs"  = sprintf("%d (mean width %s)", length(qi), meanw(qi)),
    "Target motifs" = sprintf("%d (mean width %s)", length(ti), meanw(ti))
  )
}

## Adapter: build the normalised matches table from compare_motifs() (v1)
## output and render. v1 lacks offset/strand, so recompute them.
compare_report_from_v1 <- function(fun.call, comparisons, motifs, mot.names,
                                    mot.mats, mot.bkgs, mot.nsites, alph,
                                    args, query.i, output, max.print) {
  np  <- min(as.integer(max.print), nrow(comparisons))
  top <- as.data.frame(comparisons[seq_len(np), , drop = FALSE])
  off <- .cr_offsets_via_prep(
    mot.mats, mot.bkgs, mot.nsites, alph,
    data.frame(qi = top$subject.i, ti = top$target.i),
    list(method = args$method, tryRC = args$tryRC,
         min.overlap = args$min.overlap, min.mean.ic = args$min.mean.ic,
         min.position.ic = args$min.position.ic,
         relative_entropy = args$relative_entropy,
         normalise.scores = args$normalise.scores,
         score.strat = args$score.strat))
  matches <- data.frame(
    query  = top$subject, query.i  = top$subject.i,
    target = top$target,  target.i = top$target.i,
    offset = off$offset,  strand   = off$strand,
    score  = top$score,   logPval  = top$logPval,
    Pval   = top$Pval,    Eval     = top$Eval,
    stringsAsFactors = FALSE)
  motifs <- Map(function(m, nm) { m@name <- nm; m }, motifs, mot.names)
  compare_report(
    fun.call, matches, motifs, "view_motifs",
    view.args = list(method = args$method, tryRC = args$tryRC,
                     min.overlap = args$min.overlap,
                     min.mean.ic = args$min.mean.ic,
                     relative_entropy = args$relative_entropy,
                     normalise.scores = args$normalise.scores,
                     min.position.ic = args$min.position.ic,
                     score.strat = args$score.strat,
                     ## Display as information content (bits), matching the
                     ## v2 report; args$use.type is the comparison space, not
                     ## the display type.
                     use.type = "ICM"),
    summary = c(.cr_motif_stats(motifs, query.i),
                list(Method = args$method,
                     "Reverse complement" = args$tryRC,
                     "Max P-value" = args$max.p, "Max E-value" = args$max.e,
                     "Total matches" = nrow(comparisons))),
    output = output, max.print = max.print)
}

## Adapter: build the matches table from compare_motifs2() (v2) output and
## render. v2 already carries offset/strand from the scored alignment.
compare_report_from_v2 <- function(fun.call, long, motifs, mot.names,
                                    args, query.i, output, max.print) {
  np  <- min(as.integer(max.print), nrow(long))
  top <- long[seq_len(np), , drop = FALSE]
  matches <- data.frame(
    query  = top$subject, query.i  = top$subject.i,
    target = top$target,  target.i = top$target.i,
    offset = top$offset,  strand   = top$strand,
    overlap = top$overlap, score   = top$score,
    pvalue = top$pvalue,  qvalue   = top$qvalue,
    stringsAsFactors = FALSE)
  motifs <- Map(function(m, nm) { m@name <- nm; m }, motifs, mot.names)
  compare_report(
    fun.call, matches, motifs, "view_motifs2",
    view.args = list(use.type = "ICM", tryRC = args$RC,
                     min.overlap = args$min.overlap),
    summary = c(.cr_motif_stats(motifs, query.i),
                list("Null model" = args$null,
                     "Reverse complement" = args$RC,
                     "q-value cutoff" = args$qvalue,
                     "Total matches" = nrow(long))),
    output = output, max.print = max.print)
}
