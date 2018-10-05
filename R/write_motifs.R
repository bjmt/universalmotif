#' Export motifs in universalmotif format.
#'
#' Write motifs as universalmotif objects to file. For optimal storage of
#' `universalmotif` class motifs, consider using [saveRDS()] and
#' [readRDS()]. The `universalmotif` format will not be documented,
#' as realistically the need to generate these manually/elsewhere should
#' be nonexistent.
#'
#' @param minimal `logical(1)` Only write essential motif information.
#' @param multifreq `logical(1)` Write `multifreq` slot, if present.
#'
#' @return `NULL`, invisibly.
#'
#' @family write_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams write_jaspar
#' @export
write_motifs <- function(motifs, file, minimal = FALSE, multifreq = TRUE) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file), 1, FALSE, "character")
  all_checks <- c(char_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  motifs <- convert_motifs(motifs)
  if (!is.list(motifs)) motifs <- list(motifs)

  out <- paste("# universalmotif version", packageDescription("universalmotif")$Version)

  core <- lapply(motifs, function(x) write_motifs_per_motif(x, minimal, multifreq))

  out <- c(out, do.call(c, core), "")

  writeLines(out, con <- file(file))
  close(con)

  invisible(NULL)

}

write_motifs_per_motif <- function(motif, minimal, multifreq) {

  ## mandatory slots
  out.man <- vector("character", 9)

  out.man[1] <- paste(rep("=", 80), collapse = "")
  out.man[2] <- paste("name:", motif@name)
  out.man[3] <- paste("alphabet:", motif@alphabet)
  out.man[4] <- paste("type:", motif@type)
  out.man[5] <- paste("strand:", motif@strand)
  out.man[6] <- paste("pseudocount:", motif@pseudocount)
  out.man[7] <- paste("bkg:", paste(formatC(motif@bkg, format = "f", digits = 4),
                                    collapse = " "))

  if (!minimal) {
    ## for looks
    out.man[8] <- paste("icscore:", formatC(motif@icscore, format = "f", digits = 2))
    out.man[9] <- paste("consensus:", motif@consensus)
  } else {
    out.man <- out.man[1:7]
  }

  if (minimal) {
    out.opt <- character(0)
    out.ext <- character(0)
  } else {

    out.opt <- vector("character", 7)

    ## optional slots
    if (length(motif@altname) > 0) {
      out.opt[1] <- paste("altname:", motif@altname)
    }
    if (length(motif@organism) > 0) {
      out.opt[2] <- paste("organism:", motif@organism)
    }
    if (length(motif@nsites) > 0) {
      out.opt[3] <- paste("nsites:", motif@nsites)
    }
    if (length(motif@bkgsites) > 0) {
      out.opt[4] <- paste("bkgsites:", motif@bkgsites)
    }
    if (length(motif@pval) > 0) {
      out.opt[5] <- paste("pval:", motif@pval)
    }
    if (length(motif@qval) > 0) {
      out.opt[6] <- paste("qval:", motif@qval)
    }
    if (length(motif@eval) > 0) {
      out.opt[7] <- paste("eval:", motif@eval)
    }

    out.opt <- out.opt[out.opt != ""]

    ## extrainfo
    if (length(motif@extrainfo) > 0) {
      out.ext <- vector("character", length(motif@extrainfo) + 1)
      out.ext[1] <- "extrainfo:"
      for (i in seq_along(motif@extrainfo)) {
        mot.inf <- paste0("> ", names(motif@extrainfo)[i], ": ", motif@extrainfo[i])
        out.ext[i + 1] <- mot.inf
      }
    } else out.ext <- character(0)

  }
  
  alph <- rownames(motif@motif)

  ## motif matrix
  out.mot <- vector("character", nrow(motif@motif) + 1)
  out.mot[1] <- "motif:"
  for (i in seq_len(nrow(motif@motif))) {
    mot.row <- formatC(motif@motif[i, ], format = "f", digits = 6)
    mot.row <- paste(mot.row, collapse = " ")
    out.mot[i + 1] <- paste0(alph[i], "> ", mot.row)
  }

  out.mul <- character(0)
  if (multifreq) {
    ## multifreq
    if (length(motif@multifreq) > 0) {
      out.mul <- "multifreq:"
      for (i in seq_along(motif@multifreq)) {
        out.mul <- c(out.mul, paste(">", names(motif@multifreq)[i]))
        for (j in seq_len(nrow(motif@multifreq[[i]]))) {
          mul.row <- formatC(motif@multifreq[[i]][j, ], format = "f", digits = 6)
          mul.row <- paste(mul.row, collapse = " ")
          out.mul <- c(out.mul, paste(">>", mul.row))
        }
      }
    }
  }

  c(out.man, out.opt, out.ext, out.mot, out.mul)

}
