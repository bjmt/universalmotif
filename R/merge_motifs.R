#' Merge motifs.
#'
#' DNA/RNA only. The 'extrainfo' slot will be dropped, as well as any
#' non-identical slots. 
#'
#' @param motifs List of motifs.
#' @param method Character.
#' @param min.overlap Numeric.
#' @param min.mean.ic Numeric.
#' @param tryRC Logical.
#' @param relative_entropy Logical.
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return A single motif object.
#'
#' @examples
#' motifs <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#'
#' merged.motif <- merge_motifs(motifs)
#' # compare with:
#' merged.motif2 <- merge_motifs(motifs, method = "motifStack")
#'
#' @seealso \code{\link{compare_motifs}}
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
merge_motifs <- function(motifs, method = "Pearson", min.overlap = 6,
                         min.mean.ic = 0.5, tryRC = TRUE,
                         relative_entropy = FALSE, BPPARAM = SerialParam()) {

  CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  motifs <- convert_motifs(motifs, BPPARAM = BPPARAM)
  motifs <- convert_type(motifs, "PPM", BPPARAM = BPPARAM)

  mot <- merge_mot_list(motifs, tryRC, min.overlap, min.mean.ic, method,
                        relative_entropy)

  mot <- .internal_convert(mot, unique(CLASS_IN), BPPARAM = BPPARAM)
  mot

}

merge_mot_pair <- function(mot1, mot2, weight1, weight2, ic1, ic2, tryRC,
                           min.overlap, min.mean.ic, method, relative_entropy) {

  out <- merge_motifs_internal(mot1, mot2, method, min.overlap, tryRC, ic1, ic2,
                               min.mean.ic, weight1, weight2)
  
  matrix(out[!is.na(out)], nrow = nrow(out))


}

merge_mot_list <- function(motifs, tryRC, min.overlap, min.mean.ic, method,
                           relative_entropy) {

  mot.names <- vapply(motifs, function(x) x["name"], character(1))
  mot.altnames <- do.call(c, sapply(motifs, function(x) x["altname"]))
  mot.families <- do.call(c, sapply(motifs, function(x) x["family"]))
  mot.orgs <- do.call(c, sapply(motifs, function(x) x["organism"]))
  mot.bkgsites <- do.call(c, sapply(motifs, function(x) x["bkgsites"]))
  mot.strands <- vapply(motifs, function(x) x["strand"], character(1))
  mot.extrainfo <- lapply(motifs, function(x) x["extrainfo"])

  mot.mats <- lapply(motifs, function(x) x["motif"])

  alph <- motifs[[1]]["alphabet"]

  mot.ic1 <- .pos_iscscores(motifs[[1]], mot.mats[[1]], relative_entropy)
  mot.ic2 <- .pos_iscscores(motifs[[2]], mot.mats[[2]], relative_entropy)
  new.mat <- merge_mot_pair(mot.mats[[1]], mot.mats[[2]], 1, 1, mot.ic1,
                            mot.ic2, tryRC, min.overlap, min.mean.ic, method,
                            relative_entropy)
  bkg.1 <- motifs[[1]]["bkg"]
  bkg.2 <- motifs[[2]]["bkg"]
  bkg.new <- vapply(seq_along(bkg.1), function(x) mean(c(bkg.1[x], bkg.2[x])),
                    numeric(1))
  nsites.1 <- motifs[[1]]["nsites"]
  nsites.2 <- motifs[[2]]["nsites"]
  nsites.new <- max(c(nsites.1, nsites.2))
  pseudo.1 <- motifs[[1]]["pseudocount"]
  pseudo.2 <- motifs[[2]]["pseudocount"]
  pseudo.new <- mean(c(pseudo.1, pseudo.2))

  mot.new <- create_motif(new.mat, alphabet = alph, pseudocount = pseudo.new,
                          bkg = bkg.new, nsites = nsites.new)

  if (length(motifs) > 2) {
    add.weight <- 2
    for (i in seq(3, length(motifs))) {
      mot.ic <- .pos_iscscores(motifs[[i]], mot.mats[[i]], relative_entropy)
      new.ic <- .pos_iscscores(mot.new, mot.new["motif"], relative_entropy)
      new.mat <- merge_mot_pair(new.mat, mot.mats[[i]], add.weight, 1,
                                new.ic, mot.ic, tryRC, min.overlap,
                                min.mean.ic, method, relative_entropy)
      bkg.1 <- motifs[[i]]["bkg"]
      bkg.2 <- mot.new["bkg"]
      bkg.new <- vapply(seq_along(bkg.1), function(x) mean(c(bkg.1[x], bkg.2[x])),
                      numeric(1))
      nsites.1 <- motifs[[i]]["nsites"]
      nsites.2 <- mot.new["nsites"]
      nsites.new <- max(c(nsites.1, nsites.2))
      pseudo.1 <- motifs[[i]]["pseudocount"]
      pseudo.2 <- mot.new["pseudocount"]
      pseudo.new <- mean(c(pseudo.1, pseudo.2))
      mot.new <- create_motif(new.mat, alphabet = alph, pseudocount = pseudo.new,
                              bkg = bkg.new, nsites = nsites.new)
      add.weight <- add.weight + 1
    }
  }

  new.name <- paste(mot.names, collapse = "/")
  new.altname <- paste(mot.altnames, collapse = "/")
  if (nchar(new.altname) == 0) new.altname <- character(0)
  new.family <- paste(mot.families, collapse = "/")
  if (nchar(new.family) == 0) new.family <- character(0) 
  new.organism <- paste(mot.orgs, collapse = "/")
  if (nchar(new.organism) == 0) new.organism <- character(0)
  if (length(mot.bkgsites) > 1) {
    new.bkgsites <- max(mot.bkgsites)
  } else new.bkgsites <- numeric(0)
  if (length(unique(mot.strands)) > 1) {
    new.strand <- "+-" 
  } else new.strand <- unique(mot.strands)
  new.extrainfo <- do.call(c, mot.extrainfo)

  mot.new["name"] <- new.name
  mot.new["altname"] <- new.altname
  mot.new["family"] <- new.family
  mot.new["organism"] <- new.organism
  mot.new["bkgsites"] <- new.bkgsites
  mot.new["strand"] <- new.strand
  mot.new["extrainfo"] <- new.extrainfo

  mot.new

}
