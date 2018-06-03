#' Merge motifs using msa or motifStack
#'
#' DNA/RNA only. The 'extrainfo' slot will be dropped, as well as any
#' non-identical slots. 
#'
#' @param motifs List of motifs.
#' @param newname Name of final merged motif.
#' @param method Character. 'motifStack' or 'msa'.
#' @param tryRC Logical. Only used if method = 'msa'.
#' @param RCstrategy Character. 'motif_dist' or 'motif_simil'.
#' @param bgNoise motifStack::mergeMotifs param.
#' @param cluster msa::msaClustalW param. nj (default) or upgma
#' @param substitutionMatrix msa::msaClustalW param. iub or clustalw
#' @param BPPARAM Param for bplapply. Only used if tryRC = TRUE.
#' @param ... Other settings for function msaClustalW.
#'
#' @return A single motif object.
#'
#' @details
#'    If method = 'msa', the motifs are aligned as consensus strings before
#'    they are averaged as PPMs. For method = 'motifStack', the
#'    motifStack::mergeMotifs function is used. 
#'
#'    If tryRC = FALSE, the motifs
#'    are aligned and merged as-is; if tryRC = TRUE, then one of 
#'    motif_dist or motif_simil is applied to all combinations of the motifs
#'    and their reverse complements. The highest scoring combination of motifs 
#'    is then
#'    used for alignment and merging. Please note that this option
#'    increases computational time significantly with motif number.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
merge_motifs <- function(motifs, newname = "merged motif",
                         method = "msa", tryRC = TRUE,
                         RCstrategy = "motif_dist", bgNoise = NA,
                         cluster = "default",  substitutionMatrix = "iub",
                         BPPARAM = bpparam(), ...) {

  CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  CLASS_IN <- unique(CLASS_IN)

  motifs <- convert_motifs(motifs)

  TYPE_IN <- vapply(motifs, function(x) x["type"], character(1))
  TYPE_IN <- unique(TYPE_IN)

  if (length(CLASS_IN) > 1 || length(TYPE_IN) > 1) {
    message("not all motifs are of the same class and/or type, motifs converted to class 'universalmotif' and type 'PPM'")
  }

  if (method == "msa") {

    alphabet <- vapply(motifs, function(x) x["alphabet"], character(1))
    alphabet <- unique(alphabet)
    if (length(alphabet) != 1) stop("all motifs must share the same alphabet")
    if (!alphabet %in% c("DNA", "RNA")) stop("only DNA and RNA alphabets are allowed")

    motifs <- convert_type(motifs, "PPM", pseudoweight = 0)

    if (tryRC) {
      motifs.rc <- motif_rc(motifs)
      motifs.rc <- lapply(motifs.rc, function(x) {
                          x["name"] <- paste0(x["name"], "-RC")
                          x
                         })
      num_mots <- length(motifs)
      motifs.all <- c(motifs, motifs.rc)
      if (RCstrategy == "motif_simil") {
        motifs.all.simil <- motif_simil(motifs.all)
      } else if (RCstrategy == "motif_dist") {
        motifs.all.simil <- motif_dist(motifs.all)
      } else stop("unrecognized 'RCstrategy'")
      comb.totry <- combn(num_mots * 2, num_mots)
      comb.num <- ncol(comb.totry)
      scores <- bplapply(seq_len(comb.num),
                         function(x) sum(as.dist(motifs.all.simil[comb.totry[, x],
                                                 comb.totry[, x]])),
                         BPPARAM = BPPARAM)
      if (RCstrategy == "motif_simil") {
        motifs.touse <- order(unlist(scores), decreasing = TRUE)[1]
      } else {
        motifs.touse <- order(unlist(scores))[1]
      }
      motifs.touse <- comb.totry[, motifs.touse]
      motifs <- motifs.all[motifs.touse]
    }

    motif.matrices <- lapply(motifs, function(x) x["motif"])

    if (alphabet == "DNA") {
      sequences <- lapply(motif.matrices, function(x) {
                            x <- apply(x, 2, function(x) get_consensus(x, "DNA"))
                            paste(x, collapse = "")
                           })
      sequences <- unlist(sequences)
      aln <- msaClustalW(DNAStringSet(sequences), gapOpening = 999,
                         cluster = cluster, order = "input",
                         substitutionMatrix = substitutionMatrix)
    } else if (alphabet == "RNA") {
      sequences <- lapply(motif.matrices, function(x) {
                            x <- apply(x, 2, function(x) get_consensus(x, "RNA"))
                            paste(x, collapse = "")
                           })
      sequences <- unlist(sequences)
      aln <- msaClustalW(RNAStringSet(sequences), gapOpening = 999,
                         cluster = cluster, order = "input",
                         substitutionMatrix = substitutionMatrix)
    } 

    merged.matrix <- as.matrix(aln)
    split.matrix <- lapply(seq_len(nrow(merged.matrix)),
                           function(x) merged.matrix[x, ])
    matrix.i <- lapply(split.matrix, function(x) x != "-")

    matrix.adjusted <- mapply(function(mat, index) {
                            left.add <- 0
                            for (i in seq_along(index)) {
                              if (isTRUE(index[i])) break
                              left.add <- left.add + 1
                            }
                            right.add <- 0
                            for (j in rev(seq_along(index))) {
                              if (isTRUE(index[j])) break
                              right.add <- right.add + 1
                            }
                            left.mat <- matrix(rep(NA, 4 * left.add), nrow = 4)
                            right.mat <- matrix(rep(NA, 4 * right.add), nrow = 4)
                            cbind(left.mat, mat, right.mat)
                           }, motif.matrices, matrix.i, SIMPLIFY = FALSE)

    matrix.array <- do.call(cbind, matrix.adjusted)
    matrix.array <- array(matrix.array, dim = c(dim(matrix.adjusted[[1]]),
                                                length(matrix.adjusted)))
    final.matrix <- apply(matrix.array, c(1, 2), mean, na.rm = TRUE)

    margs <- list()
    pseudoweight <- unique(vapply(motifs, function(x) x["pseudoweight"], numeric(1)))
    if (length(pseudoweight) == 1) margs <- c(margs, list(pseudoweight = pseudoweight))
    altname <- unique(unlist(sapply(motifs, function(x) x["altname"])))
    if (length(altname) == 1) margs <- c(margs, list(altname = altname))
    nsites <- unique(unlist(sapply(motifs, function(x) x["nsites"])))
    if (length(nsites) == 1) margs <- c(margs, list(nsites = nsites))
    family <- unique(unlist(sapply(motifs, function(x) x["family"])))
    if (length(family) == 1) margs <- c(margs, list(family = family))
    organism <- unique(unlist(sapply(motifs, function(x) x["organism"])))
    if (length(organism) == 1) margs <- c(margs, list(organism = organism))
    bkgsites <- unique(unlist(sapply(motifs, function(x) x["bkgsites"])))
    if (length(bkgsites) == 1) margs <- c(margs, list(bkgsites = bkgsites))
    strand <- unique(vapply(motifs, function(x) x["strand"], character(1)))
    if (length(strand) == 1) margs <- c(margs, list(strand = strand))
    
    motif <- do.call(universalmotif, c(list(motif = final.matrix,
                                            name = newname, alphabet = alphabet),
                                       margs))

  } else if (method == "motifStack") {

    if (CLASS_IN == "motifStack-pcm") {
      motifs <- convert_motifs(motifs, "motifStack-pcm")
    } else if (CLASS_IN == "motifStack-pfm") {
      motifs <- convert_motifs(motifs, "motifStack-pfm")
    } else if (TYPE_IN == "PCM") {
      motifs <- convert_motifs(motifs, "motifStack-pcm")
    } else {
      motifs <- convert_motifs(motifs, "motifStack-pfm")
    }
    motif <- mergeMotifs(motifs, bgNoise = bgNoise)
    motif <- convert_motifs(motif)
    motif["name"] <- newname
  
  } else stop("'method' must be one of 'msa', 'motifStack'")

  if (length(CLASS_IN) == 1 && length(TYPE_IN) == 1) {
    motif <- convert_type(motif, TYPE_IN)
    motif <- .internal_convert(motif, CLASS_IN)
    motif
  } else motif

}