#' Merge motifs using msaClustalW
#'
#' The 'extrainfo' slot will be dropped, as well as any non-identical slots.
#' Gaps are not allowed.
#'
#' @param motifs List of motifs.
#' @param newname Name of final merged motif.
#' @param cluster nj (defaul) or upgma
#' @param substitutionMatrix iub or clustalw
#' @param ... Other settings for function msaClustalW.
#'
#' @return A single motif object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
merge_motifs <- function(motifs, newname = "merged motif", cluster = "default", 
                         substitutionMatrix = "clustalw", ...) {


  CLASS_IN <- vapply(motifs, .internal_convert, character(1))
  CLASS_IN <- unique(CLASS_IN)

  motifs <- convert_motifs(motifs)

  TYPE_IN <- vapply(motifs, function(x) x["type"], character(1))
  TYPE_IN <- unique(TYPE_IN)

  if (length(CLASS_IN) > 1 || length(TYPE_IN) > 1) {
    message("not all motifs are of the same class and/or type, motifs converted to class 'universalmotif' and type 'PPM'")
  }

  alphabet <- vapply(motifs, function(x) x["alphabet"], character(1))
  alphabet <- unique(alphabet)
  if (length(alphabet) != 1) stop("all motifs must share the same alphabet")
  if (!alphabet %in% c("DNA", "RNA")) stop("only DNA and RNA alphabets are allowed")

  motifs <- convert_type(motifs, "PPM", pseudoweight = 0)

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
  altname <- unique(sapply(motifs, function(x) x["altname"]))
  if (length(altname) == 1) margs <- c(margs, list(altname = altname))
  nsites <- unique(unlist(sapply(motifs, function(x) x["nsites"])))
  if (length(nsites) == 1) margs <- c(margs, list(nsites = nsites))
  family <- unique(sapply(motifs, function(x) x["family"]))
  if (length(family) == 1) margs <- c(margs, list(family = family))
  organism <- unique(sapply(motifs, function(x) x["organism"]))
  if (length(organism) == 1) margs <- c(margs, list(organism = organism))
  bkgsites <- unique(unlist(sapply(motifs, function(x) x["bkgsites"])))
  if (length(bkgsites) == 1) margs <- c(margs, list(bkgsites = bkgsites))
  strand <- unique(vapply(motifs, function(x) x["strand"], character(1)))
  if (length(strand) == 1) margs <- c(margs, list(strand = strand))
  
  motif <- do.call(universalmotif, c(list(motif = final.matrix,
                                          name = newname, alphabet = alphabet),
                                     margs))

  if (length(CLASS_IN) == 1 && length(TYPE_IN) == 1) {
    motif <- .internal_convert(motif, CLASS_IN)
    motif <- convert_type(motif, TYPE_IN)
    motif
  } else motif

}
