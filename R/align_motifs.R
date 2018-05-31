#' Align motifs using msaClustalW
#'
#' Align a list of motifs then return motifs of equal length.
#'
#' @param motifs List of motifs
#' @param cluster nj (default or upgma)
#' @param substitutionMatrix iub or clustalw
#' @param ... Other settings for function msaClustalW
#'
#' @return List of motifs.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
align_motifs <- function(motifs, cluster = "default",
                         substitutionMatrix = "iub", ...) {

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
                          left.mat <- matrix(rep(0.25, 4 * left.add), nrow = 4)
                          right.mat <- matrix(rep(0.25, 4 * right.add), nrow = 4)
                          cbind(left.mat, mat, right.mat)
                         }, motif.matrices, matrix.i, SIMPLIFY = FALSE)

  motifs.adjusted <- mapply(function(x, y) {
                              x["motif"] <- y
                              if (x["alphabet"] %in% c("DNA", "RNA")) {
                                x["consensus"] <- apply(x["motif"], 2,
                                                        get_consensus,
                                                        alphabet = x["alphabet"],
                                                        type = "PPM",
                                                        pseudoweight = x["pseudoweight"])
                                colnames(x@motif) <- x["consensus"]
                              } else if (x["alphabet"] == "AA") {
                                x["consensus"] <- apply(x["motif"], 2,
                                                        get_consensusAA,
                                                        type = "PPM",
                                                        pseudoweight = x["pseudoweight"])
                                colnames(x@motif) <- x["consensus"]
                              }
                              return(x)
                            },
                            motifs, matrix.adjusted, SIMPLIFY = FALSE)

  if (length(CLASS_IN) == 1 && length(TYPE_IN) == 1) {
    motifs.adjusted <- convert_type(motifs.adjusted, TYPE_IN)
    motifs.adjusted <- .internal_convert(motifs.adjusted, CLASS_IN)
    motifs.adjusted
  } else motifs.adjusted

}
