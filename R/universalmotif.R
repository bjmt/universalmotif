#' universalmotif: Import, Modify and Export Motifs with R
#'
#' @details
#' A collection of utility functions for working with motifs.
#'
#' @docType package
#' @name universalmotif-pkg
#'
#' @import methods
#' @import ggplot2
#' @import ggseqlogo
#' @import ggtree
#' @import BiocParallel
#' @import msa
#' @importFrom stats as.dist hclust runif rnorm
#' @importFrom utils combn read.table
#' @importFrom Biostrings width consensusMatrix BString AAString matchPWM
#' @importFrom Biostrings PWM DNAStringSet BStringSet AAStringSet as.character
#' @importFrom Biostrings DNAString RNAString RNAStringSet AA_STANDARD
#' @importFrom Biostrings DNA_BASES RNA_BASES DNA_ALPHABET as.matrix
#' @importFrom ape as.phylo
#' @importFrom seqLogo seqLogo
#' @importFrom motifStack plotMotifLogo mergeMotifs
#' @importFrom MotIV motifDistances makePWM readDBScores
#' @importFrom TFBSTools PWMSimilarity PFMatrix toPWM toICM
#' @importFrom TFBSTools TFFMFirst TFFMDetail
#' @importFrom Rdpack reprompt
#' @importFrom gtools rdirichlet permutations
NULL
