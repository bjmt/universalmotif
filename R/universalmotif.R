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
#' @import msa
#' @importFrom BiocParallel bplapply bpmapply bpparam
#' @importFrom stats as.dist hclust runif rnorm
#' @importFrom utils read.table
#' @importFrom Biostrings width consensusMatrix BString AAString matchPWM
#' @importFrom Biostrings PWM DNAStringSet BStringSet AAStringSet
#' @importFrom Biostrings DNAString RNAString RNAStringSet AA_STANDARD
#' @importFrom Biostrings DNA_BASES RNA_BASES DNA_ALPHABET as.matrix
#' @importFrom Biostrings oligonucleotideTransitions trinucleotideFrequency
#' @importFrom Biostrings dinucleotideFrequency oligonucleotideFrequency
#' @importFrom Biostrings reverseComplement
#' @importFrom ape as.phylo
#' @importFrom seqLogo seqLogo
#' @importFrom motifStack plotMotifLogo mergeMotifs
#' @importFrom MotIV motifDistances makePWM readDBScores
#' @importFrom TFBSTools PWMSimilarity PFMatrix toPWM toICM
#' @importFrom TFBSTools TFFMFirst TFFMDetail
#' @importFrom Rdpack reprompt
#' @importFrom gtools rdirichlet permutations
#' @importFrom Logolas logomaker
#' @importFrom PWMEnrich PFMtoPWM
#' @importFrom Rcpp sourceCpp
#' @useDynLib universalmotif
#' @importClassesFrom MotifDb MotifList
NULL
