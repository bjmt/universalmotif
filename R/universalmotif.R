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
#' @import ggtree
#' @import ggseqlogo
#' @importFrom BiocParallel bplapply bpmapply bpparam
#' @importFrom stats as.dist hclust runif rnorm chisq.test fisher.test sd fft
#' @importFrom stats p.adjust shapiro.test t.test wilcox.test pnorm quantile
#' @importFrom stats rpois
#' @importFrom utils read.table setTxtProgressBar txtProgressBar
#' @importFrom Biostrings width consensusMatrix BString AAString matchPWM
#' @importFrom Biostrings PWM DNAStringSet BStringSet AAStringSet
#' @importFrom Biostrings DNAString RNAString RNAStringSet AA_STANDARD
#' @importFrom Biostrings DNA_BASES RNA_BASES DNA_ALPHABET as.matrix
#' @importFrom Biostrings oligonucleotideTransitions trinucleotideFrequency
#' @importFrom Biostrings dinucleotideFrequency oligonucleotideFrequency
#' @importFrom Biostrings reverseComplement
#' @importFrom ape as.phylo
#' @importFrom seqLogo seqLogo
#' @importFrom motifStack mergeMotifs
#' @importFrom MotIV makePWM
#' @importFrom TFBSTools PFMatrix toPWM toICM
#' @importFrom TFBSTools TFFMFirst getPosProb
#' @importFrom Rdpack reprompt
#' @importFrom gtools rdirichlet permutations
#' @importFrom PWMEnrich PFMtoPWM
#' @importFrom Rcpp sourceCpp
#' @importClassesFrom MotifDb MotifList
#' @importMethodsFrom BiocGenerics cbind rownames colnames ncol nrow rowSums
#' @importMethodsFrom BiocGenerics colMeans rowMeans normalize subset colSums
#' @importMethodsFrom BiocGenerics as.data.frame
#' @useDynLib universalmotif
NULL
