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
#' @importFrom ggseqlogo ggseqlogo geom_logo theme_logo
#' @importFrom ggtree ggtree geom_tiplab geom_tiplab2 groupOTU %<+%
#' @importFrom ggtree geom_tippoint
#' @importFrom BiocParallel bplapply bpmapply bpparam
#' @importFrom stats as.dist hclust runif rnorm chisq.test fisher.test sd fft
#' @importFrom stats p.adjust shapiro.test t.test wilcox.test pnorm quantile
#' @importFrom stats rpois
#' @importFrom utils read.table setTxtProgressBar txtProgressBar menu
#' @importFrom Biostrings width consensusMatrix BString AAString matchPWM
#' @importFrom Biostrings PWM DNAStringSet BStringSet AAStringSet
#' @importFrom Biostrings DNAString RNAString RNAStringSet AA_STANDARD
#' @importFrom Biostrings DNA_BASES RNA_BASES DNA_ALPHABET as.matrix
#' @importFrom Biostrings oligonucleotideTransitions trinucleotideFrequency
#' @importFrom Biostrings dinucleotideFrequency oligonucleotideFrequency
#' @importFrom Biostrings reverseComplement
#' @importFrom ape as.phylo
#' @importFrom TFBSTools PFMatrix toPWM toICM
#' @importFrom TFBSTools TFFMFirst getPosProb
#' @importFrom Rdpack reprompt
#' @importFrom gtools rdirichlet permutations
#' @importFrom Rcpp sourceCpp
#' @importMethodsFrom BiocGenerics cbind rownames colnames ncol nrow rowSums
#' @importMethodsFrom BiocGenerics colMeans rowMeans normalize subset colSums
#' @importMethodsFrom BiocGenerics as.data.frame
#' @useDynLib universalmotif
NULL
