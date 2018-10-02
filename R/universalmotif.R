#' universalmotif: Import, Modify and Export Motifs with R
#'
#' @details
#' A collection of utility functions for working with motifs.
#'
#' @docType package
#' @name universalmotif-pkg
#'
#' @importFrom methods getClass is new slot "slot<-" validObject
#' @importMethodsFrom methods initialize show
#' @importFrom ggplot2 aes element_blank facet_wrap ggplot theme ylab
#' @importFrom ggseqlogo ggseqlogo geom_logo theme_logo
#' @importFrom ggtree ggtree geom_tiplab geom_tiplab2 groupOTU %<+%
#' @importFrom ggtree geom_tippoint
#' @importFrom stats as.dist hclust runif rnorm chisq.test fisher.test sd fft
#' @importFrom stats p.adjust shapiro.test t.test wilcox.test pnorm quantile
#' @importFrom stats rpois
#' @importFrom utils read.table menu packageDescription
#' @importFrom Biostrings width consensusMatrix BString AAString matchPWM
#' @importFrom Biostrings PWM DNAStringSet BStringSet AAStringSet
#' @importFrom Biostrings DNAString RNAString RNAStringSet AA_STANDARD
#' @importFrom Biostrings DNA_BASES RNA_BASES DNA_ALPHABET as.matrix
#' @importFrom Biostrings oligonucleotideTransitions trinucleotideFrequency
#' @importFrom Biostrings dinucleotideFrequency oligonucleotideFrequency
#' @importFrom Biostrings reverseComplement writeXStringSet
#' @importFrom ape as.phylo
#' @importFrom Rdpack reprompt
#' @importFrom gtools rdirichlet permutations
#' @importFrom Rcpp sourceCpp
#' @importFrom processx run
#' @importFrom BiocGenerics cbind rownames colnames ncol nrow rowSums
#' @importFrom BiocGenerics colMeans rowMeans normalize subset colSums
#' @importFrom BiocGenerics as.data.frame
#' @useDynLib universalmotif
NULL
