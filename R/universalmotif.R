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
#' @importFrom ggplot2 aes element_blank facet_wrap ggplot theme ylab xlab
#' @importFrom ggplot2 geom_line geom_point geom_hline geom_text theme_bw
#' @importFrom ggseqlogo ggseqlogo geom_logo theme_logo
#' @importFrom ggtree ggtree geom_tiplab geom_tiplab2 groupOTU %<+%
#' @importFrom ggtree geom_tippoint
#' @importFrom stats as.dist hclust runif rnorm chisq.test fisher.test sd fft
#' @importFrom stats p.adjust shapiro.test t.test wilcox.test pnorm quantile
#' @importFrom stats rpois dnorm ecdf var qnorm plogis pweibull ks.test
#' @importFrom MASS fitdistr
#' @importFrom utils read.table menu packageDescription combn capture.output
#' @importFrom Biostrings width consensusMatrix BString AAString matchPWM
#' @importFrom Biostrings PWM DNAStringSet BStringSet AAStringSet
#' @importFrom Biostrings DNAString RNAString RNAStringSet AA_STANDARD
#' @importFrom Biostrings DNA_BASES RNA_BASES DNA_ALPHABET as.matrix
#' @importFrom Biostrings oligonucleotideTransitions trinucleotideFrequency
#' @importFrom Biostrings dinucleotideFrequency oligonucleotideFrequency
#' @importFrom Biostrings reverseComplement writeXStringSet seqtype AA_ALPHABET
#' @importFrom ape as.phylo
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @importFrom processx run
#' @importFrom BiocGenerics cbind rownames colnames ncol nrow rowSums
#' @importFrom BiocGenerics colMeans rowMeans normalize subset colSums
#' @importFrom BiocGenerics as.data.frame
#' @importFrom S4Vectors safeExplode wmsg DataFrame Rle
#' @importFrom yaml as.yaml yaml.load
#' @useDynLib universalmotif
NULL

# Notes:
#   - Printing abridged strings: S4Vectors:::sketchStr()

# Idea: create a masking motif, then find spots to mask with scan_sequences().
# Any different mask.symbol means converting to BString.
# Work with Masked*String objects?
# Replace letters: chartr, replaceAt, replaceLetterAt, subseq
# Find letters: matchPattern, matchPWM
# injectHardMask: converts MaskedXString to XString, using a character such
#                 as "+" for masked letters
