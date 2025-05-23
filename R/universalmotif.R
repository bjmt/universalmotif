#' universalmotif: Import, Modify and Export Motifs with R
#'
#' @description
#' A collection of utility functions for working with motifs.
#'
#' @docType package
#' @name universalmotif-pkg
#' @aliases universalmotif-pkg
#'
#' @importFrom methods getClass is new slot "slot<-" validObject as
#' @importMethodsFrom methods initialize show
#' @importFrom ggplot2 aes element_blank facet_wrap ggplot theme ylab xlab
#' @importFrom ggplot2 geom_line geom_point geom_hline geom_text theme_bw ylim
#' @importFrom ggplot2 margin scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 element_line element_text waiver geom_polygon
#' @importFrom ggplot2 theme_minimal scale_fill_manual theme_void
#' @importFrom grid unit
#' @importFrom rlang .data
#' @importFrom stats as.dist hclust runif rnorm chisq.test fisher.test sd fft
#' @importFrom stats p.adjust shapiro.test t.test wilcox.test pnorm quantile
#' @importFrom stats rpois dnorm ecdf var qnorm plogis pweibull ks.test pbinom
#' @importFrom stats cutree xtabs
#' @importFrom MASS fitdistr
#' @importFrom utils read.table menu packageDescription combn capture.output
#' @importFrom Biostrings width consensusMatrix BString AAString matchPWM
#' @importFrom Biostrings PWM DNAStringSet BStringSet AAStringSet
#' @importFrom Biostrings DNAString RNAString RNAStringSet AA_STANDARD
#' @importFrom Biostrings DNA_BASES RNA_BASES DNA_ALPHABET as.matrix
#' @importFrom Biostrings oligonucleotideTransitions trinucleotideFrequency
#' @importFrom Biostrings dinucleotideFrequency oligonucleotideFrequency
#' @importFrom Biostrings reverseComplement writeXStringSet seqtype AA_ALPHABET
#' @importFrom Biostrings mask injectHardMask
#' @importFrom IRanges stack IRanges findOverlaps
#' @importFrom Rcpp sourceCpp
#' @importFrom BiocGenerics cbind rownames colnames ncol nrow
#' @importFrom BiocGenerics normalize subset as.data.frame
#' @importFrom MatrixGenerics colMeans rowMeans colSums rowSums
#' @importFrom S4Vectors safeExplode wmsg DataFrame Rle aggregate
#' @importFrom S4Vectors isTRUEorFALSE
#' @importFrom yaml as.yaml yaml.load
#' @useDynLib universalmotif
"_PACKAGE"

# Notes:
#   - Printing abridged strings: S4Vectors:::sketchStr()
