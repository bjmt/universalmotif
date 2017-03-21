#' universalmotif: importing and export motifs with R
#'
#' A collection of utility functions for working with motifs.
#'
#' @docType package
#' @name universalmotif
#'
#' @import methods
#' @import utils
#' @importFrom TFBSTools PFMatrix toPWM toICM
#' @importClassesFrom TFBSTools PFMatrix PWMatrix ICMatrix
#' @importFrom seqLogo makePWM
#' @importClassesFrom seqLogo pwm
#' @importClassesFrom motifStack pcm pfm
#' @importFrom PWMEnrich PFMtoPWM
#' @importClassesFrom PWMEnrich PWM
#' @importFrom Biostrings PWM width
#' @importClassesFrom Biostrings DNAStringSet RNAStringSet DNAString RNAString
#'    XString XStringSet AAString AAStringSet
NULL

# @import seqLogo
# @import motifStack
# @import TFBSTools
# @import MotifDb
# @import PWMEnrich

# Note: the most important thing is importing packages in the description. 
# Importing them via namespace only makes the functions available without
# typing `::`, but the package is still installed (as directed by the
# description import).

# Note 2: biocCheck still wants me to NAMESPACE import them..
