#' universalmotif: importing and export motifs with R
#'
#' A collection of utility functions for working with motifs.
#'
#' @docType package
#' @name universalmotif-pkg
#'
#' @import methods
#' @import utils
#' @importFrom TFBSTools PFMatrix toPWM toICM
#' @importClassesFrom TFBSTools PFMatrix PWMatrix ICMatrix
#' @importFrom seqLogo seqLogo
#' @importClassesFrom seqLogo pwm
#' @importClassesFrom motifStack pcm pfm
#' @importFrom PWMEnrich PFMtoPWM
#' @importClassesFrom PWMEnrich PWM
#' @importFrom Biostrings PWM width
#' @importClassesFrom Biostrings DNAStringSet RNAStringSet DNAString RNAString
#'    XString XStringSet AAString AAStringSet
#' @importFrom MotIV makePWM
#' @importClassesFrom rGADEM motif
#' @importClassesFrom MotifDb MotifList
NULL

# @import seqLogo makePWM --> both seqLogo and MotIV export this function name
# @import motifStack
# @import TFBSTools
# @import MotifDb
# @import PWMEnrich
# @importClassesFrom MotifDb MotifList

# currently R CMD CHECK is complaining about not being able to find
# `pwm2-class`..
# `pwm2-class` is NOT exported by MotIV!

# Note: the most important thing is importing packages in the description. 
# Importing them via namespace only makes the functions available without
# typing `::`, but the package is still installed (as directed by the
# description import).

# Note 2: biocCheck still wants me to NAMESPACE import them..
