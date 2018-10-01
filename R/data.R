#' Example motif.
#'
#' A simple DNA motif.
#'
#' @format [universalmotif-class]
"examplemotif"

#' Another example motif.
#'
#' A simple DNA motif with a non-empty `multifreq` slot.
#'
#' @format [universalmotif-class]
"examplemotif2"

#' JASPAR2018 CORE database scores.
#'
#' For use with [compare_motifs()]. Unormalised.
#'
#' @format `list` with three `data.frames`
"JASPAR2018_CORE_DBSCORES"

#' JASPAR2018 CORE databse scores (normalised).
#'
#' For use with [compare_motifs()]. Normalised.
#'
#' @format `list` with three `data.frames`
"JASPAR2018_CORE_DBSCORES_NORM"

#' Arabidopsis promoters.
#'
#' 50 Arabidopsis promoters, each 1000 bases long.
#'
#' @format [Biostrings::DNAStringSet]
"ArabidopsisPromoters"

#' Arabidopsis motif.
#'
#' Arabidopsis motif trained from [`ArabidopsisPromoters`] using
#' MEME version 4.
#'
#' @format [universalmotif-class]
"ArabidopsisMotif"
