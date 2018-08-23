#' Example motif.
#'
#' A simple DNA motif.
#'
#' @format \linkS4class{universalmotif}
"examplemotif"

#' Another example motif.
#'
#' A simple DNA motif with a non-empty \code{multifreq} slot.
#'
#' @format \linkS4class{universalmotif}
"examplemotif2"

#' JASPAR2018 CORE database scores.
#'
#' For use with \code{\link{compare_motifs}}. Unormalised.
#'
#' @format \code{list} with three \code{data.frames}
"JASPAR2018_CORE_DBSCORES"

#' JASPAR2018 CORE databse scores [normalised].
#'
#' For use with \code{\link{compare_motifs}}. Normalised.
#'
#' @format \code{list} with three \code{data.frames}
"JASPAR2018_CORE_DBSCORES_NORM"

#' Arabidopsis promoters.
#'
#' 50 Arabidopsis promoters, each 1000 bases long.
#'
#' @format \linkS4class{DNAStringSet}
"ArabidopsisPromoters"

#' Arabidopsis motif.
#'
#' Arabidopsis motif trained from \code{\link{ArabidopsisPromoters}} using
#' MEME v4.
#'
#' @format \linkS4class{universalmotif}
"ArabidopsisMotif"
