#' Example motif in `universalmotif` format.
#'
#' A simple DNA motif. To recreate this motif:
#' `create_motif("TATAWAW", nsites = numeric())`
#'
#' @format [universalmotif-class]
"examplemotif"

#' Another example motif in `universalmotif` format.
#'
#' A simple DNA motif with a non-empty `multifreq` slot. See
#' the 'Advanced usage' vignette, section 2, for manually recreating
#' this motif.
#'
#' @format [universalmotif-class]
"examplemotif2"

#' JASPAR2018 CORE database scores.
#'
#' For use with [compare_motifs()]. The precomputed scores allow for fast P-value estimation.
#' These scores were
#' generated using [make_DBscores()] with the JASPAR2018 CORE motif set,
#' with `normalise.scores = FALSE`. This particular set of scores is organized
#' as a list, with each list item being a `data.frame` of scores for a specific comparison
#' method. In each `data.frame` is the mean and sd of scores resulting between
#' the comparisons of JASPAR2018 CORE motifs with randomized motifs of the specified
#' `subject` and `target` motif length.
#'
#' @format `list` with three `data.frames`
"JASPAR2018_CORE_DBSCORES"

#' JASPAR2018 CORE database scores (normalised).
#'
#' For use with [compare_motifs()]. The precomputed scores allow for fast P-value estimation.
#' These scores were generated
#' using [make_DBscores()] with the JASPAR2018 CORE motif set, with
#' `normalise.scores = TRUE`. This particular set of scores is organized
#' as a list, with each list item being a `data.frame` of scores for a specific comparison
#' method. In each `data.frame` is the mean and sd of scores resulting between
#' the comparisons of JASPAR2018 CORE motifs with randomized motifs of the specified
#' `subject` and `target` motif length.
#'
#' @format `list` with three `data.frames`
"JASPAR2018_CORE_DBSCORES_NORM"

#' Arabidopsis promoters as a `DNAStringSet`.
#'
#' 50 Arabidopsis promoters, each 1000 bases long. See the 'Advanced usage'
#' vignette, section 7, for an example workflow describing extracting promoter
#' sequences.
#'
#' @format \code{\link{DNAStringSet}}
"ArabidopsisPromoters"

#' Arabidopsis motif in `universalmotif` format.
#'
#' Arabidopsis motif trained from [`ArabidopsisPromoters`] using
#' MEME version 4. This motif was generated at the command line using
#' the following command: `meme promoters.fa -revcomp -nmotifs 3  -mod anr -dna`.
#'
#' @format [universalmotif-class]
"ArabidopsisMotif"
