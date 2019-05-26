#' Example motif in `universalmotif` format.
#'
#' A simple DNA motif. To recreate this motif:
#' `create_motif("TATAWAW", nsites = numeric())`
#'
#' @format [universalmotif-class]
"examplemotif"

#' Another example motif in `universalmotif` format.
#'
#' A simple DNA motif with a non-empty `multifreq` slot.
#' To recreate to this motif:
#' `add_multifreq(examplemotif, DNAStringSet(rep(c("CAAAACC", "CTTTTCC"), 3)))`
#'
#' @format [universalmotif-class]
"examplemotif2"

#' JASPAR2018 CORE database scores
#'
#' For use with [compare_motifs()]. The precomputed scores allow for fast P-value estimation.
#' These scores were
#' generated using [make_DBscores()] with the JASPAR2018 CORE motif set,
#' with `normalise.scores = FALSE`. This particular set of scores is organized
#' as a list, with each list item being a `data.frame` of scores for a specific comparison
#' method. In each `data.frame` is the location and scale of scores resulting from a
#' logistic distribution using the the comparisons of JASPAR2018 CORE motifs with
#' randomized motifs of the specified
#' `subject` and `target` motif length. Created using [make_DBscores()] from
#' universalmotif v1.4.0. The parameters used can be seen via
#' `JASPAR2018_CORE_DBSCORES$args`.
#'
#' @format `list` with ten `data.frame`s and an additional list of parameters.
"JASPAR2018_CORE_DBSCORES"

#' JASPAR2018 CORE database scores (normalised)
#'
#' For use with [compare_motifs()]. The precomputed scores allow for fast P-value estimation.
#' These scores were generated
#' using [make_DBscores()] with the JASPAR2018 CORE motif set, with
#' `normalise.scores = TRUE`. This particular set of scores is organized
#' as a list, with each list item being a `data.frame` of scores for a specific comparison
#' method. In each `data.frame` is the location and scale of scores resulting from a
#' logistic distribution using the the comparisons of JASPAR2018 CORE motifs with
#' randomized motifs of the specified
#' `subject` and `target` motif length. Created using [make_DBscores()] from
#' universalmotif v1.4.0. The parameters used can be seen via
#' `JASPAR2018_CORE_DBSCORES_NORM$args`.
#'
#' @format `list` with ten `data.frame`s and an additional list of parameters.
"JASPAR2018_CORE_DBSCORES_NORM"

#' Arabidopsis promoters as a `DNAStringSet`.
#'
#' 50 Arabidopsis promoters, each 1000 bases long. See the
#' "Sequence manipulation and scanning"
#' vignette, section 9, for an example workflow describing extracting promoter
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
