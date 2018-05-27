#' Motif distance score matrix.
#'
#' MotIV implementation of the STAMP motif analysis software.
#'
#' @param motifs List of motifs
#' @param metric Distance metric. One of PCC (Pearson Correlation Coefficient),
#'               ALLR (Average Log-Likelihood Ratio), ALLR_LL (ALLR with low
#'               score limit of 2), CS (Chi-Squared), KL (Kullback-Lieber),
#'               SSD (Sum of Squared Distances).
#' @param align Alignment strategy. One of NW (Needleman-Wunsch), 
#'              SW (Smith-Waterman with linear gaps), SWA (Smith-Waterman
#'              with affine gaps), SWU (Smith-Waterman with ungapped and
#'              alignments are extended).
#' @param gap_open Gap open penalty.
#' @param gap_extend Gap extension penaly.
#'
#' @return List of motifs or motif object.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.distances <- motif_dist(jaspar)
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
motif_dist <- function(motifs, metric = "PCC", align = "SWU",
                       gap_open = 1, gap_extend = 0.5) {

  motifs <- convert_motifs(motifs)

  motif_types <- vapply(motifs, function(x) x["type"], character(1))
  if (any(motif_types != "PWM")) {
    motifs[motif_types != "PWM"] <- lapply(motifs[motif_types != "PWM"],
                                           function(x) convert_type(x, "PWM"))
    # motifs <- lapply(motifs, function(x) convert_type(x, "PWM"))
  }

  motif_names <- vapply(motifs, function(x) x["name"], character(1))
  motifs <- lapply(motifs, function(x) x["motif"])
  names(motifs) <- motif_names

  jaspar.scores <- readDBScores(system.file("extdata",
                                            "jaspar2010_PCC_SWU.scores",
                                            package = "MotIV"))
  distances <- motifDistances(motifs, cc = metric, align = align,
                              go = gap_open, ge = gap_extend,
                              DBscores = jaspar.scores)

  as.matrix(distances)

}
