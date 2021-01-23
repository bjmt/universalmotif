#' Plot logos from numeric matrices.
#'
#' This function provides the plotting capabilities of [view_motifs()] without
#' requiring `universalmotif`-class objects. Instead, it takes a numeric matrix
#' with row names as input.
#'
#' @param x A numeric matrix with row names. The row names can be a mix of
#'    different character lengths.
#' @param colour.scheme `character` A named character vector of colour names.
#'    Provide colours for individual letters, even if the row names are made
#'    up of multiple characters.
#' @param x.spacer `numeric(1)` Add horizontal spacing between letters. The
#'    number is taken as the fraction of the width of an individual position.
#'    Increasing this value is recommended for letters made up of multiple
#'    characters.
#'
#' @return A `ggplot` object. If you wish to plot the data yourself from
#'    polygon paths, access them using `$data` on the output object.
#'
#' @examples
#' ## Feel free to mix and match row name character lengths.
#' data(examplemotif)
#' toplot <- examplemotif["motif"]
#' rownames(toplot)[1] <- "AA"
#' view_logo(toplot)
#'
#' @seealso [view_motifs()]
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @inheritParams view_motifs
#' @export
view_logo <- function(x, fontDF = NULL, colour.scheme = NULL, min.height = 0.01,
  x.spacer = 0.04, y.spacer = 0.01, sort.positions = FALSE,
  sort.positions.decreasing = TRUE, fit.to.height = NULL) {

  if (!is.matrix(x))
    stop("`x` must be a matrix")

  if (is.null(rownames(x)))
    stop("`x` must have row names")

  if (is.null(fontDF)) fontDF <- fontDFroboto
  fontDF <- as.data.frame(fontDF)

  alph <- rownames(x)
  alph <- sort_unique_cpp(unlist(lapply(alph, function(x) strsplit(x, "", TRUE)[[1]])))

  if (!is.null(colour.scheme) && any(!names(colour.scheme) %in% alph))
    stop(wmsg("colour.scheme must be a named vector with all possible row names present"))

  if (is.null(colour.scheme)) {
    colour.scheme <- structure(rep("black", length(alph)), names = alph)
  }

  plotdata <- make_matrix_polygon_data(x, fontDF, min.height, x.spacer, y.spacer,
    sort.positions, sort.positions.decreasing, fit.to.height)

  limits_x <- c(1 - 0.501, ncol(x) + 0.501)

  ggplot(plotdata, aes(.data$x, .data$y, group = .data$letter.id,
      fill = .data$group)) +
    geom_polygon() +
    scale_x_continuous(limits = limits_x, expand = c(0.02, 0)) +
    scale_fill_manual(values = colour.scheme[alph]) +
    theme_void() +
    theme(legend.position = "none")

}
