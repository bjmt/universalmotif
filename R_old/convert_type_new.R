#' Convert \linkS4class{universalmotif} type.
#'
#' @param motif Motif object.
#' @param type Character. Either PCM, PPM, PWM, ICM.
#' @param pseudoweight Numeric.
#' @param bkg Numeric.
#' @param IC_floor Logical.
#'
#' @return Motif object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
convert_type <- function(motif, type, pseudoweight, bkg, IC_floor = FALSE) {

  in_type <- motif["type"]

  if (in_type == type) return(motif)

  if (missing(pseudoweight)) pseudoweight <- motif["pseudoweight"]
  if (missing(bkg)) bkg <- motif["bkg"]

  # PCM in:
  if (in_type == "PCM") {
    if (type == "PPM") {
      motif["motif"] <- apply(motif["motif"], 2, pcm_to_ppm,
                              pseudoweight = pseudoweight)
      motif["type"] <- "PPM"
      return(motif)
    } else if (type == "PWM") {
      motif["motif"] <- apply(motif["motif"], 2, pcm_to_ppm,
                              pseudoweight = pseudoweight)
      motif["motif"] <- apply(motif["motif"], 2, ppm_to_pwm,
                              background = bkg,
                              pseudoweight = pseudoweight,
                              nsites = motif["nsites"])
      motif["type"] <- "PWM"
      return(motif)
    } else if (type == "ICM") {
      motif["motif"] <- apply(motif["motif"], 2, pcm_to_ppm,
                              pseudoweight = pseudoweight)
      motif["motif"] <- apply(motif["motif"], 2, ppm_to_icm,
                              bkg = bkg, IC_floor = IC_floor)
      motif["type"] <- "ICM"
      return(motif)
    }
  }

  # PPM in:
  if (in_type == "PPM") {
    if (type == "PCM") {
      motif["motif"] <- apply(motif["motif"], 2, ppm_to_pcm,
                              nsites = motif["nsites"])
      motif["type"] <- "PCM"
      return(motif)
    } else if (type == "PWM") {
      motif["motif"] <- apply(motif["motif"], 2, ppm_to_pwm,
                              background = bkg,
                              pseudoweight = pseudoweight,
                              nsites = motif["nsites"])
      motif["type"] <- "PWM"
      return(motif)
    } else if (type == "ICM")
  }

}
