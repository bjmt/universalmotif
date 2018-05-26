#' Convert \linkS4class{universalmotif} type.
#'
#' @param motifs Motif object or list.
#' @param type Character. Either PCM, PPM, PWM, ICM.
#' @param pseudoweight Numeric.
#' @param bkg Numeric.
#' @param IC_floor Logical.
#' @param IC_ceiling Logical.
#'
#' @return Motif object.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
convert_type <- function(motifs, type, pseudoweight, bkg, IC_floor = TRUE,
                         IC_ceiling = TRUE) {

  motif <- motifs
  if (class(motif) == "list") {
    margs <- list(type = type, IC_floor = IC_floor, IC_ceiling = IC_ceiling) 
    if (!missing(pseudoweight)) margs <- c(margs, list(pseudoweight = pseudoweight))
    if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))
    motif <- lapply(motif, function(x) do.call(convert_type,
                                               c(list(motifs = x), margs)))
    return(motif)

  }

  if (!type %in% c("PCM", "PPM", "PWM", "ICM")) {
    stop("unrecognized 'type'")
  }

  CLASS_IN <- .internal_convert(motif)
  motif <- convert_motifs(motif)

  in_type <- motif["type"]

  if (in_type == type) return(motif)

  if (missing(pseudoweight)) pseudoweight <- motif["pseudoweight"]
  if (missing(bkg)) bkg <- motif["bkg"]

  # PCM in:
  if (in_type == "PCM") {
    if (type == "PPM") {
      motif@motif <- apply(motif["motif"], 2, pcm_to_ppm,
                           pseudoweight = pseudoweight)
      motif["type"] <- "PPM"
    } else if (type == "PWM") {
      motif@motif <- apply(motif["motif"], 2, pcm_to_ppm,
                           pseudoweight = pseudoweight)
      motif@motif <- apply(motif["motif"], 2, ppm_to_pwm,
                           background = bkg,
                           smooth = any(motif["motif"] == 0),
                           pseudoweight = pseudoweight,
                           nsites = motif["nsites"])
      motif["type"] <- "PWM"
    } else if (type == "ICM") {
      motif@motif <- apply(motif["motif"], 2, pcm_to_ppm,
                           pseudoweight = pseudoweight)
      motif@motif <- apply(motif["motif"], 2, ppm_to_icm,
                           bkg = bkg, IC_floor = IC_floor,
                           IC_ceiling = IC_ceiling)
      motif["type"] <- "ICM"
    }
  }

  # PPM in:
  if (in_type == "PPM") {
    if (type == "PCM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pcm,
                           nsites = motif["nsites"])
      motif["type"] <- "PCM"
    } else if (type == "PWM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pwm,
                           background = bkg,
                           smooth = any(motif["motif"] == 0),
                           pseudoweight = pseudoweight,
                           nsites = motif["nsites"])
      motif["type"] <- "PWM"
    } else if (type == "ICM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_icm,
                           bkg = bkg, IC_floor = IC_floor,
                           IC_ceiling = IC_ceiling)
      motif["type"] <- "ICM"
    }
  }

  # PWM in:
  if (in_type == "PWM") {
    if (type == "PCM") {
      motif@motif <- apply(motif["motif"], 2, pwm_to_ppm,
                           background = bkg)
      motif@motif <- apply(motif["motif"], 2, ppm_to_pcm,
                           nsites = motif["nsites"])
      motif["type"] <- "PCM"
    } else if (type == "PPM") {
      motif@motif <- apply(motif["motif"], 2, pwm_to_ppm,
                           background = bkg)
      motif["type"] <- "PPM"
    } else if (type == "ICM") {
      motif@motif <- apply(motif["motif"], 2, pwm_to_ppm,
                           background = bkg)
      motif@motif <- apply(motif["motif"], 2, ppm_to_icm,
                           bkg = bkg, IC_floor = IC_floor,
                           IC_ceiling)
      motif["type"] <- "ICM"
    }
  }

  # ICM in:
  if (in_type == "ICM") {
    stop("motifs of type 'ICM' cannot be converted")
  }

  motif <- .internal_convert(motif, CLASS_IN)
  motif

}
