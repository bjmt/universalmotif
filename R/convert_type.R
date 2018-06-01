#' Convert \linkS4class{universalmotif} type.
#'
#' @param motifs Motif object or list.
#' @param type Character. Either PCM, PPM, PWM, ICM.
#' @param pseudoweight Numeric.
#' @param bkg Numeric.
#' @param BPPARAM Param for bplapply.
#'
#' @return Motif object.
#'
#' @examples
#' jaspar.ppm <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                       package = "universalmotif"))
#' jaspar.pwm <- convert_type(jaspar.ppm, type = "PWM")
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
convert_type <- function(motifs, type, pseudoweight, bkg,
                         BPPARAM = bpparam()) {

  IC_floor <-  TRUE
  IC_ceiling <- FALSE
  motif <- motifs
  if (class(motif) == "list") {
    margs <- list(type = type) 
    if (!missing(pseudoweight)) margs <- c(margs, list(pseudoweight = pseudoweight))
    if (!missing(bkg)) margs <- c(margs, list(bkg = bkg))
    motif <- bplapply(motif, function(x) do.call(convert_type,
                                                 c(list(motifs = x), 
                                                   margs)),
                      BPPARAM = BPPARAM)
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
                           IC_ceiling = IC_ceiling,
                           nsites = motif["nsites"],
                           schneider_correction = schneider_correction)
      motif["type"] <- "ICM"
    }
  }

  # PPM in:
  if (in_type == "PPM") {
    if (type == "PCM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pcm,
                           nsites = motif["nsites"])
      motif["type"] <- "PCM"
      if (length(motif["nsites"]) == 0) motif["nsites"] <- 100
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
                           IC_ceiling = IC_ceiling,
                           nsites = motif["nsites"],
                           schneider_correction = schneider_correction)
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
      if (length(motif["nsites"]) == 0) motif["nsites"] <- 100
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
                           IC_ceiling = IC_ceiling,
                           nsites = motif["nsites"],
                           schneider_correction = schneider_correction)
      motif["type"] <- "ICM"
    }
  }

  # ICM in:
  if (in_type == "ICM") {
    motif@motif <- apply(motif["motif"], 2, icm_to_ppm)
    if (type == "PCM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pcm,
                           nsites = motif["nsites"])
      if (length(motif["nsites"]) == 0) motif["nsites"] <- 100
      motif["type"] <- "PCM"
    } else if (type == "PPM") {
      motif["type"] <- "PPM"
    } else if (type == "PWM") {
      motif@motif <- apply(motif["motif"], 2, ppm_to_pwm,
                           background = bkg,
                           smooth = any(motif["motif"] == 0),
                           pseudoweight = pseudoweight,
                           nsites = motif["nsites"])
      motif["type"] <- "PWM"
    }
  }

  motif <- .internal_convert(motif, CLASS_IN)
  motif

}
