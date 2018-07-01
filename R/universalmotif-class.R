#' universalmotif: Motif class.
#'
#' Container for motif objects. See \code{\link{create_motif}} for creating
#' motifs.
#'
#' @slot name Character. Length 1.
#' @slot altname Character. Length 0 or 1.
#' @slot family Character. Length 0 or 1.
#' @slot organism Character. Length 0 or 1.
#' @slot motif Matrix.
#' @slot alphabet Character. Length 1.
#' @slot type Character. Length 1.
#' @slot icscore Numeric. Length 1.
#' @slot nsites Numeric. Length 0 or 1.
#' @slot pseudocount Numeric. Length 1.
#' @slot bkg Numeric. Length equal to number of letters in alphabet.
#' @slot bkgsites Numeric. Lenght 0 or 1.
#' @slot consensus Character. Length 0 or 1.
#' @slot strand Character. Length 1.
#' @slot pval Numeric. Length 0 or 1.
#' @slot qval Numeric. Length 0 or 1.
#' @slot eval Numeric. Length 0 or 1.
#' @slot hmmfirst Matrix.
#' @slot hmmsecond Matrix.
#' @slot extrainfo Character. Length 0 or more.
#'
#' @return A motif object of class \linkS4class{universalmotif}.
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name universalmotif-class
#' @rdname universalmotif-class
#' @exportClass universalmotif
universalmotif <- setClass("universalmotif",
                           slots = list(name = "character",
                                        altname = "character",
                                        family = "character",
                                        organism = "character",
                                        motif = "matrix",
                                        alphabet = "character",
                                        type = "character", icscore = "numeric",
                                        nsites = "numeric",
                                        pseudocount = "numeric",
                                        bkg = "numeric",
                                        bkgsites = "numeric",
                                        consensus = "character",
                                        strand = "character", pval = "numeric",
                                        qval = "numeric", eval = "numeric",
                                        hmmfirst = "matrix",
                                        hmmsecond = "matrix",
                                        extrainfo = "character"))

setValidity("universalmotif",
            function(object) {
            
            msg <- vector()
            valid <- TRUE

            ## character slot length checks

            # mandatory
            char_checks1 <- c("name", "type", "alphabet", "strand")
            # optional
            char_checks2 <- c("altname", "family", "consensus", "organism")

            char_results1 <- vapply(char_checks1,
                                   function(x) length(object[x]) == 1,
                                   logical(1))
            char_results2 <- vapply(char_checks1,
                                   function(x) length(object[x]) <= 1,
                                   logical(1))

            char_checks_all <- c(char_checks1, char_checks2)
            char_results_all <- c(char_results1, char_results2)

            if (any(isFALSE(char_results_all))) {
              valid <- FALSE
              for (check in char_checks_all[isFALSE(char_results_all)]) {
                msg <- c(msg, paste0("motif '", check,
                                     "' must be a character vector of length 1"))
              }
            }
            
            ## character slot specific checks

            # alphabets <- c("DNA", "RNA", "AA", "custom")
            types <- c("PCM", "PPM", "PWM", "ICM")
            strands <- c("+", "-", "+-", "-+")

            # if (!object["alphabet"] %in% alphabets) {
              # valid <- FALSE
              # msg <- c(msg,
                       # "motif 'alphabet' must be either 'DNA', 'RNA', 'AA', or 'custom'")
            # }

            if (!object["type"] %in% types) {
              valid <- FALSE
              msg <- c(msg, "motif 'type' must be either 'PCM', 'PPM', 'PWM', or 'ICM'")
            }

            if (!object["strand"] %in% strands) {
              valid <- FALSE
              msg <- c(msg, "motif 'strand' must be either '+', '-', or '+-'")
            }

            ## numeric slot length checks

            num_checks <- c("pval", "qval", "eval", "bkgsites", "pseudocount",
                            "nsites")
            num_results <- vapply(num_checks,
                                  function(x) length(object[x]) <= 1,
                                  logical(1))

            if (any(isFALSE(num_results))) {
              valid <- FALSE
              for (check in num_checks[isFALSE(num_checks)]) {
                msg <- c(msg, paste0("motif '", check,
                                     "' must be a numeric vector of length 1"))
              }
            }

            ## bkg slot check

            bkg_check1 <- object["bkg"]
            bkg_check2 <- nrow(object["motif"])
            if (length(bkg_check1) > 0) {
              if (length(bkg_check1) != bkg_check2) {
                valid <- FALSE
                msg <- c(msg, "motif 'bkg' must be a numeric vector of length equal to the number of letters in motif")
              }
            }

            ## motif slot check

            mat_type <- object["type"]
            mat <- object["motif"]
            mat_colsums <- colSums(mat)

            # PCM
            if (mat_type == "PCM") {
              mat_nsites <- object["nsites"]
              if (length(unique(mat_colsums)) > 1) {
                warning("not all positions have identical count totals")
                # if (sd(unique(mat_colsums)) > 0.1) {
                  # valid <- FALSE
                  # msg <- c(msg, "motif of type 'PCM' must have equal column sums")
                # }
              }
              if (length(mat_nsites) > 0) {
                if (length(unique(mat_colsums)) == 1) {
                  if (unique(mat_colsums) != mat_nsites) {
                    valid <- FALSE
                    msg <- c(msg, "motif of type 'PCM' must have column sums equal to 'nsites'")
                  }
                }
              }
            }

            # PPM
            if (mat_type == "PPM") {
              if (any(mat_colsums > 1.01) || any(mat_colsums < 0.99)) {
                valid <- FALSE
                msg <- c(msg, "motif of type 'PPM' must have column sums equal to 1")
              }
              if (any(mat < 0)) {
                valid <- FALSE
                msg <- c(msg, "motif of type 'PPM' can only have positive probabilities")
              }
            }

            # check it matches alphabet
            alph <- object["alphabet"]
            if (alph %in% c("DNA", "RNA")) {
              if (nrow(mat) != 4) {
                valid <- FALSE
                msg <- c(msg,
                         "motif with 'alphabet' of type 'DNA' or 'RNA' can only have 4 rows")
              }
            } else if (alph == "AA") {
              if (nrow(mat) != 20) {
                valid <- FALSE
                msg <- c(msg,
                         "motif with 'alphabet' of type 'AA' can only have 20 rows")
              }
            } else if (alph != "custom") {
              if (nrow(mat) != length(strsplit(alph, "")[[1]])) {
                valid <- FALSE
                msg <- c(msg, paste0("motif with alphabet '", "' has an incorrect number of rows"))
              }
            }

            if (valid) TRUE else msg
            
            })
