######################################################################
## Benjamin Tremblay
##
## Read UNIPROBE motifs from a text file.
##
######################################################################

#' [UNDER CONSTRUCTION] Load UNIPROBE motifs from a text file.
#'
#' Support for uniprobe-style motifs. Each motif has a header and a number
#' matrix. Only DNA.
#'
#' @param motif_file Character.
#' @param verbose Logical.
#' @param mot_length_cutoff Integer.
#' @param out_class Character.
#'
#' @return a list of motif objects of the desired class.
#'
#' @examples
#'    motifs <- system.file("extdata", "uniprobe.txt",
#'                          package = "universalmotif")
#'    rmotifs <- read_uniprobe(motifs)
#'
#' @author Benjamin Tremblay, \email{b2trembl@uwaterloo.ca}
#' @include utils.R
#' @export
read_uniprobe <- function(motif_file, verbose = FALSE,
                          mot_length_cutoff = NULL, out_class = "matrix-2") {

  # check args
  check_logi_args(as.list(environment())[2])
  check_filter_args(as.list(environment())[3])
  check_out_class(out_class)

  # read file
  uniprobe_raw <- readLines(con <- file(motif_file)); close(con)
  if (length(uniprobe_raw) == 0) stop("could not read file, or file is empty",
                                      call. = FALSE)
  names(uniprobe_raw) <- seq_along(uniprobe_raw)

  # initial motif finding:
      # idea: get all A:, C:, G:, and T:; additional sorting of energy and 
      # enrichment matricies will have be done
  uniA <- uniprobe_raw[vapply(uniprobe_raw, function(x) 
                              grepl("A:", x, fixed = TRUE), logical(1))]
  uniC <- uniprobe_raw[vapply(uniprobe_raw, function(x) 
                              grepl("C:", x, fixed = TRUE), logical(1))]
  uniG <- uniprobe_raw[vapply(uniprobe_raw, function(x) 
                              grepl("G:", x, fixed = TRUE), logical(1))]
  uniT <- uniprobe_raw[vapply(uniprobe_raw, function(x) 
                              grepl("T:", x, fixed = TRUE), logical(1))]

  test1 <- c(length(uniA), length(uniC), length(uniG), length(uniT))
  if (!diff(range(test1)) < 1) {
    stop("could not find all letters; possible incomplete motifs",
         call. = FALSE)
  }

  uniAi <- names(uniA)
  uniCi <- names(uniC)
  uniGi <- names(uniG)
  uniTi <- names(uniT)

  uni_all <- mapply(function(a, b, c, d) c(a, b, c, d),
                    uniAi, uniCi, uniGi, uniTi, SIMPLIFY = FALSE)

  uni_names_indices <- lapply(uni_all, uniprobe_names, uniprobe_raw)
  uni_names_indices <- uni_names_indices[vapply(uni_names_indices, function(x)
                                                !is.null(x), logical(1))]

  # need a function to handle names returned as "Probability matrix"

  motifs <- lapply(uni_names_indices, uniprobe_load, uniprobe_raw)

  names(motifs) <- vapply(uni_names_indices, function(x) x[1], character(1))

  return(motifs)

}

######################################################################
######################################################################

uniprobe_names <- function(uni_all, uniprobe_raw) {

  tocheck <- as.integer(uni_all[1]) - 1
  finalname <- NULL
  for (i in tocheck:1) {
    if (grepl("Enrichment score matrix", uniprobe_raw[i],
              fixed = TRUE)) return(NULL)
    if (grepl("Energy matrix", uniprobe_raw[i], fixed = TRUE)) return(NULL)
    if (grepl("T:", uniprobe_raw[i], fixed = TRUE)) return(NULL)
    if (uniprobe_raw[i] != "" && !grepl("^P0", uniprobe_raw[i]) &&
        !grepl("^PO", uniprobe_raw[i])) {
      finalname <- uniprobe_raw[i]
      return(c("thename" = finalname, "nameindex" = as.character(i), 
               "A" = uni_all[1], "C" = uni_all[2], "G" = uni_all[3],
               "T" = uni_all[4]))
    }
  }
  stop("could not find any motif names", call. = FALSE)

}

uniprobe_load <- function(uni_indicies, uniprobe_raw) {

  motif <- uniprobe_raw[as.integer(uni_indicies[3:6])] 
  motif <- as.matrix(read.table(text = motif, row.names = 1))
  rownames(motif) <- c("A", "C", "G", "T")
  return(motif)

}
