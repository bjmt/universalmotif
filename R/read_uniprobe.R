######################################################################
## Benjamin Tremblay
##
## Read UNIPROBE motifs from a text file.
##
## (The code for this function is quite ugly and inefficient; however
## this is due to how different the old and new PWM formats are.)
######################################################################

#' [UNDER CONSTRUCTION] Load UNIPROBE motifs from a text file.
#'
#' Support for uniprobe-style motifs. Each motif has a header and a number
#' matrix. Only DNA. Only PWM files.
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
#' @author Benjamin Tremblay, \email{b2trembl@@uwaterloo.ca}
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
  if (verbose) cat("Found", length(uni_names_indices), "motifs.\n")
  # the following function is to support the new version of uniprobe PWM;
  # however I am somewhat unsure as to how to code this. this will do for now.
  uni_names_indices <- lapply(uni_names_indices, uniprobe_names2,
                              uniprobe_raw)
  uni_names_indices <- mapply(function(x, y) {
                                if (!is.null(x[1])) return(x)
                                x[1] <- y
                                return(x)
                              },
                              uni_names_indices, seq_along(uni_names_indices),
                              SIMPLIFY = FALSE)

  motifs <- lapply(uni_names_indices, uniprobe_load, uniprobe_raw)

  mot_names <- vapply(uni_names_indices, function(x) x[1], character(1))

  names(motifs) <- mot_names

  motifs <- mapply(uniprobe_to_umot, motifs, mot_names, SIMPLIFY = FALSE)

  return(motifs)

}

######################################################################
######################################################################

uniprobe_names <- function(uni_all, uniprobe_raw) {

  tocheck <- as.integer(uni_all[1]) - 1
  for (i in tocheck:1) {
    if (grepl("Enrichment score matrix", uniprobe_raw[i],
              fixed = TRUE)) return(NULL)
    if (grepl("Energy matrix", uniprobe_raw[i], fixed = TRUE)) return(NULL)
    #if (grepl("T:", uniprobe_raw[i], fixed = TRUE)) return(NULL)
    if (uniprobe_raw[i] != "" && !grepl("^P0\\s+", uniprobe_raw[i]) &&
        !grepl("^PO\\s+", uniprobe_raw[i]) &&
        !grepl("A:", uniprobe_raw[i], fixed = TRUE) &&
        !grepl("C:", uniprobe_raw[i], fixed = TRUE) &&
        !grepl("G:", uniprobe_raw[i], fixed = TRUE) &&
        !grepl("T:", uniprobe_raw[i], fixed = TRUE)) {
      return(c("thename" = uniprobe_raw[i], "nameindex" = as.character(i), 
               "A" = uni_all[1], "C" = uni_all[2], "G" = uni_all[3],
               "T" = uni_all[4]))
    }
  }
  stop("could not find any motif names", call. = FALSE)

}

uniprobe_names2 <- function(uni_names_indices, uniprobe_raw) {

  if (uni_names_indices[1] != "Probability matrix") return(uni_names_indices)
  tocheck <- as.integer(uni_names_indices[2]) - 1
  for (i in tocheck:1) {
    if (grepl("Probability matrix", uniprobe_raw[i], fixed = TRUE)) {
      uni_names_indices[1] <- NULL
      return(uni_names_indices)
    } 
    # unsure whether I should cut all of this !grepl code and just have
    # grepl("#").. but then I'm not sure how to handle cases with no header
    if (!grepl("Enrichment score matrix", uniprobe_raw[i], fixed = TRUE) &&
        !grepl("Energy matrix", uniprobe_raw[i], fixed = TRUE) &&
        !grepl("^P0\\s+", uniprobe_raw[i]) &&
        !grepl("^PO\\s+", uniprobe_raw[i]) &&
        !grepl("A:", uniprobe_raw[i], fixed = TRUE) &&
        !grepl("C:", uniprobe_raw[i], fixed = TRUE) &&
        !grepl("G:", uniprobe_raw[i], fixed = TRUE) &&
        !grepl("T:", uniprobe_raw[i], fixed = TRUE) &&
        uniprobe_raw[i] != "") {
      ifelse(grepl("#", uniprobe_raw[i], fixed = TRUE), msplit <- 3, 
             msplit <- 1)
      uni_names_indices[1] <- strsplit(uniprobe_raw[i],
                                       split = "\\s+")[[1]][msplit]
      return(uni_names_indices)
    }
  }

}

uniprobe_load <- function(uni_indicies, uniprobe_raw) {

  motif <- uniprobe_raw[as.integer(uni_indicies[3:6])] 
  motif <- as.matrix(read.table(text = motif, row.names = 1))
  rownames(motif) <- c("A", "C", "G", "T")
  colnames(motif) <- NULL
  return(motif)

}

uniprobe_to_umot <- function(motif, name) {

  motif <- new("universalmotif", name = name, motif = motif, type = "PPM")

}
