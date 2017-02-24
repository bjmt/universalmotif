#######################################################################
#                          Benjamin Tremblay                          #
#                             2017-02-22                              #
#           write motifs to a file connection in MEME format          #
#######################################################################

######################
#  public functions  #
######################

write_meme <- function(motif_list, file_out, version = "4",
                       bkg = c(0.25, 0.25, 0.25, 0.25), strands = "+ -",
                       alphabet = NULL) {
  if (!is.list(motif_list)) motif_list <- list(motif_list)
  motif_class <- class_check(motif_list)
  preamble <- meme_preamble(motif_list, version, bkg, strands, motif_class)
  motifs <- write_motifs(motif_list, motif_class)
  # TODO: need to append empty line between motifs + motif titles
  motifs <- unlist(motifs)
  final <- list(preamble, motifs)
  final <- unlist(final)
  file_con <- file(file_out)
  writeLines(final, con = file_con)
  # # in this case converted_motif is a vector of chars, each representing a line
  close(file_con)
}  # don't forget to add empty lines between preamble + motifs, and between motifs
#-----------------------------------------------------------

########################
#  internal functions  #
########################

#-----------------------------------------------------------
meme_preamble <- function(motif_list, version, bkg, strands, motif_class) {
  if (!is.list(motif_list)) motif_list <- list(motif_list)
  alphabet <- get_alph(motif_list, motif_class)
  alphabet2 <- paste(alphabet, collapse = "")
  line_1 <- paste("MEME version", version)
  line_2 <- paste("ALPHABET=", alphabet2)
  if (!is.null(strands)) line_3 <- paste("strands:", strands) else {
    line_3 <- NA
  }
  if (!is.null(bkg)) {
    line_4 <- "Background letter frequencies"
    if (length(alphabet) != length(bkg)) {
      stop("Alphabet and background frequencies do not match.")
    }
    line_5 <- paste(alphabet, bkg)
    line_5 <- paste(line_5, collapse = " ")
  } else {
    line_4 <- NA
    line_5 <- NA
  }
  line_6 <- ""
  preamble <- c(line_1, line_2, line_3, line_4, line_5, line_6)
  preamble <- preamble[!is.na(preamble)]
}
#-----------------------------------------------------------
get_alph <- function(motif_list, motif_class) {
  if (motif_class == "matrix") {
    if (all(unlist(sapply(motif_list, rowSums)) < 1.01) &&
            all(unlist(sapply(motif_list, rowSums)) > 0.99)) {
      return(colnames(motif_list[[1]]))
    } else if (all(unlist(sapply(motif_list, colSums)) < 1.01) &&
                   all(unlist(sapply(motif_list, colSums)) > 0.99)) {
      return(rownames(motif_list[[1]]))
    }
    stop("Neither columns or rows add up to 1.")
  }
  if (motif_class == "pwm") {
    if (attr(attr(motif_list[[1]], "class"), "package") == "seqLogo") {
      return(rownames(attr(motif_list[[1]], "pwm")))
    }
  }
}
#-----------------------------------------------------------
write_motifs <- function(motif_list, motif_class) {
  final_list <- lapply(motif_list, write_motif, motif_class = motif_class)
}
#-----------------------------------------------------------
write_motif <- function(motif, motif_class) {
  if (motif_class == "matrix") motif <- type_matrix(motif)
  if (motif_class == "pwm") {
    if (attr(attr(motif, "class"), "package") == "seqLogo") {
      motif <- type_seqLogo_pwm(motif)
    } else {
      motif <- type_matrix(as.matrix(motif))
      warning("Unrecognized motif class detected; tried to coerce to class matrix.")
    }
  }
  return(motif)
}
#-----------------------------------------------------------
class_check <- function(motif_list) {
  class_all <- sapply(motif_list, class)
  if (length(unique(class_all)) != 1) {
    stop("Motifs are not all of the same class.")
  }
  return(unique(class_all))
}
#-----------------------------------------------------------
type_seqLogo_pwm <- function(motif) {

}
#-----------------------------------------------------------
type_matrix <- function(motif) {
  if (all(rowSums(motif) < 1.01) &&
      all(rowSums(motif) > 0.99)) type <- 1
  if (all(colSums(motif) < 1.01) &&
      all(colSums(motif) > 0.99)) type <- 2
  motif <- apply(motif, type, function(x) paste(as.vector(x), collapse = " "))
  return(motif)
}
