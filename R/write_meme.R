#######################################################################
#                          Benjamin Tremblay                          #
#                             2017-02-22                              #
#           write motifs to a file connection in MEME format            #
#######################################################################

######################
#  public functions  #
######################

write_meme <- function(motif_list, file_out, version = "4",
                       bkg = c(0.25, 0.25, 0.25, 0.25), strands = "+ -") {
  preamble <- meme_preamble(motif_list, version, bkg, strands)

  # outfile <- file(file_out)
  # writeLines(converted_motif, con = file_out)
  # # in this case converted_motif is a vector of chars, each representing a line
  # close(outfile)
}
#-----------------------------------------------------------

########################
#  internal functions  #
########################

#-----------------------------------------------------------
meme_preamble <- function(motif_list, version, bkg, strands) {
  alphabet <- get_alph(motif_list)
  alphabet2 <- paste(alphabet, collapse = "")
  line_1 <- paste0("MEME version ", version)
  line_2 <- paste0("ALPHABET= ", alphabet2)
  if (!is.null(strands)) line_3 <- paste0("strands: ", strands) else {
    line_3 <- NA
  }
  if (!is.null(bkg)) {
    line_4 <- "Background letter frequencies"
    line_5 <- paste(alphabet, bkg)
    line_5 <- paste(line_5, collapse = " ")
  } else {
    line_4 <- NA
    line_5 <- NA
  }
  preamble <- c(line_1, line_2, line_3, line_4, line_5)
  preamble <- preamble[!is.na(preamble)]
}
#-----------------------------------------------------------
get_alph <- function(motif_list) {
  if (class(motif_list[[1]]) == "matrix") {
     
    return(colnames(motif_list[[1]]))
  }
}
