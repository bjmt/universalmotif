## Benjamin Tremblay

# functions to write motifs to a file connection from various formats


write_meme <- function(motif_list, file_out, version = "4",
                       bkg = c(0.25, 0.25, 0.25, 0.25), strands = "+ -",
                       alphabet = NULL) {
  UseMethod("write_meme")
}
setGeneric("write_meme")  # for S4 classes
