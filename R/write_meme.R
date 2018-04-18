######################################################################
## Benjamin Tremblay
##
## write_meme
##
######################################################################

#' @include universalmotif.R 
#' @export
write_meme <- function(motifs, file_out, version = 4,
                       bkg = NULL, strands = NULL, alphabet = NULL) {
  motif_list <- lapply(motifs, convert_motifs)
  write_meme_(motif_list = motif_list, file_out = file_out, 
              version = version, bkg = bkg, strands = strands,
              alphabet = alphabet)
  invisible(NULL)
}

write_meme_ <- function(motif_list, file_out, alphabet, strands,
                       version, bkg) {

  # need to check that bkg freq, letters, and strands are the same
  if (is.null(alphabet)) {
    alphabets <- vapply(motif_list, function(x) motif_slots(x, "alphabet"),
                        character(1))
    alphabet <- unique(alphabets)
    if (length(alphabet) != 1) stop("All motifs must share a common alphabet.",
                                    " Set 'alphabet' to override")
    alph_letters <- rownames(motif_slots(motif_list[[1]], "motif"))
  }
  if (is.null(bkg)) {
    bkgs <- sapply(motif_list, function(x) motif_slots(x, "bkg"))
    bkg <- apply(bkgs, 1, unique)
    if (class(bkg) == "list") stop("All motifs must share common background ", 
                                   "letter frequencies. Set 'bkg' to override")
  }
  if (is.null(strands)) {
    strands <- sapply(motif_list, function(x) motif_slots(x, "strand"))
    strands <- apply(strands, 1, unique)
  }

  pre1 <- paste("MEME version", version)
  pre3 <- paste("ALPHABET=", paste(alph_letters, collapse = ""))
  pre5 <- paste("strands:", paste(strands, collapse = " "))
  pre7 <- "Background letter frequencies"
  pre8 <- paste(alph_letters, bkg, collapse = " ")
  
  preamble <- c(pre1, "", pre3, "", pre5, "", pre7, pre8)

  mots <- ""
  for (i in seq_along(motif_list)) {
    motif <- motif_list[[i]]
    mot1 <- paste("MOTIF", motif_slots(motif, "name"))
    mot2.1 <- paste("alength=", length(alph_letters))
    mot2.2 <- paste("w=", ncol(motif_slots(motif, "motif")))
    mot2.3 <- paste("n=", motif_slots(motif, "nsites"))
    evalue <- motif_slots(motif, "eval")
    if (length(evalue) == 0) evalue <- 0
    mot2.4 <- paste("E=", evalue) 
    mot2 <- paste("letter-probability matrix:", mot2.1, mot2.2, mot2.3, mot2.4)
    type <- motif_slots(motif, "type")
    if (type == "PPM") {
      mot3 <- t(motif_slots(motif, "motif"))
    } else {
      mot3 <- convert_type(motif, "PPM") 
      mot3 <- t(motif_slots(mot3, "motif"))
    }
    mot3 <- apply(mot3, 1, function(x) paste("", sprintf("%#.6f", x),
                                             collapse = " "))
    mot4 <- c(mot1, mot2, mot3, "")
    mots <- c(mots, mot4)
  }

  writeLines(c(preamble, mots[-length(mots)]), file_con <- file(file_out))
  close(file_con)

}
