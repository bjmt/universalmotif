######################################################################
## Benjamin Tremblay
##
## write_homer
##
######################################################################

#' @include universalmotif.R 
#' @export
write_homer <- function(motifs, file_out, threshold = NULL) {
  motif_list <- lapply(motifs, convert_motifs)
  write_homer_(motif_list = motif_list, file_out = file_out,
               threshold = threshold)
  invisible(NULL)
}

write_homer_ <- function(motif_list, file_out, threshold) {

  mots <- vector()
  for (i in seq_along(motif_list)) {
    motif <- motif_list[[i]]
    mot1 <- paste0(">", motif_slots(motif, "consensus"))
    mot2 <- motif_slots(motif, "name")
    if (is.null(threshold)) {
      if ("detection_threshold" == names(motif_slots(motif, "extranum"))) {
        threshold <- motif_slots(motif, "extranum")
        mot3 <- threshold
      } else {
        threshold  <- 5
        mot3 <- threshold 
      }
    } else if (length(threshold) == 1) {
      mot3 <- threshold
    } else mot3 <- threshold[i]
    mot4 <- convert_type(motif, "PPM")
    mot4 <- t(motif_slots(mot4, "motif"))
    mot4 <- apply(mot4, 1, function(x) paste(sprintf("%#.3f", x),
                                             collapse = "\t"))
    mot123 <- paste(mot1, mot2, mot3, sep = "\t")
    mots <- c(mots, mot123, mot4)
  }

  writeLines(mots, file_con <- file(file_out))
  close(file_con)

}
