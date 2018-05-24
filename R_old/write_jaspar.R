######################################################################
## Benjamin Tremblay
##
## write_jaspar
##
######################################################################

#' @include universalmotif.R
#' @export
write_jaspar <- function(motifs, file_out, matrix_id = NULL) {
  motif_list <- lapply(motifs, convert_motifs)
  write_jaspar_(motif_list = motif_list, file_out = file_out,
                matrix_id = matrix_id)
  invisible(NULL)
}

write_jaspar_ <- function(motif_list, file_out, matrix_id) {

  if (is.null(matrix_id)) matrix_id <- rep("XX0000", length(motif_list))
  mots <- vector()
  for (i in seq_along(motif_list)) {
    motif <- motif_list[[i]]
    header <- paste0(">", matrix_id[i], " ", motif_slots(motif, "name"))
    mot <- convert_type(motif, "PCM")
    mot <- motif_slots(mot, "motif")
    motA <- paste("A", "[", paste(mot[1, ], collapse = " "), "]")
    motC <- paste("C", "[", paste(mot[2, ], collapse = " "), "]")
    motG <- paste("G", "[", paste(mot[3, ], collapse = " "), "]")
    motT <- paste("T", "[", paste(mot[4, ], collapse = " "), "]")
    mots <- c(header, motA, motC, motG, motT, "")
  }

  writeLines(mots, file_con <- file(file_out))
  close(file_con)

}
