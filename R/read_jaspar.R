#' Import JASPAR motifs.
#'
#' Import JASPAR formatted motifs. See \url{http://jaspar.genereg.net/}.
#' Can be either DNA, RNA, or AA motifs.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#'
#' @references
#'
#' Khan A, Fornes O, Stigliani A, Gheorghe M, Castro-Mondragon JA,
#' van der Lee R, Bessy A, Cheneby J, Kulkarni SR, Tan G, Baranasic
#' D, Arenillas DJ, Sandelin A, Vandepoele K, Lenhard B, Ballester B,
#' Wasserman WW, Parcy F, Mathelier A (2018). “JASPAR 2018: update of
#' the open-access database of transcription factor binding profiles
#' and its web framework.” *Nucleic Acids Research*, **46**, D260-D266.
#'
#' @family read_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @inheritParams read_cisbp
#' @export
read_jaspar <- function(file, skip = 0) {

  # param check --------------------------------------------
  args <- as.list(environment())
  char_check <- check_fun_params(list(file = args$file),
                                 1, FALSE, TYPE_CHAR)
  num_check <- check_fun_params(list(skip = args$skip), 1, FALSE, TYPE_NUM)
  all_checks <- c(char_check, num_check)
  if (length(all_checks) > 0) stop(all_checks_collapse(all_checks))
  #---------------------------------------------------------

  con <- file(file)
  on.exit(close(con), add = TRUE)
  raw_lines <- readLines(con)
  if (skip > 0) raw_lines <- raw_lines[-seq_len(skip)]
  raw_lines <- raw_lines[raw_lines != ""]

  motif_names <- which(grepl("^>", raw_lines))
  if (length(motif_names) == 0) {
    stop("no motifs found in '", file, "'")
  }
  motif_starts <- motif_names + 1
  motif_stops <- c(motif_names[-1] - 1, length(raw_lines))

  if (length(unique(c(length(motif_names), length(motif_starts),
                      length(motif_stops)))) != 1) {
    stop("motifs incorrectly formatted")
  }

  motif_names <- raw_lines[motif_names]
  motif_names <- sub(">", "", motif_names)
  motif_names <- lapply(motif_names, function(x) strsplit(x, "\\s+")[[1]])

  motifs <- mapply(function(x, y) raw_lines[x:y],
                     motif_starts, motif_stops,
                     SIMPLIFY = FALSE)

  get_matrix <- function(x) {
    x <- sub("\\[", "", x)
    x <- sub("\\]", "", x)
    per_line1 <- function(x) {
      x <- strsplit(x, "\\s+")[[1]]
      x <- x[x != ""]
      as.numeric(x[-1])
    }
    per_line2 <- function(x) {
      x <- strsplit(x, "\\s+")[[1]]
      x <- x[x != ""]
      as.numeric(x)
    }
    alphabet <- vapply(x, function(x) strsplit(x, "\\s+")[[1]][1],
                       character(1))
    if (any(alphabet %in% LETTERS)) {
      x2 <- sapply(x, per_line1)
      x2 <- matrix(x2, nrow = length(x), byrow = TRUE)
      rownames(x2) <- alphabet
      x2
    } else {
      x2 <- sapply(x, per_line2)
      matrix(x2, nrow = length(x), byrow = TRUE)
    }
  }

  motifs <- lapply(motifs, get_matrix)

  jaspar2umot <- function(motif, name) {
    if (is.list(motif) || ncol(motif) == 0 || nrow(motif) == 0) {
      stop("empty motif matrix in '", file, "' for motif '", name[1], "'")
    }
    ## Labelled rows may come in any order and give the alphabet; unlabelled
    ## matrices are DNA or, with 20 rows, amino acids.
    if (!is.null(rownames(motif))) {
      motif <- motif[order(rownames(motif)), , drop = FALSE]
      lets <- rownames(motif)
      if (identical(lets, AA_STANDARD2)) {
        alphabet <- "AA"
      } else if (identical(lets, RNA_BASES)) {
        alphabet <- "RNA"
      } else if (identical(lets, DNA_BASES)) {
        alphabet <- "DNA"
      } else alphabet <- collapse_cpp(lets)
    } else if (nrow(motif) == 4) {
      alphabet <- "DNA"
    } else if (nrow(motif) == 20) {
      alphabet <- "AA"
    } else {
      stop("motif '", name[1], "' has ", nrow(motif),
           " unlabelled rows; label the rows to give its alphabet")
    }
    mot <- universalmotif_cpp(name = name[1], altname = name[2],
                   type = "PCM", alphabet = alphabet,
                   nsites = max(colSums(motif)),
                   motif = motif)
    validObject_universalmotif(mot)
    mot
  }

  motifs <- mapply(jaspar2umot, motifs, motif_names, 
                     SIMPLIFY = FALSE)

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs

}
