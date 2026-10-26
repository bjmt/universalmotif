#' Import TRANSFAC motifs.
#'
#' Import TRANSFAC formatted motifs. Assumed to be DNA motifs, type PCM.
#' See `system.file("extdata", "transfac.txt", pacakge="universalmotif")`
#' for an example motif.
#'
#' @return `list` [universalmotif-class] objects.
#'
#' @details
#'
#' A few TRANSFAC tags are recognized, including AC, ID, NA, HC and OS.
#' HC will be set to the family slot and OS to the organism slot.
#' If AC, ID and NA are present, then AC will be set as the motif name
#' and NA as the alternate name. If AC is absent, then ID is set as the
#' name. If ID is also absent, then NA is set as the motif name.
#'
#' @examples
#' transfac <- read_transfac(system.file("extdata", "transfac.txt",
#'                                       package = "universalmotif"))
#' 
#' @references
#'
#' Wingender E, Dietze P, Karas H, Knuppel R (1996). “TRANSFAC: A
#' Database on Transcription Factors and Their DNA Binding Sites.”
#' *Nucleic Acids Research*, **24**, 238-241.
#'
#' @family read_motifs
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@@uwaterloo.ca}
#' @inheritParams read_cisbp
#' @export
read_transfac <- function(file, skip = 0) {

  # From https://biopython.readthedocs.io/en/latest/chapter_motifs.html:
  #
  # AC	Accession number
  # AS	Accession numbers, secondary
  # BA	Statistical basis
  # BF	Binding factors
  # BS	Factor binding sites underlying the matrix
  # CC	Comments
  # CO	Copyright notice
  # DE	Short factor description
  # DR	External databases
  # DT	Date created/updated
  # HC	Subfamilies
  # HP	Superfamilies
  # ID	Identifier
  # NA	Name of the binding factor
  # OC	Taxonomic classification
  # OS	Species/Taxon
  # OV	Older version
  # PV	Preferred version
  # TY	Type
  # XX	Empty line; these are not stored in the Record.
  #
  # RN	Reference number
  # RA	Reference authors
  # RL	Reference data
  # RT	Reference title
  # RX	PubMed ID

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
  raw_lines <- raw_lines[!grepl("^XX", raw_lines)]
  raw_lines <- raw_lines[raw_lines != ""]

  motif_stops <- grep("^//", raw_lines)
  if (length(motif_stops) == 0) {
    stop("no motifs found in '", file, "' (no '//' separators)")
  }

  if (1 == motif_stops[1]) {
    motif_stops <- motif_stops[-1]
    raw_lines <- raw_lines[-1]
  }
  if (length(motif_stops) == 0) {
    stop("no motifs found in '", file, "'")
  }
  motif_starts <- c(1, motif_stops[-length(motif_stops)] + 1)
  motif_stops <- motif_stops - 1

  motifs <- mapply(function(x, y) raw_lines[x:y],
                     motif_starts, motif_stops,
                     SIMPLIFY = FALSE)

  ## Records without a matrix, such as the VV header block of TRANSFAC
  ## releases, are skipped.
  motifs <- motifs[vapply(motifs, function(x) any(grepl("^P[0O]", x)), logical(1))]
  if (length(motifs) == 0) {
    stop("no motifs found in '", file, "'")
  }

  get_matrix <- function(x) {
    p0 <- which(grepl("^P[0O]", x))[1]
    header <- strsplit(trimws(x[p0]), "\\s+")[[1]][-1]
    x <- x[-seq_len(p0)]
    mot_i <- vapply(x, function(x) {
                         x <- strsplit(x, "\\s+")[[1]][1]
                         x <- suppressWarnings(as.numeric(x))
                         !is.na(x)
                       }, logical(1))
    mot <- x[mot_i]
    per_line <- function(x) {
      x <- strsplit(x, "\\s+")[[1]]
      as.numeric(x[2:5])
    }
    mot <- vapply(mot, per_line, numeric(4))
    mot <- matrix(mot, ncol = 4, byrow = TRUE)
    ## Columns follow the letters of the P0 line when it names the bases.
    if (length(header) == 4 && identical(sort(header), DNA_BASES))
      mot <- mot[, match(DNA_BASES, header), drop = FALSE]
    mot
  }

  motif_matrix <- lapply(motifs, get_matrix)

  ## The name is the accession if there is one, with the identifier or
  ## factor name as the alternate name, whatever order the fields come in.
  parse_meta <- function(x) {
    key <- sub("\\s.*$", "", x)
    value <- trimws(sub("^\\S+", "", x))
    get <- function(k) {
      v <- value[key == k & value != ""]
      if (length(v)) v[1] else character()
    }
    ac <- get("AC")
    id <- get("ID")
    na <- get("NA")
    if (length(ac) && length(id)) {
      name <- ac; altname <- id
    } else if (length(ac) && length(na)) {
      name <- ac; altname <- na
    } else if (length(id) && length(na)) {
      name <- id; altname <- na
    } else {
      name <- c(ac, id, na); altname <- character()
    }
    list(name = if (length(name)) name[1] else "motif", altname = altname,
         family = get("HC"), organism = get("OS"))
  }

  motif_meta <- lapply(motifs, parse_meta)

  nsites <- lapply(motif_matrix, function(x) max(rowSums(x)))
  # Whole counts with the same sum at every position are kept exactly,
  # since converting back from probabilities takes an nsites of 1 for 100.
  counts <- lapply(motif_matrix, function(x) {
    cs <- rowSums(x)
    if (all(x == round(x)) && all(cs == cs[1])) t(x) else NULL
  })
  motif_matrix <- lapply(motif_matrix,
                         function(x) {
                           cs <- rowSums(x)
                           for (i in seq_along(cs)) {
                             x[i, ] <- x[i, ] / cs[i]
                           }
                           x
                         })

  motifs <- mapply(function(x, y, z) {
                      mot <- universalmotif_cpp(name = y$name,
                                     altname = if (length(y$altname)) y$altname else NA_character_,
                                     family = if (length(y$family)) y$family else NA_character_,
                                     organism = if (length(y$organism)) y$organism else NA_character_,
                                     motif = t(x),
                                     alphabet = "DNA",
                                     nsites = z,
                                     type = "PPM")
                      validObject_universalmotif(mot)
                      mot
                     }, motif_matrix, motif_meta, nsites)

  motifs <- mapply(function(m, x) {
    m <- convert_type_internal(m, "PCM")
    if (!is.null(x)) {
      dimnames(x) <- dimnames(m@motif)
      m@motif <- x
      validObject_universalmotif(m)
    }
    m
  }, motifs, counts)

  if (length(motifs) == 1) motifs <- motifs[[1]]
  motifs

}
