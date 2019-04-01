#' universalmotif: Motif class.
#'
#' Container for motif objects. See [create_motif()] for creating
#' motifs as well as a more detailed description of the slots. For a
#' brief description of available methods, see `examples`.
#'
#' @slot name `character(1)`
#' @slot altname `character(1)`
#' @slot family `character(1)`
#' @slot organism `character(1)`
#' @slot motif `matrix`
#' @slot alphabet `character(1)`
#' @slot type `character(1)`
#' @slot icscore `numeric(1)` Generated automatically.
#' @slot nsites `numeric(1)`
#' @slot pseudocount `numeric(1)`
#' @slot bkg `numeric` Length equal to number of letters in alphabet.
#' @slot bkgsites `numeric(1)`
#' @slot consensus `character` Generated automatically.
#' @slot strand `character(1)`
#' @slot pval `numeric(1)`
#' @slot qval `numeric(1)`
#' @slot eval `numeric(1)`
#' @slot multifreq `list`
#' @slot extrainfo `character`
#'
#' @return A motif object of class [universalmotif-class].
#'
#' @examples
#' ## [
#' ## Access the slots.
#' motif <- create_motif()
#' motif["motif"]
#' # you can also access multiple slots at once, released as a list
#' motif[c("motif", "name")]
#'
#' ## [<-
#' ## Replace the slots.
#' motif["name"] <- "new name"
#' # some slots are protected
#' # motif["consensus"] <- "AAAA"  # not allowed
#'
#' ## c
#' ## Assemble a list of motifs.
#' c(motif, motif)
#'
#' ## as.data.frame
#' ## Represent a motif as a data.frame. The actual motif matrix is lost.
#' ## Necessary for `summarise_motifs`.
#' as.data.frame(motif)
#'
#' ## subset
#' ## Subset a motif matrix by column.
#' subset(motif, 3:7)  # extract motif core
#'
#' ## normalize
#' ## Apply the pseudocount slot (or `1`, if the slot is set to zero) to the
#' ## motif matrix.
#' motif2 <- create_motif("AAAAA", nsites = 100, pseudocount = 1)
#' normalize(motif2)
#'
#' ## rowMeans
#' ## Calculate motif rowMeans.
#' rowMeans(motif)
#'
#' ## colMeans
#' ## Calculate motif colMeans.
#' colMeans(motif)
#'
#' ## colSums
#' ## Calculate motif colSums
#' colSums(motif)
#'
#' ## rowSums
#' ## Calculate motif rowSums.
#' rowSums(motif)
#'
#' ## nrow
#' ## Count motif rows.
#' nrow(motif)
#'
#' ## ncol
#' ## Count motif columns.
#' ncol(motif)
#'
#' ## colnames
#' ## Get motif colnames.
#' colnames(motif)
#'
#' ## rownames
#' ## Get motif rownames.
#' rownames(motif)
#'
#' ## cbind
#' ## Bind motifs together to create a new motif.
#' cbind(motif, motif2)
#'
#' @author Benjamin Jean-Marie Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @name universalmotif-class
#' @rdname universalmotif-class
#' @exportClass universalmotif
universalmotif <- setClass("universalmotif",
                           slots = list(name = "character",
                                        altname = "character",
                                        family = "character",
                                        organism = "character",
                                        motif = "matrix",
                                        alphabet = "character",
                                        type = "character", icscore = "numeric",
                                        nsites = "numeric",
                                        pseudocount = "numeric",
                                        bkg = "numeric",
                                        bkgsites = "numeric",
                                        consensus = "character",
                                        strand = "character", pval = "numeric",
                                        qval = "numeric", eval = "numeric",
                                        multifreq = "list",
                                        extrainfo = "character"))

setValidity("universalmotif",
            function(object) {

            msg <- vector()
            valid <- TRUE

            ## character slot length checks

            # mandatory
            char_checks1 <- c("name", "type", "alphabet", "strand")
            # optional
            char_checks2 <- c("altname", "family", "consensus", "organism")

            char_results1 <- vapply(char_checks1,
                                   function(x) length(object[x]) == 1,
                                   logical(1))
            char_results2 <- vapply(char_checks1,
                                   function(x) length(object[x]) <= 1,
                                   logical(1))

            char_checks_all <- c(char_checks1, char_checks2)
            char_results_all <- c(char_results1, char_results2)

            if (any(isFALSE(char_results_all))) {
              valid <- FALSE
              for (check in char_checks_all[isFALSE(char_results_all)]) {
                msg <- c(msg, paste0("motif '", check,
                                     "' must be a character vector of length 1"))
              }
            }

            ## character slot specific checks

            # alphabets <- c("DNA", "RNA", "AA", "custom")
            types <- c("PCM", "PPM", "PWM", "ICM")
            strands <- c("+", "-", "+-", "-+")

            if (!object["type"] %in% types) {
              valid <- FALSE
              msg <- c(msg, "motif 'type' must be either 'PCM', 'PPM', 'PWM', or 'ICM'")
            }

            if (!object["strand"] %in% strands) {
              valid <- FALSE
              msg <- c(msg, "motif 'strand' must be either '+', '-', or '+-'")
            }

            ## numeric slot length checks

            num_checks <- c("pval", "qval", "eval", "bkgsites", "pseudocount",
                            "nsites")
            num_results <- vapply(num_checks,
                                  function(x) length(object[x]) <= 1,
                                  logical(1))

            if (any(isFALSE(num_results))) {
              valid <- FALSE
              for (check in num_checks[isFALSE(num_checks)]) {
                msg <- c(msg, paste0("motif '", check,
                                     "' must be a numeric vector of length 1"))
              }
            }

            ## bkg slot check

            bkg_check1 <- object["bkg"]
            bkg_check2 <- nrow(object["motif"])
            if (length(bkg_check1) > 0) {
              if (length(bkg_check1) != bkg_check2) {
                valid <- FALSE
                msg <- c(msg, "motif 'bkg' must be a numeric vector of length equal to the number of letters in motif")
              }
            }

            ## motif slot check

            mat_type <- object["type"]
            mat <- object["motif"]
            mat_colsums <- colSums(mat)

            # PCM
            if (mat_type == "PCM") {
              mat_nsites <- object["nsites"]
              if (length(unique(mat_colsums)) > 1) {
                warning("not all positions have identical count totals")
                # if (sd(unique(mat_colsums)) > 0.1) {
                  # valid <- FALSE
                  # msg <- c(msg, "motif of type 'PCM' must have equal column sums")
                # }
              }
              if (length(mat_nsites) > 0) {
                if (length(unique(mat_colsums)) == 1) {
                  if (unique(mat_colsums) != mat_nsites) {
                    valid <- FALSE
                    msg <- c(msg, "motif of type 'PCM' must have column sums equal to 'nsites'")
                  }
                }
              }
            }

            # PPM
            if (mat_type == "PPM") {
              if (any(mat_colsums > 1.01) || any(mat_colsums < 0.99)) {
                valid <- FALSE
                msg <- c(msg, "motif of type 'PPM' must have column sums equal to 1")
              }
              if (any(mat < 0)) {
                valid <- FALSE
                msg <- c(msg, "motif of type 'PPM' can only have positive probabilities")
              }
            }

            # check it matches alphabet
            alph <- object["alphabet"]
            if (alph %in% c("DNA", "RNA")) {
              if (nrow(mat) != 4) {
                valid <- FALSE
                msg <- c(msg,
                         "motif with 'alphabet' of type 'DNA' or 'RNA' can only have 4 rows")
              }
              if (alph == "DNA" && any(rownames(object["motif"]) != DNA_BASES)) {
                valid <- FALSE
                msg <- c(msg, "DNA motif must have rownames A, C, G, T")
              }
              if (alph == "RNA" && any(rownames(object["motif"]) != RNA_BASES)) {
                valid <- FALSE
                msg <- c(msg, "RNA motif must have rownames A, C, G, U")
              }
            } else if (alph == "AA") {
              if (nrow(mat) != 20) {
                valid <- FALSE
                msg <- c(msg,
                         "motif with 'alphabet' of type 'AA' can only have 20 rows")
              }
              if (any(!rownames(object["motif"]) %in% AA_STANDARD)) {
                valid <- FALSE
                msg <- c(msg, paste("AA motif must have rownames",
                                    paste(sort(AA_STANDARD), collapse = "")))
              }
            } else if (alph != "custom") {
              if (nrow(mat) != length(safeExplode(alph))) {
                valid <- FALSE
                msg <- c(msg, paste0("motif with alphabet '",
                                     "' has an incorrect number of rows"))
              }
              if (any(!rownames(object["motif"]) %in% safeExplode(alph))) {
                valid <- FALSE
                msg <- c(msg, "motif rownames must match alphabet")
              }
            }

            # consensus
            consensus <- object["consensus"]
            if (length(consensus) != 0) {
              consensus <- safeExplode(consensus)
              if (length(consensus) != ncol(object["motif"])) {
                valid <- FALSE
                msg <- c(msg, "consensus string length does not match motif length")
              }
            }

            # check multifreq
            # multifreq <- object["multifreq"]
            # if (length(multifreq) > 0) {
              # for (i in seq_along(names(multifreq))) {
                # if (!is.matrix(multifreq[[i]])) {
                  # valid <- FALSE
                  # msg <- c(msg, "multifreq slot must be a list of matrices")
                # }
                # if (as.integer(i) != as.numeric(i) || is.na(as.numeric(i))) {
                  # valid <- FALSE
                  # msg <- c(msg, "multifreq list names must be whole numbers")
                # }
                # if (nrow(multifreq[[i]]) != 4^as.numeric(i)) {
                  # valid <- FALSE
                  # msg <- c(msg, paste0(i, "-letter multifreq matrix must have ",
                                       # 4^as.numeric(i), " rows"))
                # }
                # mot_len_check <- ncol(object["motif"]) - as.numeric(i) + 1
                # if (ncol(multifreq[[i]]) != mot_len_check) {
                  # valid <- FALSE
                  # msg <- c(msg, paste0(i, "-letter multifreq matrix must have ",
                                       # mot_len_check, " columns"))
                # }
              # }
            # }

            if (valid) TRUE else msg

            })
