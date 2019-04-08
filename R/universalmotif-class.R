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
#' @slot bkg `numeric` 0-order probabilities must be provided for all letters.
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

  slots = list(

    name = "character",
    altname = "character",
    family = "character",
    organism = "character",
    motif = "matrix",
    alphabet = "character",
    type = "character",
    icscore = "numeric",
    nsites = "numeric",
    pseudocount = "numeric",
    bkg = "numeric",
    bkgsites = "numeric",
    consensus = "character",
    strand = "character",
    pval = "numeric",
    qval = "numeric",
    eval = "numeric",
    multifreq = "list",
    extrainfo = "character"

  )

)

setValidity("universalmotif", function(object) {

  msg <- validObject_universalmotif(object, FALSE)

  if (length(msg) == 0) TRUE else msg

})
