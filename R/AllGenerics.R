#' Create a motif.
#'
#' @param input One of character vector, matrix (PCM, PPM, PWM, or ICM), or
#'              XStringSet.
#' @param alphabet Character.
#' @param type Character.
#' @param name Character.
#' @param pseudoweight Numeric.
#' @param bkg Numeric.
#' @param nsites Numeric.
#' @param altname Character.
#' @param family Character.
#' @param organism Character.
#' @param bkgsites Numeric.
#' @param strand Character.
#' @param pval Numeric.
#' @param qval Numeric.
#' @param eval Numeric.
#' @param extrainfo Character.
#'
#' @return Motif object.
#'
#' @examples
#' ##### create motifs from a single string
#'
#' # motif is initially generated as a PCM; change final type as desired
#' DNA.motif <- create_motif("TATAWAW", type = "PPM", pseudoweight = 0)
#'
#' # nsites will be set to the number of input sequences unless specified 
#' DNA.motif <- create_motif("TTTTTTT", nsites = 10)
#'
#' # if ambiguity letters are found and nsites is not specified, nsites will
#' # be set to the minimum required to respect amibiguity letters
#' DNA.motif <- create_motif("TATAWAW")
#' DNA.motif <- create_motif("NNVVWWAAWWDDN")
#'
#' # be careful about setting nsites when using ambiguity letters!
#' DNA.motif <- create_motif("NNVVWWAAWWDDN", nsites = 1)
#'
#' RNA.motif <- create_motif("UUUCCG")
#' 
#' # 'create_motif' will try to detect the alphabet type; this can be 
#' # unreliable for AA and custom alphabets as DNA and RNA alphabets are
#' # detected first
#' AA.motif <- create_motif("AVLK", alphabet = "AA")
#' 
#' custom.motif <- create_motif("QWER", alphabet = "custom")
#' # specify custom alphabet
#' custom.motif <- create_motif("QWER", alphabet = "QWERASDF")
#'
#' ###### create motifs from multiple strings of equal length
#'
#' DNA.motif <- create_motif(c("TTTT", "AAAA", "AACC", "TTGG"), type = "PPM")
#' DNA.motif <- create_motif(c("TTTT", "AAAA", "AACC", "TTGG"), nsites = 20)
#' RNA.motif <- create_motif(c("UUUU", "AAAA", "AACC", "UUGG"), type = "PWM")
#' AA.motif <- create_motif(c("ARNDCQ", "EGHILK", "ARNDCQ"), alphabet = "AA")
#' custom.motif <- create_motif(c("POIU", "LKJH", "POIU", "CVBN"),
#'                              alphabet = "custom")
#'
#' # ambiguity letters are only allowed for single consensus strings; the
#' # following fails
#' \dontrun{
#' create_motif(c("WWTT", "CCGG"))
#' create_motif(c("XXXX", "XXXX"), alphabet = "AA")
#' }
#'
#' ##### create motifs from XStringSet objects
#'
#' library(Biostrings)
#' 
#' DNA.set <- DNAStringSet(c("TTTT", "AAAA", "AACC", "TTGG"))
#' DNA.motif <- create_motif(DNA.set)
#' RNA.set <- RNAStringSet(c("UUUU", "AACC", "UUCC"))
#' RNA.motif <- create_motif(RNA.set)
#' AA.set <- AAStringSet(c("VVVLLL", "AAAIII"))
#' AA.motif <- create_motif(AA.set)
#' 
#' # custom motifs can be created from BStringSet objects
#' B.set <- BStringSet(c("QWER", "ASDF", "ZXCV", "TYUI"))
#' custom.motif <- create_motif(B.set)
#'
#' ##### create motifs from matrices
#'
#' mat <- matrix(c(1, 1, 1, 1,
#'                 2, 0, 2, 0,
#'                 0, 2, 0, 2,
#'                 0, 0, 0, 0),
#'                 nrow = 4, byrow = TRUE)
#' DNA.motif <- create_motif(mat, alphabet = "DNA")
#' RNA.motif <- create_motif(mat, alphabet = "RNA", nsites = 20)
#' custom.motif <- create_motif(mat)
#'
#' # specify custom alphabet
#' custom.motif <- create_motif(mat, alphabet = "QWER")
#'
#' # alphabet can be detected from rownames
#' rownames(mat) <- DNA_BASES
#' DNA.motif <- create_motif(mat)
#' rownames(mat) <- c("Q", "W", "E", "R")
#' custom.motif <- create_motif(mat)
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("create_motif", function(input, alphabet, type = "PCM",
                                    name = "motif", pseudoweight = 0.8,
                                    bkg, nsites, altname, family,
                                    organism, bkgsites, strand, pval, qval,
                                    eval, extrainfo)
           standardGeneric("create_motif"))

#' Convert motif class.
#'
#' @param motifs Single motif object or list.
#' @param class Desired motif class. E.g. 'motifStack-pfm'.
#' @param BPPARAM Param for bplapply.
#'
#' @return Single motif object or list.
#'
#' @examples
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.motifstack.pfm <- convert_motifs(jaspar, "motifStack-pfm")
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("convert_motifs", function(motifs,
                                      class = "universalmotif-universalmotif",
                                      BPPARAM = bpparam())
           standardGeneric("convert_motifs"))
