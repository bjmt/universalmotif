#' Create a motif.
#'
#' @param input One of character vector, matrix (PCM, PPM, PWM, or ICM), or
#'              \linkS4class{XStringSet}.
#' @param alphabet Character. 'DNA', 'RNA', 'AA', 'custom', or a combined 
#'                 string representing the letters.
#' @param type Character. 'PCM', 'PPM', 'PWM', or 'ICM'.
#' @param name Character. Motif name.
#' @param pseudoweight Numeric. Correction to be applied to prevent \code{-Inf}
#'                     from apearing in PWM matrices.
#' @param bkg Numeric. Must sum to 1 and be equal in length to the alphabet
#'            length.
#' @param nsites Numeric. Number of sites the motif was constructed from.
#' @param altname Character. Alternate motif name.
#' @param family Character. Transcription factor family.
#' @param organism Character. Species of origin.
#' @param bkgsites Numeric. Total number of sites used to find the motif.
#' @param strand Character. Whether the motif is specific to a certain strand.
#' @param pval Numeric. P-value associated with motif.
#' @param qval Numeric. Adjusted P-value associated with motif.
#' @param eval Numeric. E-value associated with motif.
#' @param extrainfo Character. Any other extra information, represented as
#'                  a named character vector.
#'
#' @return \linkS4class{universalmotif} object.
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
#' # PPM type matrices can also be used as input
#' mat.ppm <- matrix(c(0.1, 0.1, 0.1, 0.1,
#'                     0.5, 0.5, 0.5, 0.5,
#'                     0.1, 0.1, 0.1, 0.1,
#'                     0.3, 0.3, 0.3, 0.3),
#'                     nrow = 4, byrow = TRUE)
#'
#' DNA.motif <- create_motif(mat.ppm, alphabet = "DNA", type = "PPM")
#'
#' ##### create random motifs
#'
#' # these are generated as PPMs with 10 positions
#'
#' DNA.motif <- create_motif()
#' RNA.motif <- create_motif(alphabet = "RNA")
#' AA.motif <- create_motif(alphabet = "AA")
#' custom.motif <- create_motif(alphabet = "QWER")
#'
#' # the number of positions can be specified
#'
#' DNA.motif <- create_motif(5)
#'
#' # if the background frequencies are not provided, they are assumed to be
#' # uniform; if different background frequencies are used, then at each
#' # position \code{abs(rnorm(1, mean = bkg[i], sd = min(c(bkg[i], 1 - bkg[i])))}
#' # is used to generate individual frequencies for each letter
#'
#' DNA.motif <- create_motif(bkg = c(0.3, 0.2, 0.2, 0.3))
#' DNA.motif <- create_motif(10, bkg = c(0.1, 0.4, 0.4, 0.1))
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
#' @param class Desired motif class. Input as 'package-class'. If left empty,
#'              defaults to 'universalmotif-universalmotif'. (See details.)
#' @param BPPARAM See \code{\link[BiocParallel]{bpparam}}.
#'
#' @return Single motif object or list.
#'
#' @details
#'    The following packge-class combinations can be used as input:
#'    \itemize{
#'       \item MotifDb-MotifList
#'       \item TFBSTools-PFMatrix
#'       \item TFBSTools-PWMatrix
#'       \item TFBSTools-ICMatrix
#'       \item TFBSTools-PFMatrixList
#'       \item TFBSTools-PWMatrixList
#'       \item TFBSTools-ICMatrixList
#'       \item seqLogo-pwm
#'       \item motifStack-pcm
#'       \item motifStack-pfm
#'       \item PWMEnrich-PWM
#'       \item motifRG-Motif
#'       \item universalmotif-universalmotif
#'    }
#'
#'    The following package-class combinations can be output:
#'    \itemize{
#'       \item MotIV-pwm2
#'       \item TFBSTools-PFMatrix
#'       \item TFBSTools-PWMatrix
#'       \item TFBSTools-ICMatrix
#'       \item seqLogo-pwm
#'       \item motifStack-pcm
#'       \item motifStack-pfm
#'       \item PWMEnrich-PWM
#'       \item Biostrings-PWM
#'       \item rGADEM-motif
#'    }
#'
#' @examples
#' # convert from universalmotif:
#' jaspar <- read_jaspar(system.file("extdata", "jaspar.txt",
#'                                   package = "universalmotif"))
#' jaspar.motifstack.pfm <- convert_motifs(jaspar, "motifStack-pfm")
#'
#' # convert from another class to universalmotif:
#' library(TFBSTools)
#' data(MA0003.2)
#' motif <- convert_motifs(MA0003.2)
#'
#' # convert from another class to another class
#' motif <- convert_motifs(MA0003.2, "PWMEnrich-PWM")
#'
#' @references
#'    \insertRef{motiv}{universalmotif}
#'
#'    \insertRef{motifdb}{universalmotif}
#'
#'    \insertRef{tfbstools}{universalmotif}
#'
#'    \insertRef{seqlogo}{universalmotif}
#'
#'    \insertRef{motifstack}{universalmotif}
#'
#'    \insertRef{pwmenrich}{universalmotif}
#'
#'    \insertRef{motifrg}{universalmotif}
#'
#'    \insertRef{rgadem}{universalmotif}
#'
#'    \insertRef{biostrings}{universalmotif}
#'
#' @author Benjamin Tremblay, \email{b2tremblay@@uwaterloo.ca}
#' @export
setGeneric("convert_motifs", function(motifs,
                                      class = "universalmotif-universalmotif",
                                      BPPARAM = bpparam())
           standardGeneric("convert_motifs"))
