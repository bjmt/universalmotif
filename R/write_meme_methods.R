######################################################################
## Benjamin Tremblay
##
## methods for write_meme
##
######################################################################

# make sure this is loaded after respective functions:

#' @include write_functions.R
#' @export
write_meme.matrix <- function(motif_list, file_out, version = 4,
                           bkg = c(0.25, 0.25, 0.25, 0.25), strands = "+ -",
                           alphabet = NULL) {
  print("hello")
}

#' @export
write_meme.pwm <- function(motif_list, file_out, version = 4,
                           bkg = c(0.25, 0.25, 0.25, 0.25), strands = "+ -",
                           alphabet = NULL) {
  print("bye")
}

# Taking care of S4 classes:

#' @importClassesFrom seqLogo pwm
#' @export
setMethod("write_meme", "pwm", write_meme.pwm)

# testing:
# tests4 <- read_meme("~/universalmotif/inst/extdata/minimal.meme", out_format = "seqLogo-pwm")[[1]]
# tests3 <- read_meme("~/universalmotif/inst/extdata/minimal.meme")[[1]]
#
# testing <- function(test) UseMethod("testing")  # define it for use with S3 objects
# setGeneric("testing")  # define it for use with S4 objects
#
# # methods:
# testing.matrix <- function(test) print("detected matrix")
# testing.pwm <- function(test) print("detected pwm")
#
# # for S4, testing.pwm must be loaded before it can be assigned as a method
# setMethod("testing", "pwm", testing.pwm)  # assign function as method to S4 object
