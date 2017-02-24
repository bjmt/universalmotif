## Benjamin Tremblay

# methods for write_meme

# methods:
write_meme.matrix <- function()

write_meme.pwm <- function()

# Taking care of S4 classes:
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
