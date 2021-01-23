data(JASPAR2018_CORE_DBSCORES, envir = environment())
data(fontDFroboto, envir = environment())

if (is.null(getOption("meme.bin"))) options(meme.bin = "meme")

.onUnload <- function(libpath) {
  library.dynam.unload("universalmotif", libpath)
}
