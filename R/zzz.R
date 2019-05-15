data(JASPAR2018_CORE_DBSCORES, envir = environment())
data(JASPAR2018_CORE_DBSCORES_NORM, envir = environment())
data(JASPAR2018_CORE_DBSCORES_2, envir = environment())
data(JASPAR2018_CORE_DBSCORES_NORM_2, envir = environment())

if (is.null(getOption("meme.bin"))) options(meme.bin = "meme")

.onUnload <- function(libpath) {
  library.dynam.unload("universalmotif", libpath)
}
