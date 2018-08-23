data(JASPAR2018_CORE_DBSCORES, envir = environment())
data(JASPAR2018_CORE_DBSCORES_NORM, envir = environment())

.onUnload <- function(libpath) {
  library.dynam.unload("universalmotif", libpath)
}
