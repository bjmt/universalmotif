data(JASPAR2018_CORE_DBSCORES, envir = environment())

.onUnload <- function(libpath) {
  library.dynam.unload("universalmotif", libpath)
}
