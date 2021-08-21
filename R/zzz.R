data(JASPAR2018_CORE_DBSCORES, envir = environment())
data(fontDFroboto, envir = environment())

if (is.null(getOption("meme.bin"))) options(meme.bin = "meme")
if (is.null(getOption("pseudocount.warning"))) options(pseudocount.warning = TRUE)
if (is.null(getOption("universalmotif_df.warning"))) options(universalmotif_df.warning = TRUE)

.onUnload <- function(libpath) {
  library.dynam.unload("universalmotif", libpath)
}
