data(JASPAR2018_CORE_DBSCORES, envir = environment())
data(fontDFroboto, envir = environment())

if (is.null(getOption("meme.bin"))) options(meme.bin = "meme")
if (is.null(getOption("pseudocount.warning"))) options(pseudocount.warning = TRUE)
if (is.null(getOption("universalmotif_df.warning"))) options(universalmotif_df.warning = TRUE)
if (is.null(getOption("universalmotif.suggest.scan_sequences2"))) options(universalmotif.suggest.scan_sequences2 = TRUE)
if (is.null(getOption("universalmotif.suggest.compare_motifs2")))  options(universalmotif.suggest.compare_motifs2  = TRUE)

.onUnload <- function(libpath) {
  library.dynam.unload("universalmotif", libpath)
}
