data(JASPAR2018_CORE_DBSCORES, envir = environment())
data(fontDFroboto, envir = environment())

if (is.null(getOption("meme.bin"))) options(meme.bin = "meme")
if (is.null(getOption("pseudocount.warning"))) options(pseudocount.warning = TRUE)
if (is.null(getOption("universalmotif_df.warning"))) options(universalmotif_df.warning = TRUE)
if (is.null(getOption("universalmotif.suggest.scan_sequences_lite"))) options(universalmotif.suggest.scan_sequences_lite = TRUE)
if (is.null(getOption("universalmotif.suggest.compare_motifs_lite")))  options(universalmotif.suggest.compare_motifs_lite  = TRUE)
if (is.null(getOption("universalmotif.suggest.enrich_motifs_lite")))   options(universalmotif.suggest.enrich_motifs_lite   = TRUE)
if (is.null(getOption("universalmotif.suggest.merge_motifs_lite")))    options(universalmotif.suggest.merge_motifs_lite    = TRUE)
if (is.null(getOption("universalmotif.suggest.merge_similar_lite")))   options(universalmotif.suggest.merge_similar_lite   = TRUE)
if (is.null(getOption("universalmotif.suggest.motif_tree_lite")))      options(universalmotif.suggest.motif_tree_lite      = TRUE)

.onUnload <- function(libpath) {
  library.dynam.unload("universalmotif", libpath)
}
