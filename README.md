# universalmotif #

This package allows for importing most common motif types into R for use by
functions provided by other Bioconductor motif-related packages. Motifs can be 
exported into most major motif formats from various classes as defined by other
Bioconductor packages.

## Installation ##

```r
devtools::install_github("bjmt/universalmotif")
```

## Currently exported functions ##

  - `read_transfac`
  - `read_jaspar`
  - `read_cisbp`
  - `read_homer`
  - `read_uniprobe`
  - `read_meme`
  - `convert_motifs`
  - `convert_type`
  - `create_motif`
  - `filter_motifs`
  - `motif_dist`
  - `motif_simil`
  - `motif_logo`
  - `motif_tree`
  - `trim_motifs`
  - `motif_rc`
  - `write_transfac`
  - `write_jaspar`
  - `write_meme`
  - `merge_motifs`
  - `write_homer`
  - `align_motifs`

## Other todo ##

  - function documentation
  - vignette
