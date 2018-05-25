# universalmotif #

This package allows for importing most common motif types into R for use by
functions provided by other Bioconductor motif-related packages. Motifs can be
filtered during import based on parameters provided by each motif formats, as
well as during export. Motifs can be exported into most major motif formats from
various classes as defined by other Bioconductor packages.

## Installation ##

```r
devtools::install_github("bjmt/universalmotif")
```

## Currently exported functions ##

  - `read_transfac`
  - `read_jaspar`
  - `read_cisbp`
  - `convert_motifs`
  - `convert_type`
  - `create_motif`
  - `filter_motifs`
  - `motif_dist`
  - `motif_simil`
  - `motif_logo`
  - `motif_tree`
  - `trim_motifs`

## Planned exported functions ##

  - `run_meme`
  - `run_homer`
  - `write_meme`
  - `write_transfac`
  - `write_uniprobe`
  - `write_homer`
  - `write_jaspar`
  - `write_cisbp`
  - `read_homer`
  - `read_uniprobe`
  - `read_meme`

## Other todo ##

  - function documentation
  - vignette
  - to check platform: run `.Platform` for a list of info
      + for platform-specific code, have the following line:
      + `stopifnot(.Platform$OS.type == "unix")`
      + (the only other possibility is "windows")
  - add `GRanges` objects and BED format compatibility (i.e. going straight
    from these formats to MEME)
