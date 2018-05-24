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

  - `read_motifs`
  - `convert_motifs`
  - `convert_type`
  - `create_motif`
  - `filter_motifs`
  - `trim_motifs`

## Planned exported functions ##

  - `run_meme`
  - `run_tomtom`
  - `run_homer`
  - `run_stamp`
  - `write_motifs`

## Other todo ##

  - to check platform: run `.Platform` for a list of info
      + for platform-specific code, have the following line:
      + `stopifnot(.Platform$OS.type == "unix")`
      + (the only other possibility is "windows")
  - add `GRanges` objects and BED format compatibility (i.e. going straight
    from these formats to MEME)
