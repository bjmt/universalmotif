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
      + more formats still left to add!
  - `convert_motifs`
  - `convert_type`
  - `create_motif`
  - `motif_slots`
  - `seqLogo-method`

## Planned exported functions ##

  - `run_meme`
  - `run_tomtom`
  - `run_homer`
  - `run_stamp`
  - `write_motifs`
  - `filter_motifs`
  - `motifStack-method`

## Other todo ##

  - lots of code cleanup to be done!
