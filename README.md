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
      + more formats still left to add
  - `convert_motifs`
  - `convert_type`
  - `create_motif`
      + go back and use `Biostrings::consensusMatrix` for `DNAStringSet` input?
  - `motif_slots`
  - `seqLogo`
  - `filter_motifs`
  - `trim_motifs`
  - `motifStack`
  - `write_meme`

## Planned exported functions ##

  - `run_meme`
  - `run_tomtom`
  - `run_homer`
  - `run_stamp`
  - `write_motifs`
    + homer, jaspar, transfac, uniprobe
  - `PWMSimilarity` (`TFBSTools`)
  - `motifDistances`, `motifHclust`, `motifCutree` (`MotIV`)
  - `cluster_motifs`

## Other todo ##

  - lots of code cleanup to be done! (mostly for the read functions)
  - to check platform: run `.Platform` for a list of info
      + for platform-specific code, have the following line:
      + `stopifnot(.Platform$OS.type == "unix")`
      + (the only other possibility is "windows")
  - I have tests for `read_motifs`, `convert_type`, and `create_motif`;
    but I'm unsure as to whether I should write tests for `convert_motifs`
    considering it would require the respective packages being loaded..
        + edit: `testthat` can allow for tests not to be run via `R CMD CHECK`
  - add `GRanges` objects and BED format compatibility (i.e. going straight
    from these formats to MEME)
  - motif similarity clustering + heatmap generation
      + combine `PWMSimilarity`, `outer` and `pheatmap`
  - get rid of filter functionality in read functions, leave that to
    `filter_motifs`
