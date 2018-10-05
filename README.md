[![Build Status](https://travis-ci.org/bjmt/universalmotif.svg?branch=master)](https://travis-ci.org/bjmt/universalmotif) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/bjmt/universalmotif?branch=master&svg=true)](https://ci.appveyor.com/project/bjmt/universalmotif) [![codecov](https://codecov.io/gh/bjmt/universalmotif/branch/master/graph/badge.svg)](https://codecov.io/gh/bjmt/universalmotif)
# universalmotif

This package allows for importing most common motif types into R for use by
functions provided by other Bioconductor motif-related packages. Motifs can be 
exported into most major motif formats from various classes as defined by other
Bioconductor packages. Furthermore, this package allows for easy manipulation
of motifs, such as creation, trimming, shuffling, P-value calculations,
filtering, type conversion, reverse complementation, alphabet switching, random
motif site generation, and comparison. Alongside are also included functions
for interacting with sequences, such as motif scanning and enrichment, as well
as sequence creation and shuffling functions. Finally, this package implements
higher-order motifs, allowing for more accurate sequence scanning and motif
enrichment.

## Installation

If you'd like to build the vignettes locally, a few prerequisite packages not
listed in the 'Imports' section are needed:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install(c("BiocStyle", "TFBSTools", "motifStack", "Logolas",
                       "MotifDb"))
BiocManager::install("bjmt/universalmotif")
```

To install the package with pre-built vignettes:

  - Go to the 'release' tab
  - Download the latest 'universalmotif_X.X.X.tar.gz' file
  - Run `R CMD INSTALL universalmotif_X.X.X.tar.gz`
