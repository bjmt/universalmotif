[![Build Status](https://travis-ci.org/bjmt/universalmotif.svg?branch=master)](https://travis-ci.org/bjmt/universalmotif) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/bjmt/universalmotif?branch=master&svg=true)](https://ci.appveyor.com/project/bjmt/universalmotif) [![codecov](https://codecov.io/gh/bjmt/universalmotif/branch/master/graph/badge.svg)](https://codecov.io/gh/bjmt/universalmotif) [![Bioc build status](http://bioconductor.org/shields/build/release/bioc/universalmotif.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/universalmotif/) [![Bioc](http://www.bioconductor.org/shields/years-in-bioc/universalmotif.svg)](https://www.bioconductor.org/packages/devel/bioc/html/universalmotif.html#since)
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

### Bioconductor release version

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("universalmotif")
```

### Github development version

```r
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes")
remotes::install_github("bjmt/universalmotif")
```
