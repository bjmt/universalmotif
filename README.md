[![Build Status](https://travis-ci.org/bjmt/universalmotif.svg?branch=master)](https://travis-ci.org/bjmt/universalmotif)
[![codecov](https://codecov.io/gh/bjmt/universalmotif/branch/master/graph/badge.svg)](https://codecov.io/gh/bjmt/universalmotif)

# universalmotif #

This package allows for importing most common motif types into R for use by
functions provided by other Bioconductor motif-related packages. Motifs can be 
exported into most major motif formats from various classes as defined by other
Bioconductor packages.

## Installation ##

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("bjmt/universalmotif")
```

