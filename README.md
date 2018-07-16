# universalmotif #

This package allows for importing most common motif types into R for use by
functions provided by other Bioconductor motif-related packages. Motifs can be 
exported into most major motif formats from various classes as defined by other
Bioconductor packages.

## Installation ##

```r
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
devtools::install_github("bjmt/universalmotif", build_vignettes = TRUE,
                         dependencies = TRUE)
```

