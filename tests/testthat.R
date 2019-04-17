library(testthat)
library(universalmotif)

if (R.Version()$arch != "i386") test_check("universalmotif") else TRUE
