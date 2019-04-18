if (R.Version()$arch != "i386") {

  library(testthat)
  library(universalmotif)

  test_check("universalmotif") 

} else TRUE
