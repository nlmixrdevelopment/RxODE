Sys.setenv("R_TESTS" = "")
library(testthat)
library(RxODE)
test_check("RxODE")
