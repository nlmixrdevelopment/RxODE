Sys.setenv("R_TESTS" = "")
library(testthat)
library(RxODE)
rxClean();
test_check("RxODE")
