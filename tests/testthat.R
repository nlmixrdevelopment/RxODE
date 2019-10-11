options(RxODE.rtools.setup=TRUE) # Assume Rtools is setup
Sys.setenv("R_TESTS" = "")
library(testthat)
library(RxODE)
test_check("RxODE")
