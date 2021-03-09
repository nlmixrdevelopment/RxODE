library(RxODE)
library(testthat)
test_check("RxODE", stop_on_failure = FALSE,
           reporter = testthat::LocationReporter)
