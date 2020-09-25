library(RxODE)
library(testthat)
test_check("RxODE", stop_on_failure = FALSE, wrap=TRUE, 
           reporter = testthat::LocationReporter)
