require(RxODE);
context("Test Sigmoid Parsing")
require(digest)

test_that("Emax models are detected",{
    expect_equal(rxSigmoidInfo("Emax*C/(C+E50)"),"Hill(\"Emax\",\"C\",\"E50\",\"1\")");
})
