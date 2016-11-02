require(RxODE);
context("Test Plot Parsing")
require(digest)

test_that("Emax models are detected",{
    expect_equal(rxSigmoidInfo("Emax*C/(C+E50)"),"Hill(\"Emax\",\"C\",\"E50\",\"1\")");
})

test_that("IDR models are detected",{
    expect_equal(idrInfo("eff","Kin - Kout*(1-C2/(EC50+C2))*eff"), "IDR2(\"Kin\",\"Kout\",\"EC50\",\"eff\",\"C2\",edgeList,nodes)")
})


