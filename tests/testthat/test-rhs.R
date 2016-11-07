library(RxODE)
library(dplyr)
library(digest)
context("rxSolve right handed differental equations")


orhs <- RxODE("
x(0)    = 1
a       = 1.0E4
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
")


et <- eventTable();
et$add.sampling(c(4*10^seq(-1,10)))

o1 <- rxSolve(orhs,et);

print(o1);
test_that("rxSolve produces the correct results for rhs",{
    expect_warning(rxSolve(orhs,et),"The initial conditions are at t = 0.4 instead of t = 0.")
    expect_equal(digest(round(as.data.frame(o1),4)),
                 "73dd8f24d60cf8cec553cd8030b26707")
})

o1$add.sampling(0);
print(o1);
test_that("rxSolve produces the correct results for rhs",{
    expect_equal(digest(round(as.data.frame(o1),4)),
                 "f8594b891db1162357d2ccd93c687ed2")
})

rxClean();
