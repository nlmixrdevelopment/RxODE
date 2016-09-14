library(RxODE)
library(dplyr)
library(digest)
context("rxSolve objects behave as data-frames")


orhs <- RxODE("
x       = 1
a       = 1.0E4
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
", modName="rhs")


et <- eventTable();
et$add.sampling(c(0,4*10^seq(-1,10)))

o1 <- rxSolve(orhs,et);

test_that("rxSolve produces the correct results for rhs",{
    expect_equal(digest(round(as.data.frame(o1),4),"sha512"),
                 "a2554f2b4b71f38b788325f1a3b2c4171a144e5b947ee070dacc732a00e03c76b6baa9e48432300e4b6c281ed87769830ec706b242022a68b7cc9e443b9d319d")
})
