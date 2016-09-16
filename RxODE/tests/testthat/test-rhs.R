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

test_that("rxSolve produces the correct results for rhs",{
    expect_warning(rxSolve(orhs,et),"The initial conditions are at t=0.4 instead of t=0.")
    expect_equal(digest(round(as.data.frame(o1),4),"sha512"),
                 "b22d405c53b2ddba93f0ec55c606447d6ca5ae47c485db6b3db88ee2f056a6386fa92b55b784098d7b79f98bf298158e3b0a34ed128d0fd5e047521c59696178")
})

o1$add.sampling(0);

test_that("rxSolve produces the correct results for rhs",{
    expect_equal(digest(round(as.data.frame(o1),4),"sha512"),
                 "a2554f2b4b71f38b788325f1a3b2c4171a144e5b947ee070dacc732a00e03c76b6baa9e48432300e4b6c281ed87769830ec706b242022a68b7cc9e443b9d319d")
})


rxClean();
