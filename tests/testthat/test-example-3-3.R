## Example 3.3 from 
## "Solving Differential Equations in R" by Soetaert et al (2012)
## https://cran.r-project.org/web/packages/diffEq/vignettes/ODEinR.pdf Example #3
library(digest)
rigid <- RxODE("
y1(0)    = 1
y2(0)    = 0
y3(0)    = 0.9
a1       = -2
a2       = 1.25
a3       = -0.5
d/dt(y1) = a1*y2*y3
d/dt(y2) = a2*y1*y3
d/dt(y3) = a3*y1*y2
")

et <- eventTable();
et$add.sampling(seq(0,20,by=0.01))

out <- solve(rigid,et)

test_that("Test rigid body example",{
    expect_equal(digest(signif(as.data.frame(out),6),"sha512"),
                 "425700f73c63dad8a51a0b1f78a1e2dcde2cc64106609fd0528c517a52ce8c918f8ff1ed2bc2deabc5cf007adee345371e27490bb46a69953e6c12586e0b5bf6")
})
