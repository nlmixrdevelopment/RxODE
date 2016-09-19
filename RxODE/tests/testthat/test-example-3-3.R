## Example 3.3 from 
## "Solving Differential Equations in R" by Soetaert et al (2012)
library(digest)
rigid <- RxODE("
y1(0)    = 1
y2(0)    = 0
y3(0)    = 0.91
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
                 "40a260ff675b9ed41e1b677560a1b99b2f63e92d83bbdcaa4a9ca2a8063e1b2d60967c1352ca98f739df52baa723d1823cfd27f35f4eb3617a1afffa165ef9a6")
})
