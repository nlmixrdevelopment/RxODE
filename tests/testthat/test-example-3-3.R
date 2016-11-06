## Example 3.3 from 
## "Solving Differential Equations in R" by Soetaert et al (2012)
## https://cran.r-project.org/web/packages/diffEq/vignettes/ODEinR.pdf Example #3
library(digest)
context("Example 3.3");
rxClean();
rigid.txt <- "
y1(0)    = 1
y2(0)    = 0
y3(0)    = 0.9
a1       = -2
a2       = 1.25
a3       = -0.5
d/dt(y1) = a1*y2*y3
d/dt(y2) = a2*y1*y3
d/dt(y3) = a3*y1*y2
";
rigid <- RxODE(rigid.txt)

et <- eventTable();
et$add.sampling(seq(0,20,by=0.01))

out <- solve(rigid,et)

test_that("Test rigid body example",{
    expect_equal(digest(round(as.data.frame(out),3)),
                 "a3f42a57944330af983e9b78c3a69a30")
})

test_that("Different solves give same results",{
    out2 <- solve(rigid$cmpMgr,et);
    expect_equal(out,out2)
    out2 <- solve(rigid$cmpMgr$rxDll(),et);
    expect_equal(out,out2)
    out2 <- solve(rigid.txt,et);
    expect_equal(as.data.frame(out),as.data.frame(out2));
})

rxClean();
