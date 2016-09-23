## Example 3.2 from 
## "Solving Differential Equations in R" by Soetaert et al (2012)
## https://cran.r-project.org/web/packages/diffEq/vignettes/ODEinR.pdf Example #2
## Lorenz model in Vingette.

library("RxODE")
library("digest")

ode <- RxODE('
         d/dt(X) = a*X + Y*Z;
         d/dt(Y) = b*(Y - Z);
         d/dt(Z) = -X*Y + c*Y - Z;')

et <- eventTable()   # default time units
et$add.sampling(seq(from=0, to=100, by=0.01))

out <- 
    ode$solve(
        params = c(a=-8/3, b=-10, c=28), 
        events = et,
        inits = c(X=1, Y=1, Z=1)
    )

out <- round(out,6)


test_that("Runs example 3.2 correctly",{
    expect_equal(digest(out,"sha512"),
                 "3365655136f0563f22565dd9ead11e62a817f5eaf49d900311bf36cbf8e98eef684baf5173bfd5af740ec3b10e8970e061937a883c020fdfd3f80c846a7b617f");
})


print(head(out, n = 15), digits = 6)

## old <- par(mfrow=c(2,2))

## plot(out[,"time"], out[,"X"], xlab="time", ylab="X", type="l")
## plot(out[,"time"], out[,"Y"], xlab="time", ylab="Y", type="l")
## plot(out[,"time"], out[,"Z"], xlab="time", ylab="Z", type="l")
## plot(out[,"X"], out[,"Y"], xlab="X", ylab="Y", type="l")

## par(old)


rxClean()
