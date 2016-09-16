## Example 3.2 from 
## "Solving Differential Equations in R" by Soetaert et al (2012)

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
                 "d9f5a3d1b1bcf32140f823bb42f1910a2d181db2f5d3ee6eb73e42a4fd97150027574f9ef2742c0da7d8ce19f4ced3bca680324d38983fb25d7b04da3700a8e1");
})


print(head(out, n = 15), digits = 6)

## old <- par(mfrow=c(2,2))
## 
## plot(out[,"time"], out[,"X"], xlab="time", ylab="X", type="l")
## plot(out[,"time"], out[,"Y"], xlab="time", ylab="Y", type="l")
## plot(out[,"time"], out[,"Z"], xlab="time", ylab="Z", type="l")
## plot(out[,"X"], out[,"Y"], xlab="X", ylab="Y", type="l")
## 
## par(old)


rxClean()
