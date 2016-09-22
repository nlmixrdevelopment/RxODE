## Example 3.2 from 
## "Solving Differential Equations in R" by Soetaert et al (2012)
## https://cran.r-project.org/web/packages/diffEq/vignettes/ODEinR.pdf Example #1
library(digest)
old <- options(digits = 6)

ode <- 
    RxODE(model = 'd/dt(y) = r * y * (1.0 - y/K);')

## create event table with times at which we observe the system
et <- eventTable(time.units= NA)
et$add.sampling(seq(from = 0, to = 20, by = 0.2))

## same model, different initial values (2 and 12)
out1 <- ode$solve(params=c(r=1, K=10), events=et, inits=c(y=2))
out2 <- ode$solve(params=c(r=1, K=10), events=et, inits=c(y=12))

##matplot(x = out1[,1], y = cbind(out1[,2], out2[,2]), type = "l",
##    main = "logistic growth", xlab="time", ylab="")

## Now use a non-stiff solver

out1.ns <- 
    ode$solve(params=c(r=1, K=10), events=et, inits=c(y=2), stiff=FALSE)
out2.ns <- 
    ode$solve(params=c(r=1, K=10), events=et, inits=c(y=12), stiff=FALSE)

df <-round(cbind(out1, out2, out1.ns, out2.ns),6);

test_that("Runs example 3.1 correctly",{
    expect_equal(digest(df,"sha512"),
                 "968de0db2de21958f278749e77987a9a5e976dd1c3eafa51d8cb40dd128e122d9d6dc80bd93e39f44f5339503c1390cbb65ae75b3fbf72d4cbfa0a3b49545cc5");
})

head(cbind(out1, out2, out1.ns, out2.ns), n = 15)

options(old)

rxClean()
