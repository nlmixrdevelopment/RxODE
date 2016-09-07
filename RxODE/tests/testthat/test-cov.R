library(RxODE)

ode <- 
    RxODE('d/dt(X) = a*X + Y*Z;
         d/dt(Y) = b*(Y - Z);
           d/dt(Z) = -X*Y + c*Y - Z')

et <- eventTable()   # default time units
et$add.sampling(seq(from=0, to=100, by=0.01))

cov <- data.frame(c=et$get.sampling()$time/100*24+0.1);

out <- 
    rxSolve(ode,
            params = c(a=-8/3, b=-10), 
            events = et,
            inits = c(X=1, Y=1, Z=1),
            covs = cov);

test_that("time varying covariates output covariate in data frame",{
    expect_equal(cov$c,out$c);
})

rxClean()
