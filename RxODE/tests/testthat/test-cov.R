library(RxODE)

ode <- 
    RxODE('d/dt(X) = a*X + Y*Z;
         d/dt(Y) = b*(Y - Z);
         d/dt(Z) = -X*Y + c*Y - Z
if (t < 0.02 | t > 99.98){
print
}
')

et <- eventTable()   # default time units
et$add.sampling(seq(from=0, to=100, by=0.01))

cov <- data.frame(c=et$get.sampling()$time+1);

out <- rxSolve(ode,
               params = c(a=-8/3, b=-10), 
               events = et,
               inits = c(X=1, Y=1, Z=1),
               covs = cov);

test_that("time varying covariates output covariate in data frame",{
    expect_equal(cov$c,out$c);
})

out <- as.data.frame(out);
out <- out[,names(out) != "c"];

out1 <- 
    rxSolve(ode,
            params = c(a=-8/3, b=-10, c = 0), 
            events = et,
            inits = c(X=1, Y=1, Z=1))

out1 <- as.data.frame(out1);

test_that("time varying covariates produce different outputs",{
    expect_false(isTRUE(all.equal(out,out1)));
})


cov <- data.frame(c=et$get.sampling()$time+1,a=-et$get.sampling()$time/100);

out <- rxSolve(ode,
               params = c(a=-8/3, b=-10), 
               events = et,
               inits = c(X=1, Y=1, Z=1),
               covs = cov)

test_that("time varying covariates output covariate in data frame",{
    expect_equal(cov$c,out$c);
    expect_equal(cov$a,out$a);
})


rxClean()
