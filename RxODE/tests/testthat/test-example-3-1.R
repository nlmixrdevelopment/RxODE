library(digest)
old <- options(digits = 6)

ode <- 
    RxODE(
        model = 'd/dt(y) = r * y * (1.0 - y/K);', 
        modName = "test_3-1")

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
                 "9bd5ee859f6ef217cd5e047b4010a5018f1d15202db807edb351d7a5e3be7af46de5db7fb843ced3f10a8980443719eacbbb2b6268a3713f287f138b08f35d5b");
})

head(cbind(out1, out2, out1.ns, out2.ns), n = 15)

options(old)

rxClean()
