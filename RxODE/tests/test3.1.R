# Example 3.1 from the book
# "Solving Differential Equations in R" by Soetaert et al (2012)
# Use both the default stiff solver and the non-stiff solver

library("RxODE")

old <- options(digits = 6)

ode <- 
   RxODE(
      model = 'd/dt(y) = r * y * (1.0 - y/K);', 
      modName = "test_3-1",
      wd = tempdir()                   # don't pollute "./tests"
   )

# create event table with times at which we observe the system
et <- eventTable(time.units= NA)
et$add.sampling(seq(from = 0, to = 20, by = 0.2))

# same model, different initial values (2 and 12)
out1 <- ode$solve(params=c(r=1, K=10), events=et, inits=c(y=2))
out2 <- ode$solve(params=c(r=1, K=10), events=et, inits=c(y=12))

#matplot(x = out1[,1], y = cbind(out1[,2], out2[,2]), type = "l",
#    main = "logistic growth", xlab="time", ylab="")

# Now use a non-stiff solver

out1.ns <- 
   ode$solve(params=c(r=1, K=10), events=et, inits=c(y=2), stiff=FALSE)
out2.ns <- 
   ode$solve(params=c(r=1, K=10), events=et, inits=c(y=12), stiff=FALSE)

head(cbind(out1, out2, out1.ns, out2.ns), n = 15)

options(old)
