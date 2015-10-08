# Example 3.2 from 
# "Solving Differential Equations in R" by Soetaert et al (2012)

library("RxODE")
        
ode <- 
   RxODE(
      model = '
         d/dt(X) = a*X + Y*Z;
         d/dt(Y) = b*(Y - Z);
         d/dt(Z) = -X*Y + c*Y - Z;',
      modName = "test_3-2",
      wd = tempdir()                   # don't pollute "./tests"
   )

et <- eventTable()   # default time units
et$add.sampling(seq(from=0, to=100, by=0.01))

out <- 
   ode$solve(
      params = c(a=-8/3, b=-10, c=28), 
      events = et,
      inits = c(X=1, Y=1, Z=1)
   )
   
print(head(out, n = 15), digits = 6)


# old <- par(mfrow=c(2,2))
# 
# plot(out[,"time"], out[,"X"], xlab="time", ylab="X", type="l")
# plot(out[,"time"], out[,"Y"], xlab="time", ylab="Y", type="l")
# plot(out[,"time"], out[,"Z"], xlab="time", ylab="Z", type="l")
# plot(out[,"X"], out[,"Y"], xlab="X", ylab="Y", type="l")
# 
# par(old)


