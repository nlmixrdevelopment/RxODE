# (Regression) test 3 multiple instances of RxODE objects to ensure
# C symbols and operations don't conflict. 

library("RxODE")

test.dir <- tempfile("Rxmult-")

# RxODE instance 1 
m1 <- 
   RxODE(
      model = '
         C2 = centr/V2;
         C3 = peri/V3;
         d/dt(depot) =-KA*depot;
         d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
         d/dt(peri)  =                    Q*C2 - Q*C3;
         d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;',
      modName = "inst1",
      wd = test.dir
   )
et1 <- eventTable(amount.units="ug", time.units = "hours")
et1$add.dosing(dose=10000, nbr.doses=5, dosing.interval = 24)
et1$add.sampling(0:24)
et1$add.sampling(seq(from = 24+8, to = 5*24, by = 8))

o1.first <- 
   m1$solve(
      params = c(KA=.291, CL=18.6, V2=40.2, Q=10.5, V3=297.0,
         Kin=1.0, Kout=1.0, EC50=200.0),
      events = et1, 
      inits = c(0, 0, 0, 1)      
   )

# RxODE instance 2 (complete example)
m2 <- 
   RxODE(
      model = 'd/dt(y) = r * y * (1.0 - y/K);', 
      modName = "inst2",
      wd = test.dir
   )
et2 <- eventTable(time.units= NA)
et2$add.sampling(seq(from = 0, to = 20, by = 0.2))
o2.s <- 
   m2$solve(params=c(r=1, K=10), events=et2, inits=c(y=2), stiff = TRUE)

# RxODE instance 3 (complete example)
m3 <- 
   RxODE(
      model = '
         d/dt(X) = a*X + Y*Z;
         d/dt(Y) = b*(Y - Z);
         d/dt(Z) = -X*Y + c*Y - Z;',
      modName = "inst3",
      wd = tempdir()                   # don't pollute "./tests"
   )
et3 <- eventTable()   # default time units
et3$add.sampling(seq(from=0, to=100, by=0.01))
o3 <- 
   m3$solve( params = c(a=-8/3, b=-10, c=28), 
      events = et3,
      inits = c(X=1, Y=1, Z=1)
   )

# Now go back to model 1 for (same) integration
o1.second <- 
   m1$solve(
      params = c(KA=.291, CL=18.6, V2=40.2, Q=10.5, V3=297.0,
         Kin=1.0, Kout=1.0, EC50=200.0),
      events = et1, 
      inits = c(0, 0, 0, 1)      
   )

all.equal(o1.first, o1.second)

# and go back to model 2 for different ode solver
o2.ns <-
   m2$solve(params=c(r=1, K=10), events=et2, inits=c(y=12), stiff = FALSE)

# Inspect the internal compilation manager in each of the RxODE objects
prt <- 
function(obj)
{
   cat(
      sprintf('Model name: %s\nDyn lib: %s\node_solver symbol: %s\n', 
         obj$modName, obj$cmpMgr$dllfile, obj$cmpMgr$ode_solver
      )
   )
}

prt(m1)
prt(m2)
prt(m3)

# Unload object code (this is platform-dependent, as documented in the 
# "Note" section of help("dyn.load"). Remove test.dir.
m1$dynUnload()
m2$dynUnload()
m3$dynUnload()
print(unlink(test.dir, recursive = TRUE)) # 0==success, 1==failed


