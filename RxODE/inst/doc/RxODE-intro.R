## ----setup, include = FALSE----------------------------------------------
library(knitr)

## ------------------------------------------------------------------------
ode <- "
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
"

## ------------------------------------------------------------------------
library(RxODE)
work <- tempfile("Rx_intro-")
mod1 <- RxODE(model = ode, modName = "mod1", wd = work)

## ------------------------------------------------------------------------
theta <- 
   c(KA=2.94E-01, CL=1.86E+01, V2=4.02E+01, # central 
     Q=1.05E+01,  V3=2.97E+02,              # peripheral
     Kin=1, Kout=1, EC50=200)               # effects  

## ------------------------------------------------------------------------
inits <- c(depot=0, centr=0, peri=0, eff=1)    

## ------------------------------------------------------------------------
ev <- eventTable(amount.units='mg', time.units='hours')

## ------------------------------------------------------------------------
ev$add.dosing(dose=10000, nbr.doses=10, dosing.interval=12)
ev$add.dosing(dose=20000, nbr.doses=5, start.time=120, dosing.interval=24)
ev$add.sampling(0:240)

## ------------------------------------------------------------------------
head(ev$get.dosing())

## ------------------------------------------------------------------------
head(ev$get.sampling())

## ------------------------------------------------------------------------
x <- mod1$solve(theta, ev, inits)
head(x)

## ----fig.width=10--------------------------------------------------------
par(mfrow=c(1,2))
matplot(x[,"C2"], type="l", ylab="Central Concentration")
matplot(x[,"eff"], type="l", ylab = "Effect")

## ------------------------------------------------------------------------
nsub <- 100						  #number of subproblems
CL <- 1.86E+01*exp(rnorm(nsub,0,.4^2))
theta.all <- 
	cbind(KA=2.94E-01, CL=CL, V2=4.02E+01,  # central 
	Q=1.05E+01, V3=2.97E+02,                # peripheral
	Kin=1, Kout=1, EC50=200)                # effects  
head(theta.all)

## ----fig.width=10--------------------------------------------------------
nobs <- ev$get.nobs()
cp.all <- matrix(NA, nobs, nsub)
for (i in 1:nsub)
{
	theta <- theta.all[i,]
	x <- mod1$solve(theta, ev, inits=inits)
	cp.all[, i] <- x[, "C2"]
}

matplot(cp.all, type="l", ylab="Central Concentration")


## ----fig.width=10--------------------------------------------------------
cp.q <- apply(cp.all, 1, quantile, prob = c(0.05, 0.50, 0.95))
matplot(t(cp.q), type="l", lty=c(2,1,2), col=c(2,1,2), ylab="Central Concentration")


## ----cleanup, include = FALSE--------------------------------------------
mod1$dynUnload()   # otherwise we can't unlink(work)
unlink(work, recursive = TRUE)

