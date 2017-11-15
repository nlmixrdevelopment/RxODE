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
inits <- c(0, 0, 0, 1)

## ------------------------------------------------------------------------
inits <- c(eff=1);

## ------------------------------------------------------------------------
ev <- eventTable(amount.units='mg', time.units='hours')

## ------------------------------------------------------------------------
ev$add.dosing(dose=10000, nbr.doses=10, dosing.interval=12)
ev$add.dosing(dose=20000, nbr.doses=5, start.time=120, dosing.interval=24)
ev$add.sampling(0:240)

## ------------------------------------------------------------------------
ev <- eventTable(amount.units="mg", time.units="hours") %>%
    add.dosing(dose=10000, nbr.doses=10, dosing.interval=12) %>%
    add.dosing(dose=20000, nbr.doses=5, start.time=120,dosing.interval=24) %>%
    add.sampling(0:240);

## ------------------------------------------------------------------------
head(ev$get.dosing())

## ------------------------------------------------------------------------
head(ev$get.sampling())

## ------------------------------------------------------------------------
x <- mod1$solve(theta, ev, inits)
head(x)

## ------------------------------------------------------------------------
par(mfrow=c(1,2))
matplot(x[,"C2"], type="l", ylab="Central Concentration")
matplot(x[,"eff"], type="l", ylab = "Effect")

## ------------------------------------------------------------------------
x <- predict(mod1,theta, ev, inits)
print(x)

## ------------------------------------------------------------------------
x <- solve(mod1,theta, ev, inits)
print(x)

## ------------------------------------------------------------------------
x <- mod1 %>% solve(theta, ev, inits)
print(x)

## ------------------------------------------------------------------------
library(dplyr)
x <- mod1 %>% solve(theta,ev,inits) %>%  filter(time <=3)
x

## ------------------------------------------------------------------------
x <- mod1 %>% solve(theta,ev,inits);

## ------------------------------------------------------------------------
x$eff0

## ------------------------------------------------------------------------
x$eff0 <- 2
x

## ------------------------------------------------------------------------
x$t <- seq(0,5,length.out=20)
x

## ------------------------------------------------------------------------
x$KA

## ------------------------------------------------------------------------
x$KA <- 1;
x

## ------------------------------------------------------------------------
mod2 <- RxODE({
    ## the order of variables do not matter, the type of compartmental
    ## model is determined by the parameters specified.
    C2   = linCmt(KA, CL, V2, Q, V3);
    eff(0) = 1  ## This specifies that the effect compartment starts at 1.
    d/dt(eff) =  Kin - Kout*(1-C2/(EC50+C2))*eff;
})

## ------------------------------------------------------------------------
ev <- eventTable(amount.units='mg', time.units='hours') %>%
    add.dosing(dose=10000, nbr.doses=10, dosing.interval=12,dosing.to=2) %>%
    add.dosing(dose=20000, nbr.doses=5, start.time=120,dosing.interval=24,dosing.to=2) %>%
    add.sampling(0:240);

## ------------------------------------------------------------------------
(x <- mod2 %>%  solve(theta, ev))

## ------------------------------------------------------------------------
(x <- mod2 %>%  solve(theta, ev,c(eff=2)))

## ------------------------------------------------------------------------
mod3 <- RxODE({
    KA=2.94E-01;
    CL=1.86E+01; 
    V2=4.02E+01; 
    Q=1.05E+01;
    V3=2.97E+02; 
    Kin=1;
    Kout=1;
    EC50=200;
    ## The linCmt() picks up the variables from above
    C2   = linCmt();
    eff(0) = 1  ## This specifies that the effect compartment starts at 1.
    d/dt(eff) =  Kin - Kout*(1-C2/(EC50+C2))*eff;
})

(x <- mod3 %>%  solve(ev))

## ------------------------------------------------------------------------
(x <- mod3 %>%  solve(c(KA=10),ev))

## ------------------------------------------------------------------------
mod3 <- RxODE({
    KA=2.94E-01;
    CL=1.86E+01; 
    V2=4.02E+01; 
    Q=1.05E+01;
    V3=2.97E+02; 
    Kin0=1;
    Kout=1;
    EC50=200;
    ## The linCmt() picks up the variables from above
    C2   = linCmt();
    Tz= 8
    amp=0.1
    eff(0) = 1  ## This specifies that the effect compartment starts at 1.
    ## Kin changes based on time of day (like cortosol)
    Kin =   Kin0 +amp *cos(2*pi*(ctime-Tz)/24)
    d/dt(eff) =  Kin - Kout*(1-C2/(EC50+C2))*eff;
})


ev <- eventTable(amount.units="mg", time.units="hours") %>%
    add.dosing(dose=10000, nbr.doses=1, dosing.to=2) %>%
    add.sampling(seq(0,48,length.out=100));


 ## Create data frame of  8 am dosing for the first dose
cov.df  <- data.frame(ctime =(seq(0,48,length.out=100)+8) %% 24);


## ------------------------------------------------------------------------

(r1 <- solve(mod3, ev, covs=cov.df,covs_interpolation="linear"))

## ----out.width="100%"----------------------------------------------------
par(mfrow=c(1,2))
matplot(r1[,"C2"], type="l", ylab="Central Concentration")
matplot(r1[,"eff"], type="l", ylab = "Effect")

## ------------------------------------------------------------------------
(r2 <- solve(mod3, ev, covs=cov.df,covs_interpolation="constant"))

## ----out.width="100%"----------------------------------------------------
par(mfrow=c(1,2))
matplot(r2[,"C2"], type="l", ylab="Central Concentration")
matplot(r2[,"eff"], type="l", ylab = "Effect")

## ----out.width="100%"----------------------------------------------------

mod <- RxODE({
    ## Table 3 from Savic 2007
    cl = 17.2 # (L/hr)
    vc = 45.1 # L
    ka = 0.38 # 1/hr
    mtt = 0.37 # hr
    bio=1
    n = 20.1
    k = cl/vc
    ktr = (n+1)/mtt
    ## note that lgammafn is the same as lgamma in R.
    d/dt(depot) = exp(log(bio*podo)+log(ktr)+n*log(ktr*t)-ktr*t-lgammafn(n+1))-ka*depot
    d/dt(cen) = ka*depot-k*cen
})

et <- eventTable();
et$add.sampling(seq(0, 7, length.out=200));
et$add.dosing(20, start.time=0);

transit <- rxSolve(mod, et, transit_abs=TRUE)

par(mfrow=c(1,1))
with(transit,matplot(time,cen, type="l", ylab="Central Concentration", xlab=""))


## ----out.width="100%"----------------------------------------------------

mod <- RxODE({
    ## Table 3 from Savic 2007
    cl = 17.2 # (L/hr)
    vc = 45.1 # L
    ka = 0.38 # 1/hr
    mtt = 0.37 # hr
    bio=1
    n = 20.1
    k = cl/vc
    ktr = (n+1)/mtt
    d/dt(depot) = transit(n,mtt,bio)-ka*depot
    d/dt(cen) = ka*depot-k*cen
})

et <- eventTable();
et$add.sampling(seq(0, 7, length.out=200));
et$add.dosing(20, start.time=0);

transit <- rxSolve(mod, et)

par(mfrow=c(1,1))
with(transit,matplot(time,cen, type="l", ylab="Central Concentration", xlab=""))

## ------------------------------------------------------------------------
Vtpol2 <- RxODE({
    d/dt(y)  = dy
    d/dt(dy) = mu*(1-y^2)*dy - y
    ## Jacobian
    df(y)/dy(dy)  = 1
    df(dy)/dy(y)  = -2*dy*mu*y - 1
    df(dy)/dy(dy) = mu*(1-y^2)
    ## Initial conditions
    y(0) = 2
    dy(0) = 0
    ## mu
    mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
})

et <- eventTable();
et$add.sampling(seq(0, 10, length.out=200));
et$add.dosing(20, start.time=0);

(s1 <- Vtpol2 %>%  solve(et))

## ------------------------------------------------------------------------
(s2 <- Vtpol2 %>%  solve(c(mu=1000), et))

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
set.seed(1)
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

