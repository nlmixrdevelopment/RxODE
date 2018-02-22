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

(s1 <- Vtpol2 %>%  solve(et, method="lsoda"))

## ------------------------------------------------------------------------
(s2 <- Vtpol2 %>%  solve(c(mu=1000), et))

## ------------------------------------------------------------------------
set.seed(42);
mod <- RxODE({
    eff(0) = 1
    C2 = centr/V2;
    C3 = peri/V3;
    CL =  TCl*exp(eta.Cl) ## This is coded as a variable in the model
    d/dt(depot) =-KA*depot;
    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
    d/dt(peri)  =                    Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
})

theta <- c(KA=2.94E-01, TCl=1.86E+01, V2=4.02E+01,  # central 
               Q=1.05E+01, V3=2.97E+02,                # peripheral
               Kin=1, Kout=1, EC50=200)                # effects  

## ----fig.width=10--------------------------------------------------------
## the column names of the omega matrix need to match the parameters specified by RxODE
omega <- matrix(0.4^2,dimnames=list(NULL,c("eta.Cl")))

ev <- eventTable(amount.units="mg", time.units="hours") %>%
    add.dosing(dose=10000, nbr.doses=1, dosing.to=2) %>%
    add.sampling(seq(0,48,length.out=100));

sim  <- rxSolve(mod,theta,ev,omega=omega,nSub=100)

library(ggplot2)
library(gridExtra)

p1 <- ggplot(sim,aes(time,centr,color=factor(sim.id))) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)") + guides(color=FALSE)

p2 <-ggplot(sim,aes(time,eff,color=factor(sim.id))) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

grid.arrange(p1,p2,nrow=2)

## ----fig.width=10--------------------------------------------------------

library(dplyr)

p <- c(0.05, 0.5, 0.95);
s <-sim %>% group_by(time) %>%
    do(data.frame(p=p, eff=quantile(.$eff, probs=p), 
                  eff.n = length(.$eff), eff.avg = mean(.$eff),
                  centr=quantile(.$centr, probs=p),
                  centr.n=length(.$centr),centr.avg = mean(.$centr))) %>%
    mutate(Percentile=factor(sprintf("%d%%",p*100),levels=c("5%","50%","95%")))

p1 <- ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

p2 <-ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

grid.arrange(p1,p2,nrow=2)

## ------------------------------------------------------------------------
head(sim$param)

## ------------------------------------------------------------------------
theta <- sim$param;
(sim  <- rxSolve(mod,theta,ev))

## ------------------------------------------------------------------------
sim$eff0 <- 100
sim

## ------------------------------------------------------------------------
mod <- RxODE({
    eff(0) = 1
    C2 = centr/V2;
    C3 = peri/V3;
    CL =  TCl*exp(eta.Cl) ## This is coded as a variable in the model
    d/dt(depot) =-KA*depot;
    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
    d/dt(peri)  =                    Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    e = eff+eff.err
    cp = centr*(1+cp.err)
})

theta <- c(KA=2.94E-01, TCl=1.86E+01, V2=4.02E+01,  # central 
           Q=1.05E+01, V3=2.97E+02,                # peripheral
           Kin=1, Kout=1, EC50=200)                # effects  

sigma <- diag(2)*0.1
dimnames(sigma) <- list(NULL, c("eff.err","cp.err"))


sim  <- rxSolve(mod, theta, ev, omega=omega, nSub=100, sigma=sigma)

p <- c(0.05, 0.5, 0.95);
s <-sim %>% group_by(time) %>%
    do(data.frame(p=p, eff=quantile(.$e, probs=p), 
                  eff.n = length(.$e), eff.avg = mean(.$e),
                  centr=quantile(.$cp, probs=p),
                  centr.n=length(.$cp),centr.avg = mean(.$cp))) %>%
    mutate(Percentile=factor(sprintf("%d%%",p*100),levels=c("5%","50%","95%")))

p1 <- ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

p2 <-ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

grid.arrange(p1,p2,nrow=2)

## ------------------------------------------------------------------------

ev1 <- eventTable(amount.units="mg", time.units="hours") %>%
    add.dosing(dose=10000, nbr.doses=1, dosing.to=2) %>%
    add.sampling(seq(0,48,length.out=10));

ev2 <- eventTable(amount.units="mg", time.units="hours") %>%
    add.dosing(dose=5000, nbr.doses=1, dosing.to=2) %>%
    add.sampling(seq(0,48,length.out=8));

dat <- rbind(data.frame(ID=1, ev1$get.EventTable()),
             data.frame(ID=2, ev2$get.EventTable()))


## Note the number of subject is not needed since it is determined by the data
sim  <- rxSolve(mod, theta, dat, omega=omega, sigma=sigma)

sim %>% select(id, time, e, cp)

## ------------------------------------------------------------------------

## Creating covariance matrix
tmp <- matrix(rnorm(8^2), 8, 8)
tMat <- tcrossprod(tmp, tmp) / (8 ^ 2)
dimnames(tMat) <- list(NULL, names(theta))

sim  <- rxSolve(mod, theta, ev, omega=omega, nSub=100, sigma=sigma, thetaMat=tMat, nStud=10)

p <- c(0.05, 0.5, 0.95);
s <-sim %>% group_by(time) %>%
    do(data.frame(p=p, eff=quantile(.$e, probs=p), 
                  eff.n = length(.$e), eff.avg = mean(.$e),
                  centr=quantile(.$cp, probs=p),
                  centr.n=length(.$cp),centr.avg = mean(.$cp))) %>%
    mutate(Percentile=factor(sprintf("%d%%",p*100),levels=c("5%","50%","95%")))

p1 <- ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

p2 <-ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

grid.arrange(p1,p2,nrow=2)


## ------------------------------------------------------------------------
sim$omega.list

## ------------------------------------------------------------------------
sim$sigma.list

## ------------------------------------------------------------------------
sim  <- rxSolve(mod, theta, ev, omega=omega, nSub=100, sigma=sigma, thetaMat=tMat, nStud=10,
                simVariability=FALSE);

p <- c(0.05, 0.5, 0.95);
s <-sim %>% group_by(time) %>%
    do(data.frame(p=p, eff=quantile(.$e, probs=p), 
                  eff.n = length(.$e), eff.avg = mean(.$e),
                  centr=quantile(.$cp, probs=p),
                  centr.n=length(.$cp),centr.avg = mean(.$cp))) %>%
    mutate(Percentile=factor(sprintf("%d%%",p*100),levels=c("5%","50%","95%")))


p1 <- ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

p2 <-ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

grid.arrange(p1,p2,nrow=2)

## ------------------------------------------------------------------------
sim  <- rxSolve(mod, theta, ev, omega=omega, nSub=100, sigma=sigma, thetaMat=tMat, nStud=10,
                nCoresRV=2);

p <- c(0.05, 0.5, 0.95);
s <-sim %>% group_by(time) %>%
    do(data.frame(p=p, eff=quantile(.$e, probs=p), 
                  eff.n = length(.$e), eff.avg = mean(.$e),
                  centr=quantile(.$cp, probs=p),
                  centr.n=length(.$cp),centr.avg = mean(.$cp))) %>%
    mutate(Percentile=factor(sprintf("%d%%",p*100),levels=c("5%","50%","95%")))


p1 <- ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

p2 <-ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

grid.arrange(p1,p2,nrow=2)

## ----cleanup, include = FALSE--------------------------------------------
mod1$dynUnload()   # otherwise we can't unlink(work)
unlink(work, recursive = TRUE)

