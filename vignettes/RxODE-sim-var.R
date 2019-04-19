## ---- echo=FALSE---------------------------------------------------------
options(cli.unicode=FALSE, crayon.enabled=FALSE);
options(knitr.table.format = "html")
htmltools::img(src = knitr::image_uri("logo.png"), 
               alt = 'RxODE', 
               style = 'position:absolute; top:0; right:0; padding:10px; border: 0;')

## ------------------------------------------------------------------------
set.seed(42);
library(RxODE)
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

## ------------------------------------------------------------------------
## the column names of the omega matrix need to match the parameters specified by RxODE
omega <- matrix(0.4^2,dimnames=list(NULL,c("eta.Cl")))

ev <- eventTable(amount.units="mg", time.units="hours") %>%
    add.dosing(dose=10000, nbr.doses=1, dosing.to=2) %>%
    add.sampling(seq(0,48,length.out=100));

sim  <- rxSolve(mod,theta,ev,omega=omega,nSub=100)

library(ggplot2)

ggplot(sim,aes(time,centr,color=factor(sim.id))) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)") + guides(color=FALSE)

## ------------------------------------------------------------------------
ggplot(sim,aes(time,eff,color=factor(sim.id))) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

## ------------------------------------------------------------------------

library(dplyr)

p <- c(0.05, 0.5, 0.95);
s <-sim %>% group_by(time) %>%
    do(data.frame(p=p, eff=quantile(.$eff, probs=p), 
                  eff.n = length(.$eff), eff.avg = mean(.$eff),
                  centr=quantile(.$centr, probs=p),
                  centr.n=length(.$centr),centr.avg = mean(.$centr))) %>%
    mutate(Percentile=factor(sprintf("%d%%",p*100),levels=c("5%","50%","95%")))

ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

## ------------------------------------------------------------------------
ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)


## ------------------------------------------------------------------------
head(sim$param)

## ----results="asis"------------------------------------------------------
theta <- sim$param;
sim  <- rxSolve(mod,theta,ev)
rxHtml(sim)

## ----results="asis"------------------------------------------------------
sim$eff0 <- 100
rxHtml(sim)

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

ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

## ------------------------------------------------------------------------
ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

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

sim  <- rxSolve(mod, theta, ev, omega=omega, nSub=100, sigma=sigma, thetaMat=tMat, nStud=10,
     dfSub=10, dfObs=100)

p <- c(0.05, 0.5, 0.95);
s <-sim %>% group_by(time) %>%
    do(data.frame(p=p, eff=quantile(.$e, probs=p), 
                  eff.n = length(.$e), eff.avg = mean(.$e),
                  centr=quantile(.$cp, probs=p),
                  centr.n=length(.$cp),centr.avg = mean(.$cp))) %>%
    mutate(Percentile=factor(sprintf("%d%%",p*100),levels=c("5%","50%","95%")))

ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

## ------------------------------------------------------------------------
ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

## ------------------------------------------------------------------------
head(sim$omega.list)

## ------------------------------------------------------------------------
head(sim$sigma.list)

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


ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")


## ------------------------------------------------------------------------
ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

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


ggplot(s,aes(time,centr,color=Percentile)) + geom_line(size=1) + coord_trans(y = "log10") + ylab("Central Concentration") +
    xlab("Time (hr)")

## ------------------------------------------------------------------------
ggplot(s,aes(time,eff,color=Percentile)) + geom_line(size=1) + ylab("Effect") +
    xlab("Time (hr)") + guides(color=FALSE)

