## ---- echo=FALSE--------------------------------------------------------------
options(cli.unicode=FALSE, crayon.enabled=FALSE);
options(knitr.table.format = "html")
knitr::opts_chunk$set( comment = "#>")
htmltools::img(src = knitr::image_uri("logo.png"), 
               alt = 'RxODE', 
               style = 'position:absolute; top:0; right:0; padding:10px; border: 0;')
options(width=80)
Sys.setenv(RSTUDIO_CONSOLE_WIDTH=80)

## -----------------------------------------------------------------------------
library(RxODE)
(ev <- eventTable());

## -----------------------------------------------------------------------------
(ev <- et());

## -----------------------------------------------------------------------------
## Model from RxODE tutorial
m1 <-RxODE({
    KA=2.94E-01;
    CL=1.86E+01;
    V2=4.02E+01;
    Q=1.05E+01;
    V3=2.97E+02;
    Kin=1;
    Kout=1;
    EC50=200;
    ## Added modeled bioavaiblity, duration and rate
    fdepot = 1;
    durDepot = 8;
    rateDepot = 1250;
    C2 = centr/V2;
    C3 = peri/V3;
    d/dt(depot) =-KA*depot;
    f(depot) = fdepot
    dur(depot) = durDepot
    rate(depot) = rateDepot
    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
    d/dt(peri)  =                    Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
    eff(0) = 1
});

## -----------------------------------------------------------------------------
ev <- eventTable(amount.units="mg", time.units="hr")

## The methods ar attached to the event table, so you can use them
## directly
ev$add.dosing(dose=10000, nbr.doses = 3)# loading doses
## Starts at time 0; Default dosing interval is 24

## You can also pipe the event tables to these methods.
ev <- ev %>%
    add.dosing(dose=5000, nbr.doses=14, dosing.interval=12)# maintenance

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(C2)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, until = set_units(3, days), ii=12) # loading doses

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(C2)

## -----------------------------------------------------------------------------
ev <- eventTable(amount.units="mg", time.units="hr")

## The methods ar attached to the event table, so you can use them
## directly
ev$add.dosing(dose=10000, nbr.doses = 3)# loading doses

ev$add.sampling(seq(0,24,by=4))

ev

## -----------------------------------------------------------------------------
solve(m1, ev) %>% plot(C2)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, until = set_units(3, days), ii=12) %>% # loading doses
    et(seq(0,24,by=4))

ev

## -----------------------------------------------------------------------------
solve(m1, ev) %>% plot(C2)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, until = set_units(3, days), ii=12) %>% # loading doses
    et(seq(0,48,length.out=200)) %>%
    et(id=1:4)

ev

## -----------------------------------------------------------------------------
set.seed(42)
solve(m1, ev,
      params=data.frame(KA=0.294*exp(rnorm(4)), 18.6*exp(rnorm(4)))) %>%
    plot(C2)

## -----------------------------------------------------------------------------
set.seed(42)
ev <- et(timeUnits="hr") %>%
    et(time=c(0,6), amt=10000, until = set_units(2, days), ii=12) %>% # loading doses
    et(id=1:4)

ev

## -----------------------------------------------------------------------------
ev <- ev %>% et(seq(0,48,length.out=200))

solve(m1, ev, params=data.frame(KA=0.294*exp(rnorm(4)), 18.6*exp(rnorm(4)))) %>% plot(C2)

## -----------------------------------------------------------------------------
set.seed(42)
ev <- et(timeUnits="hr") %>%
    et(time=c(0,2), amt=10000, until = set_units(2, days), ii=12) %>% # loading doses
    et(id=1:4) %>%
    et(seq(0,48,length.out=200))

solve(m1, ev, params=data.frame(KA=0.294*exp(rnorm(4)), 18.6*exp(rnorm(4)))) %>% plot(C2)

## -----------------------------------------------------------------------------
set.seed(42)
ev <- et(timeUnits="hr") %>%
    et(time=c(0,2), amt=10000, until = set_units(2, days), ii=12) %>% # loading doses
    et(id=1:4)

## Create 20 samples in the first 24 hours and 20 samples in the second 24 hours
samples <- c(lapply(1:20, function(...){c(0,24)}),
             lapply(1:20, function(...){c(20,48)}))

## Add the random collection to the event table
ev <- ev %>% et(samples)

library(ggplot2)
solve(m1, ev, params=data.frame(KA=0.294*exp(rnorm(4)), 18.6*exp(rnorm(4)))) %>% plot(C2) + geom_point()

## -----------------------------------------------------------------------------
## bid for 5 days
bid <- et(timeUnits="hr") %>%
       et(amt=10000,ii=12,until=set_units(5, "days"))

## qd for 5 days
qd <- et(timeUnits="hr") %>%
      et(amt=20000,ii=24,until=set_units(5, "days"))

## bid for 5 days followed by qd for 5 days
et <- seq(bid,qd) %>% et(seq(0,11*24,length.out=100));

rxSolve(m1, et) %>% plot(C2)

## -----------------------------------------------------------------------------
## bid for 5 days followed by qd for 5 days
et <- seq(bid,set_units(1, "week"), qd) %>%
    et(seq(0,18*24,length.out=100));

rxSolve(m1, et) %>% plot(C2)

## -----------------------------------------------------------------------------
## bid for 5 days followed by qd for 5 days
et <- seq(bid,set_units(1, "week"), qd,wait="+ii") %>%
    et(seq(0,18*24,length.out=100));

rxSolve(m1, et) %>% plot(C2)

## -----------------------------------------------------------------------------
qd <-et(timeUnits = "hr") %>% et(amt=10000, ii=24, until=set_units(2, "weeks"), cmt="depot")

et <- rep(qd, times=4, wait=set_units(1,"weeks")) %>%
      add.sampling(set_units(seq(0, 12.5,by=0.005),weeks))

rxSolve(m1, et)  %>% plot(C2)

## -----------------------------------------------------------------------------
## bid for 5 days
bid <- et(timeUnits="hr") %>%
       et(amt=10000,ii=12,until=set_units(5, "days"))

## qd for 5 days
qd <- et(timeUnits="hr") %>%
      et(amt=20000,ii=24,until=set_units(5, "days"))

et <- seq(bid,qd) %>%
    et(seq(0,18*24,length.out=500));

rxSolve(m1, et) %>% plot(C2)

## -----------------------------------------------------------------------------
## bid for 5 days
et <- rbind(bid,qd) %>%
    et(seq(0,18*24,length.out=500));

rxSolve(m1, et) %>% plot(C2)

## -----------------------------------------------------------------------------
et <- rbind(bid,wait=set_units(10,days),qd) %>%
    et(seq(0,18*24,length.out=500));

rxSolve(m1, et) %>% plot(C2)

## -----------------------------------------------------------------------------
## bid for 5 days
et <- etRbind(bid,qd, id="unique") %>%
    et(seq(0,150,length.out=500));

library(ggplot2)
rxSolve(m1, et) %>% plot(C2) + facet_wrap( ~ id)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12,until=24) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(C2) +
    xlab("Time")

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12,until=24, dur=8) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(depot, C2) +
    xlab("Time")

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12,until=24, rate=10000/8) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(depot, C2) +
    xlab("Time")

## -----------------------------------------------------------------------------
rxSolve(m1, ev, c(fdepot=0.25)) %>% plot(depot, C2) +
    xlab("Time")

## -----------------------------------------------------------------------------
rxSolve(m1, ev, c(fdepot=1.25)) %>% plot(depot, C2) +
    xlab("Time")

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12,until=24, dur=8) %>%
    et(seq(0, 24, length.out=100))


## -----------------------------------------------------------------------------
library(ggplot2)
library(gridExtra)

p1 <- rxSolve(m1, ev, c(fdepot=1.25)) %>% plot(depot) +
    xlab("Time") + ylim(0,5000)

p2 <- rxSolve(m1, ev, c(fdepot=0.25)) %>% plot(depot) +
    xlab("Time")+ ylim(0,5000)

grid.arrange(p1,p2, nrow=1)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12,until=24, dur=model) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev, c(durDepot=7)) %>% plot(depot, C2) +
    xlab("Time")


## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12,until=24, rate=model) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev, c(rateDepot=10000/3)) %>% plot(depot, C2) +
    xlab("Time")

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12, ss=1) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(C2)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=24, ss=1) %>%
    et(time=12, amt=15000, ii=24, ss=2) %>%
    et(time=24, amt=10000, ii=24, addl=3) %>%
    et(time=36, amt=15000, ii=24, addl=3) %>%
    et(seq(0, 64, length.out=500))

library(ggplot2)

rxSolve(m1, ev,maxsteps=10000) %>% plot(C2) +
    annotate("rect", xmin=0, xmax=24, ymin=-Inf, ymax=Inf, alpha=0.2) +
    annotate("text", x=12.5, y=7, label="Initial Steady State Period") +
    annotate("text", x=44,   y=7, label="Steady State AM/PM dosing")


## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12, addl=3) %>%
    et(time=6, evid=reset) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(depot,C2, eff)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12, addl=3) %>%
    et(time=6, amt=10000, evid=4) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(depot,C2, eff)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12, addl=3) %>%
    et(time=6, cmt="-depot", evid=2) %>%
    et(seq(0, 24, length.out=100))

ev

## -----------------------------------------------------------------------------
rxSolve(m1, ev) %>% plot(depot,C2, eff)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12, addl=3) %>%
    et(time=6, cmt="-eff", evid=2) %>%
    et(seq(0, 24, length.out=100))

rxSolve(m1, ev) %>% plot(depot,C2, eff)

## -----------------------------------------------------------------------------
ev <- et(timeUnits="hr") %>%
    et(amt=10000, ii=12, addl=3) %>%
    et(time=6, cmt="-eff", evid=2) %>%
    et(time=12,cmt="eff",evid=2) %>%
    et(seq(0, 24, length.out=100))

rxSolve(m1, ev) %>% plot(depot,C2, eff)


