## ---- echo=FALSE---------------------------------------------------------
options(cli.unicode=FALSE, crayon.enabled=FALSE);
options(knitr.table.format = "html")
htmltools::img(src = knitr::image_uri("logo.png"), 
               alt = 'RxODE', 
               style = 'position:absolute; top:0; right:0; padding:10px; border: 0;')

## ------------------------------------------------------------------------
library(RxODE)
library(units)

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
ev$ctime <- (ev$time+set_units(8,hr)) %% 24


## ----results="asis"------------------------------------------------------

r1 <- solve(mod3, ev, covs_interpolation="linear")
rxHtml(r1)

## ------------------------------------------------------------------------
library(ggplot2)
plot(r1,C2) +ylab("Central Concentration") + xlab("Time");

## ------------------------------------------------------------------------
plot(r1,eff) + ylab("Effect") + xlab("Time");

## ----results="asis"------------------------------------------------------
r1 <- solve(mod3, ev,covs_interpolation="constant")
rxHtml(r1)

## ------------------------------------------------------------------------
plot(r1,C2) + ylab("Central Concentration") + xlab("Time");

## ------------------------------------------------------------------------
plot(r1,eff) + ylab("Effect") + xlab("Time");

## ----results="asis"------------------------------------------------------
r1 <- solve(mod3, ev,covs_interpolation="nocb")
rxHtml(r1)

## ------------------------------------------------------------------------
plot(r1,C2) + ylab("Central Concentration") + xlab("Time");

## ------------------------------------------------------------------------
plot(r1,eff) + ylab("Effect") + xlab("Time");

