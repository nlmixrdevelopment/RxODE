## ---- echo=FALSE---------------------------------------------------------
options(cli.unicode=FALSE, crayon.enabled=FALSE);
options(knitr.table.format = "html")
htmltools::img(src = knitr::image_uri("logo.png"), 
               alt = 'RxODE', 
               style = 'position:absolute; top:0; right:0; padding:10px; border: 0;')

## ----results="asis"------------------------------------------------------

library(RxODE)
library(units)

## Setup example model
mod1 <-RxODE({
    C2 = centr/V2;
    C3 = peri/V3;
    d/dt(depot) =-KA*depot;
    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
    d/dt(peri)  =                    Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
});

## Seup parameters and initial conditions

theta <- 
    c(KA=2.94E-01, CL=1.86E+01, V2=4.02E+01, # central 
      Q=1.05E+01,  V3=2.97E+02,              # peripheral
      Kin=1, Kout=1, EC50=200)               # effects

inits <- c(eff=1);

## Setup dosing event information
ev <- eventTable(amount.units="mg", time.units="hours") %>%
    add.dosing(dose=10000, nbr.doses=10, dosing.interval=12) %>%
    add.dosing(dose=20000, nbr.doses=5, start.time=120,dosing.interval=24) %>%
    add.sampling(0:240);


## Now solve
x <- predict(mod1,theta, ev, inits)
rxHtml(x)

## ----results="asis"------------------------------------------------------
x <- solve(mod1,theta, ev, inits)
rxHtml(x)

## ----results="asis"------------------------------------------------------
x <- mod1 %>% solve(theta, ev, inits)
rxHtml(x)

## ------------------------------------------------------------------------
library(dplyr)
## You can  drop units for comparisons and filtering
x <- mod1 %>% solve(theta,ev,inits) %>% drop_units %>% filter(time <= 3) %>% as.tbl
## or keep them and compare with the proper units.
x <- mod1 %>% solve(theta,ev,inits) %>% filter(time <= set_units(3, hr)) %>% as.tbl
x

## ------------------------------------------------------------------------
x <- mod1 %>% solve(theta,ev,inits);

## ------------------------------------------------------------------------
x$eff0

## ---- results="asis"-----------------------------------------------------
x$eff0 <- 2
rxHtml(x)

## ----results="asis"------------------------------------------------------
x$t <- seq(0,5,length.out=20)
rxHtml(x)

## ------------------------------------------------------------------------
x$KA

## ----results="asis"------------------------------------------------------
x$KA <- 1;
rxHtml(x)

