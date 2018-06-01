## ---- echo=FALSE---------------------------------------------------------
options(cli.unicode=FALSE, crayon.enabled=FALSE);
options(knitr.table.format = "html")
htmltools::img(src = knitr::image_uri("logo.png"), 
               alt = 'RxODE', 
               style = 'position:absolute; top:0; right:0; padding:10px; border: 0;')

## ------------------------------------------------------------------------
library(RxODE)
mod1 <-RxODE({
    C2 = centr/V2;
    C3 = peri/V3;
    d/dt(depot) =-KA*depot;
    d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
    d/dt(peri)  =                    Q*C2 - Q*C3;
    d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;
});

## ------------------------------------------------------------------------
theta <- 
   c(KA=2.94E-01, CL=1.86E+01, V2=4.02E+01, # central 
     Q=1.05E+01,  V3=2.97E+02,              # peripheral
     Kin=1, Kout=1, EC50=200)               # effects

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
knitr::kable(head(ev$get.dosing()))

## ------------------------------------------------------------------------
knitr::kable(head(ev$get.sampling()))

## ---- results="asis"-----------------------------------------------------
x <- solve(mod1,theta, ev, inits);
rxHtml(x)

## ------------------------------------------------------------------------
x <- mod1$solve(theta, ev, inits)
knitr::kable(head(x))

## ------------------------------------------------------------------------
library(ggplot2)
x <- as.data.frame(x)
ggplot(x,aes(time,C2)) + geom_line() + ylab("Central Concentration") + xlab("Time");

## ------------------------------------------------------------------------
ggplot(x,aes(time,eff)) + geom_line() + ylab("Effect") + xlab("Time");

