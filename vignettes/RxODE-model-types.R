## ---- echo=FALSE---------------------------------------------------------
options(cli.unicode=FALSE, crayon.enabled=FALSE);
options(knitr.table.format = "html")
htmltools::img(src = knitr::image_uri("logo.png"),
               alt = 'RxODE',
               style = 'position:absolute; top:0; right:0; padding:10px; border: 0;')

## ------------------------------------------------------------------------
library(RxODE)
mod <- RxODE({
    ipre <- 10 * exp(-ke * t);
})
mod

## ------------------------------------------------------------------------
et  <- et(seq(0,24,length.out=50))
cmt1 <- rxSolve(mod,et,params=c(ke=0.5))
cmt1

## ------------------------------------------------------------------------
mod <- RxODE({
    ke <- 0.5
    V <- 1
    ipre <- linCmt();
})
mod

## ------------------------------------------------------------------------
et  <- et(amt=10,time=0,cmt=depot) %>%
    et(seq(0,24,length.out=50))
cmt1 <- rxSolve(mod,et,params=c(ke=0.5))
cmt1

## ------------------------------------------------------------------------
library(RxODE)
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

## Setup a mixed solved/ode system:
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
x <- mod2 %>%  solve(theta, ev)
print(x)

## ------------------------------------------------------------------------
x <- mod2 %>%  solve(theta, ev,c(eff=2))
print(x)

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

x <- mod3 %>%  solve(ev)
print(x)

## ------------------------------------------------------------------------
x <- mod3 %>%  solve(c(KA=10),ev)
print(x)

