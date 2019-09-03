context("Test inductive linearization")

test_that("Matrix exponential alone works", {

    mod <- RxODE("
a = 6
b = 0.6
d/dt(intestine) = -a*intestine
d/dt(blood)     = a*intestine - b*blood
", indLin=TRUE)


    et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1 / 24))
    et$add.dosing(dose=2/24,rate=2,start.time=0,
                  nbr.doses=10,dosing.interval=1)

    pk <- rxSolve(mod,et, method="indLin");

    pk2 <- rxSolve(mod,et, method="liblsoda");

    expect_equal(as.data.frame(pk), as.data.frame(pk2))

    mmModel <- RxODE({
        ka = 0.2
        Vc = 4.7
        Vmax <- 7
        Km = 5.7
        d/dt(depot) = -ka * depot
        Cp = center / Vc
        d/dt(center) = ka * depot - Vmax/(Km + Cp) * Cp
    }, indLin=TRUE)

    et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1/24))
    et$add.dosing(dose=2,start.time=0,
                  nbr.doses=10,dosing.interval=6)

    pk <- rxSolve(mmModel,et, method="indLin");

    pk2 <- rxSolve(mmModel,et, method="liblsoda");

    expect_equal(as.data.frame(pk), as.data.frame(pk2))

})

## iSec <- RxODE({
##     d/dt(Ga) = -ka * Ga
##     d/dt(Gt) = ka * Ga - ka * Gt
##     Gprod = Gss * (Clg + Clgi * Iss)
##     d/dt(Gc) = ka * Gt - Gprod + Q / Vp * Gp - (Clg + Clgi * Ie + Q) / Vg * Gc
##     Gc(0) = Gss * Vg
##     d/dt(Gp) = -Q / Vp * Gp + Q / Vg * Gc
##     d/dt(Ge) = Gc * Kge - Ge * Kge
##     d/dt(I) = (Iss * Cli) * (1 + Sincr * Gt) * (Ge / Gss) ^ IPRG - Cli / Vi * I
##     I(0) = Iss * Vi
##     d/dt(Ie) = kie * I - kie * Ie
## }, indLin=TRUE)

## iVL <- RxODE({
##     d/dt(Tni) = lambda - (1 - INHrti) * gamma * Tni * Vin - dni * Tni
##     Tni0 = da * dv * (alphaL + dL) / (gamma * p * (alphaL + fr * dL))
##     Tni(0) = Tni0
##     d/dt(Ta) = fr * (1 - INHrti) * gamma * Tni * Vin + alphaL * Tl - da * Ta
##     Ta(0) = dv * Vin0/
## })
