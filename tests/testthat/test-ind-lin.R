context("Test inductive linearization")

test_that("Matrix exponential alone works", {

    ## Case 1 ME alone from wikipedia
    mod <- RxODE({
        d/dt(x) = 2 * x - y + z
        d/dt(y) = 3 * y - 1 * z
        d/dt(z) = 2 * x + y + 3 * z
        x(0) = 0.1
        y(0) = 0.1
        z(0) = 0.1
    }, indLin=TRUE)

    m <- rxSolve(mod, et(seq(0, 24, length.out=50)), method="indLin")
    m2 <- rxSolve(mod, et(seq(0, 24, length.out=50)), method="lsoda")

    expect_equal(as.data.frame(m), as.data.frame(m2), tol=1e-5)

    ## Now do without indLin in the RxODE

    mod <- RxODE({
        d/dt(x) = 2 * x - y + z
        d/dt(y) = 3 * y - 1 * z
        d/dt(z) = 2 * x + y + 3 * z
        x(0) = 0.1
        y(0) = 0.1
        z(0) = 0.1
    })

    m <- rxSolve(mod, et(seq(0, 24, length.out=50)), method="indLin")
    m2 <- rxSolve(mod, et(seq(0, 24, length.out=50)), method="lsoda")

    expect_equal(as.data.frame(m), as.data.frame(m2), tol=1e-5)


    ## Case 2 ME alone with inhomogenous systems

    mod <- RxODE({
        d/dt(x) = 2 * x - y + z + exp(-2 * t)
        d/dt(y) = 3 * y - 1 * z
        d/dt(z) = 2 * x + y + 3 * z + exp(-2 * t)
        x(0) = 0.1
        y(0) = 0.1
        z(0) = 0.1
    }, indLin=TRUE)

    m <- rxSolve(mod, et(seq(0, 24, length.out=50)), method="indLin")
    m2 <- rxSolve(mod, et(seq(0, 24, length.out=50)), method="lsoda")

    ## gridExtra::grid.arrange(plot(m), plot(m2))

    ## FIXME?
    ## expect_equal(as.data.frame(m), as.data.frame(m2), tol=1e-5)

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

    expect_equal(as.data.frame(pk), as.data.frame(pk2), tol=1e-5)

    ## plot(microbenchmark::microbenchmark(rxSolve(mod,et, method="indLin",indLinMatExpType=1L),rxSolve(mod,et, method="indLin",indLinMatExpType=2L), rxSolve(mod,et, method="indLin",indLinMatExpType=3L), rxSolve(mod,et, method="lsoda")), log="y")

    et2 <- eventTable(time.units="days")
    et2$add.sampling(seq(0,10,by=1 / 24))
    et2$add.dosing(dose=2,start.time=0,
                   nbr.doses=10,dosing.interval=1)

    pk <- rxSolve(mod,et2, method="indLin");

    pk2 <- rxSolve(mod,et2, method="liblsoda");

    expect_equal(as.data.frame(pk), as.data.frame(pk2), tol=1e-5)

    ## Inductive linearization
    mmModel <- RxODE({
        ka = 1
        Vc = 1
        Vmax <- 0.00734
        Km = 0.3672
        Cp = center / Vc
        d/dt(center) =  - Vmax/(Km + Cp) * Cp
    }, indLin=TRUE)


    mmModel <- RxODE({
        ka = 1
        Vc = 1
        Vmax <- 0.00734
        Km = 0.3672
        Cp = center / Vc
        d/dt(center) =  - Vmax/(Km + Cp) * Cp + exp(-10 * t)
    }, indLin=TRUE)


    ## Inductive + 1x1 matrix
    ## FIXME this should be inductive too...
    mmModel <- RxODE({
        ka = 1
        Vc = 1
        Vmax <- 0.00734
        Km = 0.3672
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - Vmax/(Km + Cp) * Cp
        Cp = center / Vc
    }, indLin=TRUE)

    ## This is inductive
    mmModel <- RxODE({
        ka = 1
        Vc = 1
        Vmax <- 0.00734
        Km = 0.3672
        d/dt(depot) = -ka * depot
        Cp = center / Vc
        d/dt(center) = ka * depot - Vmax/(Km + Cp) * Cp
    }, indLin=TRUE)

    mmModel <- RxODE({
        ka = 1
        Vc = 1
        Vmax <- 0.00734
        Km = 0.3672
        V4 <- 4.3
        Q  <- 1.5
        K12<- Q/Vc
        K21<- Q/Vp
        Cp = center / Vc
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - Vmax/(Km + Cp) * Cp + K21*periph-K12*center
        d/dt(periph) =-K21*periph+K12*center;
    }, indLin=TRUE)


    ## Inductive linearization
    mmModel <- RxODE({
        ka = 1
        Vc = 1
        Vmax <- 0.00734
        Km = 0.3672
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

    ## gridExtra::grid.arrange(plot(pk), plot(pk2))

    expect_equal(as.data.frame(pk), as.data.frame(pk2), tol=7e-5)

    mmModel <- RxODE({
        ka = 1
        Vc = 1
        Vmax <- 0.00734
        Km = 0.3672
        d/dt(depot) = -ka * depot
        Cp = center / Vc
        d/dt(center) = ka * depot - Vmax/(Km + Cp) * Cp + 5 * exp(-0.5 * t)
    }, indLin=TRUE)

    pk <- rxSolve(mmModel,et, method="indLin");
    pk2 <- rxSolve(mmModel,et, method="lsoda");

    ## gridExtra::grid.arrange(plot(pk), plot(pk2))
    ## These are not equal...
    ## expect_equal(as.data.frame(pk), as.data.frame(pk2), tol=7e-5)

    ## plot(microbenchmark::microbenchmark(rxSolve(mmModel,et, method="indLin",indLinMatExpType=1L),rxSolve(mmModel,et, method="indLin",indLinMatExpType=2L), rxSolve(mmModel,et, method="indLin",indLinMatExpType=3L), rxSolve(mmModel,et, method="lsoda")), log="y")

    ## Van der Pol Equation
    ## mu = 1000 stiff
    ## me = 1 non-stiff
    ## van <- RxODE({
    ##     y(0) = 2
    ##     d/dt(y)  = dy
    ##     d/dt(dy) = mu*(1-y^2)*dy - y
    ## }, indLin=TRUE)

    ## et <- eventTable();
    ## et$add.sampling(seq(0, 3000, length.out=200));

    ## s1 <- rxSolve(van, et, c(mu=1000), method="lsoda")
    ## s2 <- rxSolve(van, et, c(mu=1000), method="indLin")
    ## s3 <- rxSolve(van, et, c(mu=1000), method="dop853")



    ## gridExtra::grid.arrange(plot(s1), plot(s2), plot(s3))

    ## expect_equal(as.data.frame(s1), as.data.frame(s2))

    ## s1 <- rxSolve(van, et, c(mu=1), method="lsoda")
    ## s2 <- rxSolve(van, et, c(mu=1), method="indLin")
    ## s3 <- rxSolve(van, et, c(mu=1), method="dop853")

    ## gridExtra::grid.arrange(plot(s1), plot(s2), plot(s3))

    ## expect_equal(as.data.frame(s1), as.data.frame(s2))

    ## microbenchmark::microbenchmark(rxSolve(mmModel,et, method="indLin"),
    ##                                rxSolve(mmModel,et, method="liblsoda"))

})

iSec <- RxODE({
    d/dt(Ga) = -ka * Ga
    d/dt(Gt) = ka * Ga - ka * Gt
    Gprod = Gss * (Clg + Clgi * Iss)
    d/dt(Gc) = ka * Gt - Gprod + Q / Vp * Gp - (Clg + Clgi * Ie + Q) / Vg * Gc
    Gc(0) = Gss * Vg
    d/dt(Gp) = -Q / Vp * Gp + Q / Vg * Gc
    d/dt(Ge) = Gc * Kge - Ge * Kge
    d/dt(I) = (Iss * Cli) * (1 + Sincr * Gt) * (Ge / Gss) ^ IPRG - Cli / Vi * I
    I(0) = Iss * Vi
    d/dt(Ie) = kie * I - kie * Ie
}, indLin=TRUE)

## iVL <- RxODE({
##     d/dt(Tni) = lambda - (1 - INHrti) * gamma * Tni * Vin - dni * Tni
##     Tni0 = da * dv * (alphaL + dL) / (gamma * p * (alphaL + fr * dL))
##     Tni(0) = Tni0
##     d/dt(Ta) = fr * (1 - INHrti) * gamma * Tni * Vin + alphaL * Tl - da * Ta
##     Ta(0) = dv * Vin0/
## })
