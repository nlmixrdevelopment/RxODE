rxPermissive({
    context("Bad dosing")
    ## https://cran.r-project.org/web/packages/diffEq/vignettes/ODEinR.pdf p11

    ## 6.1

    mod <- RxODE("
a = 6
b = 0.6
d/dt(intestine) = -a*intestine
d/dt(blood)     = a*intestine - b*blood
")

    et <- eventTable(time.units="days") %>%
        add.sampling(seq(0,10,by=1/24)) %>%
        add.dosing(dose=2/24,rate=2,start.time=0, nbr.doses=10,dosing.interval=1, dosing.to=3);

    test_that("Warning for bad dose", {
        expect_warning(solve(mod, et), rex::rex("Dose to Compartment 3 ignored (not in ODE; id=1)"));
    })

    et <- eventTable(time.units="days") %>%
        add.sampling(seq(0,10,by=1/24)) %>%
        add.dosing(dose=2/24,rate=2,start.time=0, nbr.doses=10,dosing.interval=1, dosing.to=1)

    test_that("No Warning for good dose", {
        expect_warning(solve(mod, et), NA);
    })

}, silent=TRUE);
