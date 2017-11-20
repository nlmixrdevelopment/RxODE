rxPermissive({
    context("Parameter and Data translation")
    test_that("par and data", {
        mod <- RxODE({
            a = 6
            b = 0.6
            d/dt(intestine) = -a*intestine
            d/dt(blood)     = a*intestine - b*blood
        })

        et <- eventTable(time.units="days")
        et$add.sampling(seq(0,10,by=1/24))
        et$add.dosing(dose=2/24,rate=2,strt.time=0,
                      nbr.doses=10,dosing.interval=1)

        tmp2 <- rxDataParSetup(mod, et);
        expect_equal(class(tmp2), c("RxODE.par.data", "RxODE.multi.data"));
        expect_equal(tmp2$pars, c(6, 0.6))
        expect_equal(tmp2$nsim, 1L)
        expect_equal(tmp2$n.pars, 2L)
        expect_equal(tmp2$inits, structure(c(0, 0), .Names = c("intestine", "blood")));
    })
}, cran=FALSE)
