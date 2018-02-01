context("Test Parallel Solve")
rxPermissive({

    mod <- RxODE({
        d/dt(intestine) = -a*intestine
        d/dt(blood)     = a*intestine - b*blood
    });

    et <- eventTable(time.units="days")
    et$add.sampling(seq(0, 10, length.out=50))
    et$add.dosing(dose=2/24,rate=2,strt.time=0,
                  nbr.doses=10,dosing.interval=1)

    p <- data.frame(a=6,b=seq(0.4,0.9,length.out=4));

    pk1 <- rxSolve(mod,p,et,cores=1)

    pk1 <- rxSolve(mod,p,et)

    pk2 <- rxSolve(mod,p,et, cores=2); # CRAN requirement of at most 2 cores.
    test_that("Parallel Solve gives same results a single threaded solve", {
        expect_equal(as.data.frame(pk1),as.data.frame(pk2))
    })

}, silent=TRUE, cran=TRUE)
