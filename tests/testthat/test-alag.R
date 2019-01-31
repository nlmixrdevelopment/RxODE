## Tests for absorption lag time.
rxPermissive({

    context("Test absorption lag-time with IV dosing")

    ## 6.1
    mod <- RxODE({
        a = 6
        b = 0.6
        d/dt(intestine) = -a*intestine
        d/dt(blood)     = a*intestine - b*blood
    })

    et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1/24))
    et$add.dosing(dose=2/24,start.time=0,
                  nbr.doses=10,dosing.interval=1)

    solve1 <- solve(mod,et)

    mod2 <- RxODE({
        a = 6
        b = 0.6
        d/dt(intestine) = -a*intestine
        alag(intestine)    = 2
        d/dt(blood)     = a*intestine - b*blood
    })

    solve2 <- solve(mod2,et)

    test_that("Solves with lag times are different", {
        expect_false(all(solve1$intestine ==solve2$intestine))
        expect_false(all(solve1$blood ==solve2$blood))
    })

    et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1/24))
    et$add.dosing(dose=2/24,start.time=2,
                  nbr.doses=10,dosing.interval=1)

    solve3 <- solve(mod, et)

    test_that("Absorption lag shifts event by 2", {
        expect_equal(solve3$intestine, solve2$intestine)
        expect_equal(solve3$blood, solve2$blood)
    })

})
