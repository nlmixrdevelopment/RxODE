library(RxODE)
library(data.table);
library(digest)
rxPermissive({

    context("Make sure solved rxSolve object can be converted to a data.table easily.")

    ## RxODE instance 1
    m1 <-
        RxODE(
            model = '
         C2 = centr/V2;
         C3 = peri/V3;
         d/dt(depot) =-KA*depot;
         d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
         d/dt(peri)  =                    Q*C2 - Q*C3;
         d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff;')

    test_that("RxODE instance 1 is created",{
        expect_equal(class(m1),"RxODE");
    })

    et1 <- eventTable(amount.units="ug", time.units = "hours")
    et1$add.dosing(dose=10000, nbr.doses=5, dosing.interval = 24)
    et1$add.sampling(0:24)
    et1$add.sampling(24);
    et1$add.sampling(seq(from = 24+8, to = 5*24, by = 8))


    test_that("RxODE event table 1 was created",{
        expect_equal(class(et1), "EventTable")
        expect_equal(et1$get.nobs(),38);
        expect_equal(length(et1$get.dosing()[,1]), 5);
    })

    o1.first <- rxSolve(m1, params = c(KA=.291, CL=18.6, V2=40.2, Q=10.5, V3=297.0,
                                       Kin=1.0, Kout=1.0, EC50=200.0),
                        events = et1,
                        inits = c(0, 0, 0, 1))

    test_that("Can convert to data.table using as.data.table", {
        expect_equal(class(as.data.table(o1.first)), c("data.table", "data.frame"))
    });
})
