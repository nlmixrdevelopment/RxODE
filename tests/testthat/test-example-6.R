require(RxODE);
context("Example 6")
require(digest)
rxPermissive({
    ## https://cran.r-project.org/web/packages/diffEq/vignettes/ODEinR.pdf p11

    ## 6.1

    mod <- RxODE("
a = 6
b = 0.6
d/dt(intestine) = -a*intestine
d/dt(blood)     = a*intestine - b*blood
")

    et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1/24))
    et$add.dosing(dose=2/24,rate=2,strt.time=0,
                  nbr.doses=10,dosing.interval=1)

    pk <- solve(mod,et)

    ## plot(pk$time,pk$intestine,type="l")
    ## plot(pk$time,pk$blood,type="l")

    test_that("infusion model works.",{
        expect_equal(digest(round(as.data.frame(pk),4)),"e76e55f01b5da4f3dadaa23448981b4f")
    })

    ## Change QD dosing to BID dosing.
    et2 <- eventTable(time.units="days")
    et2$add.sampling(seq(0,10,by=1/24))
    et2$add.dosing(dose=2/24,rate=2,strt.time=0,
                   nbr.doses=10,dosing.interval=0.5)

    pk2 <- pk + et2;

    ## plot(pk2$time,pk2$intestine,type="l")
    ## plot(pk2$time,pk2$blood,type="l")

    test_that("bid infusion works.",{
        expect_equal(digest(round(as.data.frame(pk2),4)),"e206c7451d6f9939d029f4f1d12b8e88")
    })


    ## Bolus 6.2
    mod <- RxODE("
b = 0.6
d/dt(blood) = - b*blood
")

    et <- eventTable(time.units="days")
    et$add.sampling(seq(0,10,by=1/24));
    et$add.dosing(start.time=0,dose=40,nbr.doses=20,dosing.interval=1);

    pk <- solve(mod,et)

    ## plot(pk$time,pk$blood,type="l")

    test_that("bolus model works",{
        expect_equal(digest(round(as.data.frame(pk),4)),
                     "51d19054249014c5b925dce413cc0f9b")
    })

    et2 <- eventTable(time.units="days")
    et2$add.sampling(seq(0,10,by=1/24));
    et2$add.dosing(start.time=0,dose=40,nbr.doses=20,dosing.interval=0.5);

    pk2 <- rxChain(pk,et2);

    ## plot(pk2$time,pk2$blood,type="l")

    test_that("bolus model bid works",{
        expect_equal(digest(round(as.data.frame(pk2),4)),
                     "43990fa93e861e6194a4014fd53d2c1c")
    })

    rxClean()
}, silent=TRUE);
