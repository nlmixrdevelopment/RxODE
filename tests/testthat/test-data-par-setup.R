rxPermissive({
    context("Parameter and Data translation")
    test_that("named par and single event table", {

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

        ## Now try numeric vector setup
        tmp2 <- rxDataParSetup(mod, et, c(a=60));
        expect_equal(tmp2$pars, c(60, 0.6))

        tmp2 <- rxDataParSetup(mod, c(a=60), et);
        expect_equal(tmp2$pars, c(60, 0.6))


        tmp2 <- rxDataParSetup(mod, et, c(b=60));
        expect_equal(tmp2$pars, c(6, 60))

        tmp2 <- rxDataParSetup(mod, c(b=60), et);
        expect_equal(tmp2$pars, c(6, 60))

        ## Now try parameters that don't exist
        tmp2 <- rxDataParSetup(mod, et, c(c=60));
        expect_equal(tmp2$pars, c(6, 0.6))

        tmp2 <- rxDataParSetup(mod, c(c=60), et);
        expect_equal(tmp2$pars, c(6, 0.6))

        ## Now do a data frame.
        tmp2 <- rxDataParSetup(mod, data.frame(a=c(6, 7), b=c(8, 9)), et);
        expect_equal(tmp2$pars, c(6, 8, 7, 9))

        tmp2 <- rxDataParSetup(mod, data.frame(a=c(6, 7, 8), b=c(8, 9, 10)), et);
        expect_equal(tmp2$pars, c(6, 8, 7, 9, 8, 10))

        tmp2 <- rxDataParSetup(mod, et, data.frame(a=c(6, 7, 8), b=c(8, 9, 10)));
        expect_equal(tmp2$pars, c(6, 8, 7, 9, 8, 10))

        tmp2 <- rxDataParSetup(mod, et, data.frame(a=c(6, 7, 8)));
        expect_equal(tmp2$pars, c(6, 0.6, 7, 0.6, 8, 0.6))

        tmp2 <- rxDataParSetup(mod, et, data.frame(b=c(6, 7, 8)));
        expect_equal(tmp2$pars, c(6, 6, 6, 7, 6, 8))

        tmp2 <- rxDataParSetup(mod, et, data.frame(c=c(6, 7, 8)));
        expect_equal(tmp2$pars, c(6, 0.6, 6, 0.6, 6, 0.6))

        ## Now matrices
        tmp2 <- rxDataParSetup(mod, as.matrix(data.frame(a=c(6, 7), b=c(8, 9))), et);
        expect_equal(tmp2$pars, c(6, 8, 7, 9))

        tmp2 <- rxDataParSetup(mod, as.matrix(data.frame(a=c(6, 7, 8), b=c(8, 9, 10))), et);
        expect_equal(tmp2$pars, c(6, 8, 7, 9, 8, 10))

        tmp2 <- rxDataParSetup(mod, et, as.matrix(data.frame(a=c(6, 7, 8), b=c(8, 9, 10))));
        expect_equal(tmp2$pars, c(6, 8, 7, 9, 8, 10))

        tmp2 <- rxDataParSetup(mod, et, as.matrix(data.frame(a=c(6, 7, 8))));
        expect_equal(tmp2$pars, c(6, 0.6, 7, 0.6, 8, 0.6))

        tmp2 <- rxDataParSetup(mod, et, as.matrix(data.frame(b=c(6, 7, 8))));
        expect_equal(tmp2$pars, c(6, 6, 6, 7, 6, 8))

        tmp2 <- rxDataParSetup(mod, et, as.matrix(data.frame(c=c(6, 7, 8))));
        expect_equal(tmp2$pars, c(6, 0.6, 6, 0.6, 6, 0.6))
    })

    test_that("named par and single event table", {

        mod <- RxODE({
            a = 6
            b = 0.6
            d/dt(intestine) = -a*intestine
            d/dt(blood)     = a*intestine - b*blood
        })

        if (file.exists("test-data-setup.Rdata")){
            load("test-data-setup.Rdata")
        } else {
            tmp <- try(devtools::package_file("tests/testthat/test-data-setup.Rdata"));
            if (!inherits(tmp, "try-error") && file.exists(tmp)){
                load(tmp)
            } else {
                skip("Can't load test dataset.")
            }
        }

        pv <- as.vector(t(data.frame(a=seq(2, 60, length.out=120),
                                     b=seq(0.6, 36, length.out=120))))
        tmp2 <- rxDataParSetup(mod, dat,
                               data.frame(a=seq(2, 60, length.out=120),
                                          b=seq(0.6, 36, length.out=120)));

        expect_equal(tmp2$pars, pv);

        tmp2 <- rxDataParSetup(mod, dat,
                               data.frame(b=seq(0.6, 36, length.out=120),
                                          a=seq(2, 60, length.out=120)))

        expect_equal(tmp2$pars, pv);


    })
}, cran=FALSE)
