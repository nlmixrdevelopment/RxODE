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

        expect_equal(class(RxODE:::rxSolvingData(mod,tmp2)), "externalptr")

        expect_equal(class(tmp2), c("RxODE.par.data", "RxODE.multi.data"));
        expect_equal(tmp2$pars, c(6, 0.6, 0))
        expect_equal(tmp2$nsim, 1L)
        expect_equal(tmp2$n.pars, 2L)
        expect_equal(tmp2$inits, structure(c(0, 0), .Names = c("intestine", "blood")));

        ## Now try numeric vector setup
        tmp2 <- rxDataParSetup(mod, et, c(a=60));
        expect_equal(tmp2$pars, c(60, 0.6, 0))

        tmp2 <- rxDataParSetup(mod, c(a=60), et);
        expect_equal(tmp2$pars, c(60, 0.6, 0))

        tmp2 <- rxDataParSetup(mod, et, c(b=60));
        expect_equal(tmp2$pars, c(6, 60, 0))

        tmp2 <- rxDataParSetup(mod, c(b=60), et);
        expect_equal(tmp2$pars, c(6, 60, 0))

        ## Now try parameters that don't exist
        tmp2 <- rxDataParSetup(mod, et, c(c=60));
        expect_equal(tmp2$pars, c(6, 0.6, 0))

        tmp2 <- rxDataParSetup(mod, c(c=60), et);
        expect_equal(tmp2$pars, c(6, 0.6, 0))

        ## Now do a data frame.
        tmp2 <- rxDataParSetup(mod, data.frame(a=c(6, 7), b=c(8, 9)), et);
        expect_equal(tmp2$pars, c(6, 8, 0, 7, 9, 4.94065645841247e-324))

        tmp2 <- rxDataParSetup(mod, data.frame(a=c(6, 7, 8), b=c(8, 9, 10)), et);
        expect_equal(tmp2$pars, c(6, 8, 0, 7, 9, 4.94065645841247e-324, 8, 10, 9.88131291682493e-324))

        tmp2 <- rxDataParSetup(mod, et, data.frame(a=c(6, 7, 8), b=c(8, 9, 10)));
        expect_equal(tmp2$pars, c(6, 8, 0, 7, 9, 4.94065645841247e-324, 8, 10, 9.88131291682493e-324))

        tmp2 <- rxDataParSetup(mod, et, data.frame(a=c(6, 7, 8)));
        expect_equal(tmp2$pars, c(6, 0.6, 0, 7, 0.6, 4.94065645841247e-324, 8, 0.6, 9.88131291682493e-324))

        tmp2 <- rxDataParSetup(mod, et, data.frame(b=c(6, 7, 8)));
        expect_equal(tmp2$pars, c(6, 6, 0, 6, 7, 4.94065645841247e-324, 6, 8, 9.88131291682493e-324))

        tmp2 <- rxDataParSetup(mod, et, data.frame(c=c(6, 7, 8)));
        expect_equal(tmp2$pars, c(6, 0.6, 0, 6, 0.6, 4.94065645841247e-324, 6, 0.6, 9.88131291682493e-324))

        ## Now matrices
        tmp2 <- rxDataParSetup(mod, as.matrix(data.frame(a=c(6, 7), b=c(8, 9))), et);
        expect_equal(tmp2$pars, c(6, 8, 0, 7, 9, 4.94065645841247e-324))

        tmp2 <- rxDataParSetup(mod, as.matrix(data.frame(a=c(6, 7, 8), b=c(8, 9, 10))), et);
        expect_equal(tmp2$pars, c(6, 8, 0, 7, 9, 4.94065645841247e-324, 8, 10, 9.88131291682493e-324))

        tmp2 <- rxDataParSetup(mod, et, as.matrix(data.frame(a=c(6, 7, 8), b=c(8, 9, 10))));
        expect_equal(tmp2$pars, c(6, 8, 0, 7, 9, 4.94065645841247e-324, 8, 10, 9.88131291682493e-324))

        tmp2 <- rxDataParSetup(mod, et, as.matrix(data.frame(a=c(6, 7, 8))));
        expect_equal(tmp2$pars, c(6, 0.6, 0, 7, 0.6, 4.94065645841247e-324, 8, 0.6, 9.88131291682493e-324))

        tmp2 <- rxDataParSetup(mod, et, as.matrix(data.frame(b=c(6, 7, 8))));
        expect_equal(tmp2$pars, c(6, 6, 0, 6, 7, 4.94065645841247e-324, 6, 8, 9.88131291682493e-324))

        tmp2 <- rxDataParSetup(mod, et, as.matrix(data.frame(c=c(6, 7, 8))));
        expect_equal(tmp2$pars, c(6, 0.6, 0, 6, 0.6, 4.94065645841247e-324, 6, 0.6, 9.88131291682493e-324))
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
        id0 <- c(0, 4.94065645841247e-324, 9.88131291682493e-324, 1.48219693752374e-323, 1.97626258336499e-323, 2.47032822920623e-323, 2.96439387504748e-323, 3.45845952088873e-323, 3.95252516672997e-323, 4.44659081257122e-323, 4.94065645841247e-323, 5.43472210425371e-323, 5.92878775009496e-323, 6.4228533959362e-323, 6.91691904177745e-323, 7.4109846876187e-323, 7.90505033345994e-323, 8.39911597930119e-323, 8.89318162514244e-323, 9.38724727098368e-323, 9.88131291682493e-323, 1.03753785626662e-322, 1.08694442085074e-322, 1.13635098543487e-322, 1.18575755001899e-322, 1.23516411460312e-322, 1.28457067918724e-322, 1.33397724377137e-322, 1.38338380835549e-322, 1.43279037293961e-322, 1.48219693752374e-322, 1.53160350210786e-322, 1.58101006669199e-322, 1.63041663127611e-322, 1.67982319586024e-322, 1.72922976044436e-322, 1.77863632502849e-322, 1.82804288961261e-322, 1.87744945419674e-322, 1.92685601878086e-322, 1.97626258336499e-322, 2.02566914794911e-322, 2.07507571253324e-322, 2.12448227711736e-322, 2.17388884170148e-322, 2.22329540628561e-322, 2.27270197086973e-322, 2.32210853545386e-322, 2.37151510003798e-322, 2.42092166462211e-322, 2.47032822920623e-322, 2.51973479379036e-322, 2.56914135837448e-322, 2.61854792295861e-322, 2.66795448754273e-322, 2.71736105212686e-322, 2.76676761671098e-322, 2.81617418129511e-322, 2.86558074587923e-322, 2.91498731046335e-322, 2.96439387504748e-322, 3.0138004396316e-322, 3.06320700421573e-322, 3.11261356879985e-322, 3.16202013338398e-322, 3.2114266979681e-322, 3.26083326255223e-322, 3.31023982713635e-322, 3.35964639172048e-322, 3.4090529563046e-322, 3.45845952088873e-322, 3.50786608547285e-322, 3.55727265005698e-322, 3.6066792146411e-322, 3.65608577922522e-322, 3.70549234380935e-322, 3.75489890839347e-322, 3.8043054729776e-322, 3.85371203756172e-322, 3.90311860214585e-322, 3.95252516672997e-322, 4.0019317313141e-322, 4.05133829589822e-322, 4.10074486048235e-322, 4.15015142506647e-322, 4.1995579896506e-322, 4.24896455423472e-322, 4.29837111881884e-322, 4.34777768340297e-322, 4.39718424798709e-322, 4.44659081257122e-322, 4.49599737715534e-322, 4.54540394173947e-322, 4.59481050632359e-322, 4.64421707090772e-322, 4.69362363549184e-322, 4.74303020007597e-322, 4.79243676466009e-322, 4.84184332924422e-322, 4.89124989382834e-322, 4.94065645841247e-322, 4.99006302299659e-322, 5.03946958758071e-322, 5.08887615216484e-322, 5.13828271674896e-322, 5.18768928133309e-322, 5.23709584591721e-322, 5.28650241050134e-322, 5.33590897508546e-322, 5.38531553966959e-322, 5.43472210425371e-322, 5.48412866883784e-322, 5.53353523342196e-322, 5.58294179800609e-322, 5.63234836259021e-322, 5.68175492717434e-322, 5.73116149175846e-322, 5.78056805634258e-322, 5.82997462092671e-322, 5.87938118551083e-322);
        pv <- as.vector(t(data.frame(a=seq(2, 60, length.out=120),
                                     b=seq(0.6, 36, length.out=120),
                                     id=id0)))

        tmp2 <- rxDataParSetup(mod, dat,
                               data.frame(a=seq(2, 60, length.out=120),
                                          b=seq(0.6, 36, length.out=120)));

        expect_equal(tmp2$pars, pv);
        tmp2 <- rxDataParSetup(mod, data.frame(b=seq(0.6, 36, length.out=120),
                                               a=seq(2, 60, length.out=120)), dat)
        expect_equal(tmp2$pars, pv);

        ##
        tmp2 <- rxDataParSetup(mod, dat,
                               as.matrix(data.frame(a=seq(2, 60, length.out=120),
                                                    b=seq(0.6, 36, length.out=120))));

        expect_equal(tmp2$pars, pv);
        tmp2 <- rxDataParSetup(mod, as.matrix(data.frame(b=seq(0.6, 36, length.out=120),
                                                         a=seq(2, 60, length.out=120))), dat)
        expect_equal(tmp2$pars, pv);

        ## Now only include /some/ of the parameters

        pvB <- as.vector(t(data.frame(a=6, b=seq(0.6, 36, length.out=120), id=id0)))

        tmp2 <- rxDataParSetup(mod, dat,
                               data.frame(c=seq(2, 60, length.out=120),
                                          b=seq(0.6, 36, length.out=120)));

        expect_equal(tmp2$pars, pvB);
        tmp2 <- rxDataParSetup(mod, data.frame(b=seq(0.6, 36, length.out=120),
                                               c=seq(2, 60, length.out=120)), dat)
        expect_equal(tmp2$pars, pvB);

        ##
        tmp2 <- rxDataParSetup(mod, dat,
                               as.matrix(data.frame(c=seq(2, 60, length.out=120),
                                                    b=seq(0.6, 36, length.out=120))));

        expect_equal(tmp2$pars, pvB);
        tmp2 <- rxDataParSetup(mod, as.matrix(data.frame(b=seq(0.6, 36, length.out=120),
                                                         c=seq(2, 60, length.out=120))), dat)
        expect_equal(tmp2$pars, pvB);

    })
}, cran=FALSE)
