library(RxODE)
rxPermissive({
    rxClean()
    for (meth in c("liblsoda", "lsoda")){ ## Dop is very close but doesn't match precisely.
        context(sprintf("Simple test for time-varying covariates (%s)", meth))

        ode <- RxODE({
            b       = -1
            d/dt(X) = a*X + Y*Z;
            d/dt(Y) = b*(Y - Z);
            d/dt(Z) = -X*Y + c*Y - Z
            printf("%.10f,%.10f\n",t,c);
        })

        et <- eventTable()   # default time units
        et$add.sampling(seq(from=0, to=100, by=0.01))

        cov <- data.frame(c=et$get.sampling()$time+1);

        cov.lin <- approxfun(et$get.sampling()$time, cov$c, yleft=cov$c[1], yright=cov$c[length(cov$c)]);

        sink("temp.csv");
        cat("t,c\n");
        out <- rxSolve(ode,
                       params = c(a=-8/3, b=-10),
                       events = et,
                       inits = c(X=1, Y=1, Z=1),

                       covs = cov, add.cov=TRUE,
                       covsInterpolation="linear",
                       method=meth);
        sink();
        lin.interp <- read.csv("temp.csv");
        unlink("temp.csv");

        lin.interp$c2 <- cov.lin(lin.interp$t);

        test_that("time varying covariates output covariate in data frame",{
            expect_equal(cov$c,out$c);
        })

        test_that("Linear Approximation matches approxfun.", {
            expect_equal(lin.interp$c, lin.interp$c2);
        })

        ## NONMEM interpolation
        sink("temp.csv");
        cat("t,c\n");
        out <- rxSolve(ode,
                       params = c(a=-8/3, b=-10),
                       events = et,
                       inits = c(X=1, Y=1, Z=1),
                       covs = cov,
                       covs_interpolation="NOCB", add.cov=TRUE,
                       method=meth);
        sink()
        lin.interp <- read.csv("temp.csv");
        unlink("temp.csv");

        cov.lin <- approxfun(out$time, out$c, yleft=cov$c[1], yright=cov$c[length(cov$c)],
                             method="constant", f=1);
        lin.interp$c2 <- cov.lin(lin.interp$t);

        test_that("NOCB Approximation similar to approxfun.", {
            expect_equal(lin.interp$c, lin.interp$c2);
        })


        ## midpoint interpolation
        sink("temp.csv");
        cat("t,c\n");
        out <- rxSolve(ode,
                       params = c(a=-8/3, b=-10),
                       events = et,
                       inits = c(X=1, Y=1, Z=1),
                       covs = cov,
                       covs_interpolation="midpoint", add.cov=TRUE,
                       method=meth);
        sink()
        lin.interp <- read.csv("temp.csv");
        unlink("temp.csv");

        cov.lin <- approxfun(out$time, out$c, yleft=cov$c[1], yright=cov$c[length(cov$c)],
                             method="constant", f=0.5);

        lin.interp$c2 <- cov.lin(lin.interp$t);

        test_that("midpoint Approximation similar to approxfun.", {
            expect_equal(lin.interp$c, lin.interp$c2);
        })


        ## covs_interpolation
        sink("temp.csv");
        cat("t,c\n");
        out <- rxSolve(ode,
                       params = c(a=-8/3, b=-10),
                       events = et,
                       inits = c(X=1, Y=1, Z=1),
                       covs = cov,
                       covs_interpolation="constant", add.cov=TRUE,
                       method=meth);
        sink()
        lin.interp <- read.csv("temp.csv");
        unlink("temp.csv");

        cov.lin <- approxfun(out$time, out$c, yleft=cov$c[1], yright=cov$c[length(cov$c)],
                             method="constant");

        lin.interp$c2 <- cov.lin(lin.interp$t);

        test_that("Constant Approximation similar to approxfun.", {
            expect_equal(lin.interp$c, lin.interp$c2);
        })


        out <- as.data.frame(out);
        out <- out[,names(out) != "c"];

        sink("temp");
        out1 <-
            rxSolve(ode,
                    params = c(a=-8/3, b=-10, c = 0),
                    events = et,
                    inits = c(X=1, Y=1, Z=1), add.cov=TRUE,
                    method=meth)
        sink();
        unlink("temp");

        out1 <- as.data.frame(out1);

        test_that("time varying covariates produce different outputs",{
            expect_false(isTRUE(all.equal(out,out1)));
        })

        cov <- data.frame(c=et$get.sampling()$time+1,a=-et$get.sampling()$time/100);

        sink("temp")
        out <- rxSolve(ode,
                       params = c(a=-8/3, b=-10),
                       events = et,
                       inits = c(X=1, Y=1, Z=1),
                       covs = cov, add.cov=TRUE,
                       method=meth)

        out3 <- rxSolve(ode,
                        params = c(a=-8/3, b=-10),
                        events = et,
                        inits = c(X=1, Y=1, Z=1),
                        covs = cov, add.cov=TRUE,
                        method=meth)
        sink();
        unlink("temp");

        test_that("time varying covariates output covariate in data frame",{
            expect_equal(cov$c,out$c);
            expect_equal(cov$a,out$a);
        })


        cov <- data.frame(c=et$get.sampling()$time+1);

        sink("temp");
        out2 <- rxSolve(ode,
                        params = c(a=-8/3, b=-10),
                        events = et,
                        inits = c(X=1, Y=1, Z=1),
                        covs = cov, add.cov=TRUE,
                        method=meth)
        sink();
        unlink("temp");

        test_that("Before assinging the time varying to -8/3, out and out2 should be different",{
            expect_false(isTRUE(all.equal(out,out2)));
        })

        context(sprintf("Test First Assignment (%s)", meth))
        ## Assign a time-varying to a simple parameter
        sink("temp");
        out$a <- -8/3
        sink();
        unlink("temp");

        test_that("The out$a=-8/3 works.",{
            expect_equal(as.data.frame(out),as.data.frame(out2));
        })


        context(sprintf("Test Second Assignment (%s)", meth))
        sink("temp")
        out$a <- out3$a
        sink();
        unlink("temp");

        test_that("the out$a = time varying covariate works.",{
            expect_equal(as.data.frame(out),as.data.frame(out3));
        })

        context(sprintf("Covariate solve with data frame event table (%s)", meth))

        ## Covariate solve for data frame
        d3 <- structure(list(TIME = c(0, 0, 2.99270072992701, 192, 336, 456),
                             AMT = c(137L, 0L, -137L, 0L, 0L, 0L),
                             V2I = c(909L, 909L, 909L, 909L, 909L, 909L),
                             V1I = c(545L, 545L, 545L, 545L, 545L, 545L),
                             CLI = c(471L, 471L, 471L, 471L, 471L, 471L),
                             EVID = c(10101L, 0L, 10101L, 0L, 0L, 0L)),
                        class = "data.frame",
                        row.names = c(NA,  -6L),
                        .Names = c("TIME", "AMT", "V2I", "V1I", "CLI", "EVID"))
        mod1 = RxODE({
            d/dt(A_centr)=-A_centr*(CLI/V1I+204/V1I)+204*A_periph/V2I;
            d/dt(A_periph)=204*A_centr/V1I-204*A_periph/V2I;
            d/dt(A_circ)=-4*A_circ*exp(-ETA[2]-THETA[2])+4*A_tr3*exp(-ETA[2]-THETA[2]);
            A_circ(0)=exp(ETA[1]+THETA[1]);
            d/dt(A_prol)=4*A_prol*Rx_pow(exp(ETA[1]+THETA[1])/A_circ,exp(THETA[4]))*(-A_centr*exp(ETA[3]+THETA[3])/V1I+1)*exp(-ETA[2]-THETA[2])-4*A_prol*exp(-ETA[2]-THETA[2]);
            A_prol(0)=exp(ETA[1]+THETA[1]);
            d/dt(A_tr1)=4*A_prol*exp(-ETA[2]-THETA[2])-4*A_tr1*exp(-ETA[2]-THETA[2]);
            A_tr1(0)=exp(ETA[1]+THETA[1]);
            d/dt(A_tr2)=4*A_tr1*exp(-ETA[2]-THETA[2])-4*A_tr2*exp(-ETA[2]-THETA[2]);
            A_tr2(0)=exp(ETA[1]+THETA[1]);
            d/dt(A_tr3)=4*A_tr2*exp(-ETA[2]-THETA[2])-4*A_tr3*exp(-ETA[2]-THETA[2]);
            A_tr3(0)=exp(ETA[1]+THETA[1]);
        })
        tmp <- rxSolve(mod1, d3, structure(c(2.02103, 4.839305, 3.518676, -1.391113, 0.108127023, -0.064170725, 0.087765769),
                                           .Names=c(sprintf("THETA[%d]", 1:4), sprintf("ETA[%d]", 1:3))), add.cov=TRUE,
                       method=meth)

        library(dplyr)
        test_that("Data Frame single subject solve", {
            expect_equal(tmp %>% select(CLI, V1I, V2I) %>% as.data.frame, d3 %>% filter(EVID == 0) %>% select(CLI,V1I,V2I) %>% as.data.frame)
            expect_equal(names(tmp$params), mod1$params[-(1:3)])
        })

        d3 <- structure(list(ID = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                             TIME = c(0, 0, 2.99270072992701, 192, 336, 456, 0, 0, 3.07272727272727, 432),
                             AMT = c(137L, 0L, -137L, 0L, 0L, 0L, 110L, 0L, -110L, 0L),
                             V2I = c(909L, 909L, 909L, 909L, 909L, 909L, 942L, 942L, 942L, 942L),
                             V1I = c(545L, 545L, 545L, 545L, 545L, 545L, 306L, 306L, 306L, 306L),
                             CLI = c(471L, 471L, 471L, 471L, 471L, 471L, 405L, 405L, 405L, 405L),
                             EVID = c(10101L, 0L, 10101L, 0L, 0L, 0L, 10101L, 0L, 10101L, 0L)),
                        class = "data.frame", row.names = c(NA,  -10L),
                        .Names = c("ID", "TIME", "AMT", "V2I", "V1I", "CLI", "EVID"))


        par2 <- matrix(c(2.02103, 4.839305, 3.518676, -1.391113, 0.108127023, -0.064170725, 0.087765769,
                         2.02103, 4.839305, 3.518676, -1.391113,  -0.064170725, 0.087765769,0.108127023), nrow=2, byrow=T,
                       dimnames=list(NULL, c(sprintf("THETA[%d]", 1:4), sprintf("ETA[%d]", 1:3))));

        tmp <- rxSolve(mod1, d3, par2, add.cov=TRUE, cores=2, method=meth)

        test_that("Data Frame multi subject solve", {
            expect_equal(tmp %>% select(CLI, V1I, V2I) %>% as.data.frame, d3 %>% filter(EVID == 0) %>% select(CLI,V1I,V2I) %>% as.data.frame)
            expect_equal(names(tmp$params)[-1], mod1$params[-(1:3)])
        })

        ## Now check missing covariate values.

        d3na <- structure(list(ID = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                             TIME = c(0, 0, 2.99270072992701, 192, 336, 456, 0, 0, 3.07272727272727, 432),
                             AMT = c(137L, 0L, -137L, 0L, 0L, 0L, 110L, 0L, -110L, 0L),
                             V2I = c(909L, NA_integer_, 909L, 909L, 909L, 909L, 942L, 942L, 942L, 942L),
                             V1I = c(545L, 545L, 545L, 545L, 545L, 545L, 306L, 306L, 306L, NA_integer_),
                             CLI = c(471L, 471L, 471L, 471L, NA_integer_, 471L, 405L, 405L, 405L, 405L),
                             EVID = c(10101L, 0L, 10101L, 0L, 0L, 0L, 10101L, 0L, 10101L, 0L)),
                        class = "data.frame", row.names = c(NA,  -10L),
                        .Names = c("ID", "TIME", "AMT", "V2I", "V1I", "CLI", "EVID"))

        tmp <- rxSolve(mod1, d3na, par2, add.cov=TRUE, cores=2, method=meth)

        tmp2 <- rxSolve(mod1, d3, par2, add.cov=TRUE, cores=2, method=meth)

        context(sprintf("Test NA extrapolation for %s solving", meth))
        test_that("NA solve is the same", {
            for (i in c("id", "time", "A_centr", "A_periph", "A_circ", "A_prol", "A_tr1", "A_tr2", "A_tr3")){
                expect_equal(tmp[[i]],tmp2[[i]])
            }
        })

        d3na <- structure(list(ID = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                               TIME = c(0, 0, 2.99270072992701, 192, 336, 456, 0, 0, 3.07272727272727, 432),
                               AMT = c(137L, 0L, -137L, 0L, 0L, 0L, 110L, 0L, -110L, 0L),
                               V2I = c(909L, NA_integer_, 909L, 909L, 909L, 909L, 942L, 942L, 942L, 942L),
                               V1I = c(545L, 545L, 545L, 545L, 545L, 545L, NA_integer_, NA_integer_, NA_integer_, NA_integer_),
                               CLI = c(471L, 471L, 471L, 471L, NA_integer_, 471L, 405L, 405L, 405L, 405L),
                               EVID = c(10101L, 0L, 10101L, 0L, 0L, 0L, 10101L, 0L, 10101L, 0L)),
                          class = "data.frame", row.names = c(NA,  -10L),
                          .Names = c("ID", "TIME", "AMT", "V2I", "V1I", "CLI", "EVID"))

        test_that("All covariates are NA give a warning",{
                    expect_warning(rxSolve(mod1, d3na, par2, add.cov=TRUE, cores=2, method=meth),"One or more covariates were all NA for subject id=2")
        })

        ## Now test non time-varying covariates

        mod <- RxODE({
            tWt = WT
            k10 = (CL*(WT/70)^0.75)/V2
            k12 = Q/V2
            k21 = Q/V3
            d/dt(depot) =-KA*depot;
            d/dt(centr) = KA*depot - k10*centr - k12*centr + k21*peri;
            d/dt(peri)  =                        k12*centr - k21*peri;
            C2 = centr/V2;
            C3 = peri/V3;
            cp = C2 #*(1+cp.err)
        })

        d <- structure(list(ID = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                                   2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
                                   3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L),
                            TIME = c(0, 0.302, 2.908,
                                     3.144, 9.943, 18.021, 24, 48, 72, 96, 120, 144, 168, 168.009,
                                     169.948, 177.614, 188.3, 0, 0.438, 2.529, 4.904, 7.498, 21.243,
                                     24, 48, 72, 96, 120, 144, 168, 168.958, 169.062, 171.551, 191.689,
                                     0, 0.438, 2.529, 4.904, 7.498, 21.243, 24, 48, 72, 96, 120, 144,
                                     168, 168.958, 169.062, 171.551, 191.689),
                            DV = c(0, 155.6, 325.2,
                                   346.4, 166.9, 101, 0, 0, 0, 0, 0, 0, 0, 286, 647.2, 444.5, 354.8,
                                   0, 203.5, 382.9, 280.6, 167.3, 73.2, 0, 0, 0, 0, 0, 0, 0, 486.8,
                                   545.1, 523.5, 214, 0, 196.5, 462.9, 394, 247.4, 93.8, 0, 0, 0,
                                   0, 0, 0, 0, 680.7, 761, 828.2, 449.1),
                            WT = c(52L, 52L, 52L,
                                   52L, 52L, 52L, 52L, 52L, 52L, 52L, 52L, 52L, 52L, 52L, 52L, 52L,
                                   52L, 73L, 73L, 73L, 73L, 73L, 73L, 73L, 73L, 73L, 73L, 73L, 73L,
                                   73L, 73L, 73L, 73L, 73L, 53L, 53L, 53L, 53L, 53L, 53L, 53L, 53L,
                                   53L, 53L, 53L, 53L, 53L, 53L, 53L, 53L, 53L),
                            SEX = c(0L, 0L,
                                    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L,
                                    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                    1L),
                            AMT = c(1200L, 0L, 0L, 0L, 0L, 0L, 1200L, 1200L, 1200L,
                                    1200L, 1200L, 1200L, 1200L, 0L, 0L, 0L, 0L, 1200L, 0L, 0L, 0L,
                                    0L, 0L, 1200L, 1200L, 1200L, 1200L, 1200L, 1200L, 1200L, 0L,
                                    0L, 0L, 0L, 1200L, 0L, 0L, 0L, 0L, 0L, 1200L, 1200L, 1200L, 1200L,
                                    1200L, 1200L, 1200L, 0L, 0L, 0L, 0L),
                            EVID = c(101L, 0L, 0L,
                                     0L, 0L, 0L, 101L, 101L, 101L, 101L, 101L, 101L, 101L, 0L, 0L,
                                     0L, 0L, 101L, 0L, 0L, 0L, 0L, 0L, 101L, 101L, 101L, 101L, 101L,
                                     101L, 101L, 0L, 0L, 0L, 0L, 101L, 0L, 0L, 0L, 0L, 0L, 101L, 101L,
                                     101L, 101L, 101L, 101L, 101L, 0L, 0L, 0L, 0L)),
                       row.names = c(NA, -51L), class = "data.frame")

        tmp <- RxODE:::evTrans(d, mod)



        tmp2 <- rxSolve(mod,
                        c(KA=1.05, CL=0.121, V2=1.939,
                          Q=0.282, V3=5.65), tmp) %>%
            as.data.frame

        tmp3 <- rxSolve(mod,
                        data.frame(KA=1.05, CL=0.121, V2=1.939,
                                   Q=0.282, V3=5.65,WT=c(52, 73, 53)), d) %>%
            as.data.frame

        test_that("non time-varying covariates are the same as supplying 1 parameter for each id",{
            expect_equal(tmp2, tmp3)
        })
    }
    ## devtools::install();library(RxODE);rxTest("cov")
    rxClean()

}, silent=TRUE, cran=TRUE)
