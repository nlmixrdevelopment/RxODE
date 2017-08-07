library(RxODE)
rxPermissive({
    rxClean()

    context("Simple test for time-varying covariates")

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
                   covs = cov, add.cov=TRUE);
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

    ## covs_interpolation
    sink("temp.csv");
    cat("t,c\n");
    out <- rxSolve(ode,
                   params = c(a=-8/3, b=-10),
                   events = et,
                   inits = c(X=1, Y=1, Z=1),
                   covs = cov,
                   covs_interpolation="constant", add.cov=TRUE);
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
                inits = c(X=1, Y=1, Z=1), add.cov=TRUE)
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
                   covs = cov, add.cov=TRUE)

    out3 <- rxSolve(ode,
                    params = c(a=-8/3, b=-10),
                    events = et,
                    inits = c(X=1, Y=1, Z=1),
                    covs = cov, add.cov=TRUE)
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
                    covs = cov, add.cov=TRUE)
    sink();
    unlink("temp");

    test_that("Before assinging the time varying to -8/3, out and out2 should be different",{
        expect_false(isTRUE(all.equal(out,out2)));
    })

    context("Test First Assignment")
    ## Assign a time-varying to a simple parameter
    sink("temp");
    out$a <- -8/3
    sink();
    unlink("temp");

    test_that("The out$a=-8/3 works.",{
        expect_equal(as.data.frame(out),as.data.frame(out2));
    })


    context("Test Second Assignment")
    sink("temp")
    out$a <- out3$a
    sink();
    unlink("temp");

    test_that("the out$a = time varying covariate works.",{
        expect_equal(as.data.frame(out),as.data.frame(out3));
    })

    rxClean()
}, silent=TRUE, cran=TRUE)
