library(RxODE);
context("Test Inis");
rxPermissive({
    out <- RxODE("ini = 1; fun_ini = 2; fun = 4; addit = ini + fun_ini + pi + no_ini")

    test_that("Initial constants are correct",{
        expect_equal(rxInits(out)["pi"],c(pi=pi));
        expect_equal(rxInits(out)["ini"],c(ini=1));
        expect_equal(rxInits(out)["fun_ini"],c(fun_ini=2));
        expect_equal(rxInits(out)["fun"],c(fun=4));
    })

    test_that("Constants are not included in get.modelVars()",{
        expect_equal(out$get.modelVars()$params,c("no_ini"));
    })
    options(RxODE.syntax.allow.ini=FALSE);
    out2 <- RxODE("ini = 1; fun_ini = 2; fun = 4; addit = ini + fun_ini + pi + no_ini")

    test_that("Initial constants only include pi.", {
        expect_equal(names(rxInits(out2)), "pi");
        expect_true(rxDll(out) != rxDll(out2));
        expect_equal(out2$get.modelVars()$params,c("no_ini"));
    })
    options(RxODE.syntax.allow.ini=TRUE);

    ## out <- RxODE({
    ##     theta[1] = 3
    ##     eta[1] = 2
    ##     k = exp(theta[1] + eta[1])
    ##     d / dt(central) = -theta[1] * central
    ## })

    ## test_that("Allow THETA[#] and ETA[#]s.", {
    ## })
}, silent=TRUE)
