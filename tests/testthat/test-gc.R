library(RxODE)
library(testthat)

rxPermissive({
    context("Garbage collection")
    test_that("Check garbage collection unloads DLLs", {

        ode <- RxODE({
            b       = -1
            d/dt(X) = a*X + Y*Z;
            d/dt(Y) = b*(Y - Z);
            d/dt(Z) = -X*Y + c*Y - Z
        })

        name <- basename(rxDll(ode))
        name <- substr(name, 0, nchar(name) - nchar(.Platform$dynlib.ext))

        expect_false(is.null(getLoadedDLLs()[[name]]))
        Sys.sleep(0.5)
        gc()
        Sys.sleep(0.5)
        expect_false(is.null(getLoadedDLLs()[[name]]))
        rm(ode)
        Sys.sleep(0.5)
        gc()
        Sys.sleep(0.5)
        expect_true(is.null(getLoadedDLLs()[[name]]))

        ## Now test locking and unlocking.

        ode <- RxODE({
            b       = -1
            d/dt(X) = a*X + Y*Z;
            d/dt(Y) = b*(Y - Z);
            d/dt(Z) = -X*Y + c*Y - Z
        })

        dll <- rxDll(ode);
        rxDynProtect(dll);
        name <- basename(rxDll(ode))
        name <- substr(name, 0, nchar(name) - nchar(.Platform$dynlib.ext))

        expect_false(is.null(getLoadedDLLs()[[name]]))
        Sys.sleep(0.5)
        gc()
        Sys.sleep(0.5)
        expect_false(is.null(getLoadedDLLs()[[name]]))
        rm(ode)
        Sys.sleep(0.5)
        gc()
        Sys.sleep(0.5)
        expect_false(is.null(getLoadedDLLs()[[name]]))

        ode2 <- RxODE({
            b       = -1
            d/dt(X) = a*X + Y*Z + X;
            d/dt(Y) = b*(Y - Z);
            d/dt(Z) = -X*Y + c*Y - Z
        })

        ## Now unlock protection
        rxDynProtect("");

        dll <- rxDll(ode2);
        name2 <- basename(rxDll(ode2))
        name2 <- substr(name2, 0, nchar(name2) - nchar(.Platform$dynlib.ext))

        expect_false(is.null(getLoadedDLLs()[[name2]]))
        Sys.sleep(0.5)
        gc()
        Sys.sleep(0.5)
        expect_false(is.null(getLoadedDLLs()[[name2]]))
        rm(ode2)
        Sys.sleep(0.5)
        gc()
        Sys.sleep(0.5)
        ## Unloads BOTH dlls since it has been unprotected
        expect_true(is.null(getLoadedDLLs()[[name]]))
        expect_true(is.null(getLoadedDLLs()[[name2]]))

    })

})
