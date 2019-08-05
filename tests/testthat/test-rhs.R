library(RxODE)
library(dplyr)
library(digest)
rxPermissive({
    context("rxSolve right handed differental equations")

    rxSetIni0(FALSE)

    orhs <- RxODE("
x(0)    = 1
a       = 1.0E4
d/dt(x) = a*y*z - 0.04*x
d/dt(z) = 3.0e7*y^2
d/dt(y) = -1.0*(d/dt(x)+d/dt(z))
")


    et <- eventTable();
    et$add.sampling(c(4*10^seq(-1,10)))

    ## print(o1);
    test_that("rxSolve produces the correct results for rhs",{

        o1 <- expect_warning(rxSolve(orhs,et));

        ## expect_warning(rxSolve(orhs,et),"The initial conditions are at t = 0.4 instead of t = 0.")
        expect_equal(round(as.data.frame(o1),4),
                     structure(list(time = c(0.4, 4, 40, 400, 4000, 40000, 4e+05, 4e+06, 4e+07, 4e+08, 4e+09, 4e+10), x = c(1, 0.9117, 0.7168, 0.4506, 0.1832, 0.039, 0.0049, 5e-04, 1e-04, 0, 0, 0), z = c(0, 0.0882, 0.2831, 0.5494, 0.8168, 0.961, 0.9951, 0.9995, 0.9999, 1, 1, 1), y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), .Names = c("time", "x", "z", "y"), row.names = c(NA, -12L), class = "data.frame"))
        o1$add.sampling(0);

        expect_equal(round(as.data.frame(o1),4),
                     structure(list(time = c(0, 0.4, 4, 40, 400, 4000, 40000, 4e+05, 4e+06, 4e+07, 4e+08, 4e+09, 4e+10), x = c(1, 0.9852, 0.9055, 0.7158, 0.4505, 0.1832, 0.039, 0.0049, 5e-04, 1e-04, 0, 0, 0),     z = c(0, 0.0148, 0.0945, 0.2842, 0.5495, 0.8168, 0.961, 0.9951,     0.9995, 0.9999, 1, 1, 1), y = c(0, 0, 0, 0, 0, 0, 0, 0, 0,     0, 0, 0, 0)), .Names = c("time", "x", "z", "y"), row.names = c(NA, -13L), class = "data.frame"))

    })

    rxSetIni0(TRUE)

}, silent=TRUE)
