rxPermissive({

    context("Infusion Syntax test")
    rx1 <- RxODE({
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z;
    })

    test_that("No infusions are in the ODE.", {
        expect_equal(structure("__DDtStateVar__[0] = _InfusionRate(0) + a*X+Y*Z;\n__DDtStateVar__[1] = _InfusionRate(1) + b*(Y-Z);\n__DDtStateVar__[2] = _InfusionRate(2) + -X*Y+c*Y-Z;\n", .Names = "parseModel"),
                     rxModelVars(rx1)$model["parseModel"])
    })

    rx2 <- RxODE({
        d/dt(X) = a*X + Y*Z + rate(X);
        d/dt(Y) = b*(Y - Z) + rate(Y);
        d/dt(Z) = -X*Y + c*Y - Z + rate(Z);
    })

    test_that("No infusions are in the ODE.", {
        expect_equal(structure("__DDtStateVar__[0] = _InfusionRate(0) + a*X+Y*Z+0.0;\n__DDtStateVar__[1] = _InfusionRate(1) + b*(Y-Z)+0.0;\n__DDtStateVar__[2] = _InfusionRate(2) + -X*Y+c*Y-Z+0.0;\n", .Names = "parseModel"),
                     rxModelVars(rx2)$model["parseModel"])
    })

    rx3 <- RxODE({
        d/dt(X) = a*X + Y*Z + rate(X);
        d/dt(Y) = b*(Y - Z) + rate(X);
        d/dt(Z) = -X*Y + c*Y - Z + rate(X);
    })

    test_that("Can reference prior infusions in current d/dt()", {
        expect_equal(structure("__DDtStateVar__[0] = _InfusionRate(0) + a*X+Y*Z+0.0;\n__DDtStateVar__[1] = _InfusionRate(1) + b*(Y-Z)+_InfusionRate(0);\n__DDtStateVar__[2] = _InfusionRate(2) + -X*Y+c*Y-Z+_InfusionRate(0);\n", .Names = "parseModel"), rxModelVars(rx3)$model["parseModel"])
    })

    test_that("rate(Y) before d/dt(Y) throws errors",
    {
        expect_error(RxODE({
            d/dt(X) = a*X + Y*Z + rate(Y);
            d/dt(Y) = b*(Y - Z) + rate(X);
            d/dt(Z) = -X*Y + c*Y - Z + rate(X);
        }))})

}, silent=TRUE)
