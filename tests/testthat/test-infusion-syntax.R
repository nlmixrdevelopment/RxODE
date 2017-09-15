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
        d/dt(X) = a*X + Y*Z + rxRate(X);
        d/dt(Y) = b*(Y - Z) + rxRate(Y);
        d/dt(Z) = -X*Y + c*Y - Z + rxRate(Z);
    })

    test_that("No infusions are in the ODE.", {
        expect_equal(structure("__DDtStateVar__[0] = _InfusionRate(0) + a*X+Y*Z+0.0;\n__DDtStateVar__[1] = _InfusionRate(1) + b*(Y-Z)+0.0;\n__DDtStateVar__[2] = _InfusionRate(2) + -X*Y+c*Y-Z+0.0;\n", .Names = "parseModel"),
                     rxModelVars(rx2)$model["parseModel"])
    })

    rx3 <- RxODE({
        d/dt(X) = a*X + Y*Z + rxRate(X);
        d/dt(Y) = b*(Y - Z) + rxRate(X);
        d/dt(Z) = -X*Y + c*Y - Z + rxRate(X);
    })

    test_that("Can reference prior infusions in current d/dt()", {
        expect_equal(structure("__DDtStateVar__[0] = _InfusionRate(0) + a*X+Y*Z+0.0;\n__DDtStateVar__[1] = _InfusionRate(1) + b*(Y-Z)+_InfusionRate(0);\n__DDtStateVar__[2] = _InfusionRate(2) + -X*Y+c*Y-Z+_InfusionRate(0);\n", .Names = "parseModel"), rxModelVars(rx3)$model["parseModel"])
    })

    test_that("rate(Y) before d/dt(Y) throws errors",
    {
        expect_error(RxODE({
            d/dt(X) = a*X + Y*Z + rxRate(Y);
            d/dt(Y) = b*(Y - Z) + rxRate(X);
            d/dt(Z) = -X*Y + c*Y - Z + rxRate(X);
        }))})

    ## Now test DSL handling
    test_that("rate(Y) translates to python/sympy correctly.", {
        expect_equal(structure("rx__d_dt_depot__ = rxRate(depot) - ka * depot", .Names = "rx__d_dt_depot__"),
                     rxToSymPy("d/dt(depot)=rxRate(depot)-ka*depot"));
        expect_equal(structure("d/dt(depot) = rxRate(depot) - ka * depot", .Names = "d/dt(depot)"),
                     rxFromSymPy("rx__d_dt_depot__ = rxRate(depot) - ka * depot"))
    })

}, silent=TRUE, on.validate=TRUE);
