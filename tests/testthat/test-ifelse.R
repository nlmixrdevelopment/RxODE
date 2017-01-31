context("If/Else expansion tests")
rxPermissive({

    rx1 <- RxODE({
        b       = -1
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z
        if (t < 0.02 | t > 99.98){
            print
        } else {
            jac_print;
        }});

    test_that("Simple if/else expansion works on model", {
        expect_equal(rxExpandIfElse(rx1), c("(!(t<0.02||t>99.98))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\njac_print;",
                                            "(t<0.02||t>99.98)"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;"))
    });

    test_that("Keep initilizations if needed.", {
        expect_equal(rxExpandIfElse(rx1, FALSE), c("(!(t<0.02||t>99.98))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\njac_print;",
                                                   "(t<0.02||t>99.98)"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;"))
    })

    rx2 <- RxODE({
        b       = -1
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z
        if (t < 0.02 | t > 99.98){
            print
        }})

    test_that("Simple if expansion works on model", {
        expect_equal(rxExpandIfElse(rx2), c("(!(t<0.02||t>99.98))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                                            "(t<0.02||t>99.98)"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;"))
    });

    test_that("Keep initilizations if needed.", {
        expect_equal(rxExpandIfElse(rx2, FALSE), c("(!(t<0.02||t>99.98))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                                                   "(t<0.02||t>99.98)"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;"))
    });


    ## Compound if/then
    rx3 <- RxODE({
        b       = -1
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z
        if (t < 0.02 | t > 99.98){
            print
            if (t > 100){
                jac_print;
            } else {
                ode_print;
            }
        }});

    test_that("nested if/then", {
        expect_equal(rxExpandIfElse(rx3), c("(!(t<0.02||t>99.98))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                                            "(t<0.02||t>99.98) && (!(t>100))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\node_print;",
                                            "(t<0.02||t>99.98) && (t>100)"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;"))
    });

    test_that("Keep initilizations", {
        expect_equal(rxExpandIfElse(rx3, FALSE), c("(!(t<0.02||t>99.98))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                                                   "(t<0.02||t>99.98) && (!(t>100))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\node_print;",
                                                   "(t<0.02||t>99.98) && (t>100)"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;"))
    });

    rx4 <- RxODE({
        b       = -1
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z
        if (t < 0.02 | t > 99.98){
            print
            if (t > 100){
                jac_print;
            }
        }})

    test_that("nested if/then #2", {
        expect_equal(rxExpandIfElse(rx4), c("(!(t<0.02||t>99.98))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                                            "(t<0.02||t>99.98) && (!(t>100))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;",
                                            "(t<0.02||t>99.98) && (t>100)"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;"))
    })

    test_that("nested if/then #2 w/ini", {
        expect_equal(rxExpandIfElse(rx4, FALSE), c("(!(t<0.02||t>99.98))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                                                   "(t<0.02||t>99.98) && (!(t>100))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;",
                                                   "(t<0.02||t>99.98) && (t>100)"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;"))
    })

    ## Compound if/then
    rx5 <- RxODE({
        b       = -1
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z
        if (t < 0.02 | t > 99.98){
            print
            if (t > 100){
                jac_print;
            } else {
                ode_print;
                if (t < 0.01){
                    printf("Less than <0.001\n");
                } else {
                    printf("Between 0.01 and 0.02.\n");
                }
            }
        }
    })

    test_that("nested if/then #3", {
        expect_equal(rxExpandIfElse(rx5),
                     c("(!(t<0.02||t>99.98))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                       "(t<0.02||t>99.98) && (t>100)"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;",
                       "(t<0.02||t>99.98) && (!(t>100)) && (!(t<0.01))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\node_print;\nprintf(\"Between 0.01 and 0.02.\\n\");",
                       "(t<0.02||t>99.98) && (!(t>100)) && (t<0.01)"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\node_print;\nprintf(\"Less than <0.001\\n\");"))
    })

    test_that("nested if/then #3 w/ini", {
        expect_equal(rxExpandIfElse(rx5, FALSE),
                     c("(!(t<0.02||t>99.98))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                       "(t<0.02||t>99.98) && (t>100)"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;",
                       "(t<0.02||t>99.98) && (!(t>100)) && (!(t<0.01))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\node_print;\nprintf(\"Between 0.01 and 0.02.\\n\");",
                       "(t<0.02||t>99.98) && (!(t>100)) && (t<0.01)"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\node_print;\nprintf(\"Less than <0.001\\n\");"))
    })

    rx6 <- RxODE({
        b       = -1
        d/dt(X) = a*X + Y*Z;
        d/dt(Y) = b*(Y - Z);
        d/dt(Z) = -X*Y + c*Y - Z
        if (t < 0.02 | t > 99.98){
            print
            if (t > 100){
                jac_print;
                if (t < 150){
                    printf("Less than 150\n");
                } else {
                    printf(">= 150\n");
                }
            } else {
                ode_print;
            }
        }
    })

    test_that("nested if/then #4", {
        expect_equal(rxExpandIfElse(rx6),
                     c( "(!(t<0.02||t>99.98))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                       "(t<0.02||t>99.98) && (!(t>100))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\node_print;",
                       "(t<0.02||t>99.98) && (t>100) && (!(t<150))"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;\nprintf(\">= 150\\n\");",
                       "(t<0.02||t>99.98) && (t>100) && (t<150)"="d/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;\nprintf(\"Less than 150\\n\");")
                     )
    })

    test_that("nested if/then #4 w/ini", {
        expect_equal(rxExpandIfElse(rx6, FALSE),
                     c( "(!(t<0.02||t>99.98))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;",
                       "(t<0.02||t>99.98) && (!(t>100))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\node_print;",
                       "(t<0.02||t>99.98) && (t>100) && (!(t<150))"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;\nprintf(\">= 150\\n\");",
                       "(t<0.02||t>99.98) && (t>100) && (t<150)"="b=-1;\nd/dt(X)=a*X+Y*Z;\nd/dt(Y)=b*(Y-Z);\nd/dt(Z)=-X*Y+c*Y-Z;\nprint;\njac_print;\nprintf(\"Less than 150\\n\");")
                     )
    })

})


