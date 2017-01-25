require(RxODE);
context("Test Jacobian (df/dy) parsing")

rxPermissive({

    Vtpol2 <- RxODE("
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Jacobian
df(y)/dy(dy)  = 1
df(dy)/dy(y)  = -2*dy*mu*y - 1
df(dy)/dy(dy) = mu*(1-y^2)
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
");

    test_that("Specified jacobian is captured", {
        expect_true(any(rxDfdy(Vtpol2) == "df(y)/dy(dy)"))
        expect_true(any(rxDfdy(Vtpol2) == "df(dy)/dy(y)"))
        expect_true(any(rxDfdy(Vtpol2) == "df(dy)/dy(dy)"))
    })

    ## fixme multiple jacobian definitions raise an error.

    tmp <- RxODE("
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Jacobian
df(y)/dy(dy)  = 1
df(dy)/dy(y)  = -2*dy*mu*y - 1
df(dy)/dy(dy) = mu*(1-y^2)
df(dy)/dy(dy) = mu*(1-y^2)
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
")
    test_that("Doubled jacobain will compile correctly", {
        expect_equal(class(tmp), "RxODE");
    })

    norm <- RxODE("
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
");

    test_that("Jacobian and sensitivity not specified.", {
        expect_false(norm$calcJac);
        expect_false(norm$calcSens);
    })

    jac <- RxODE("
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
", calcJac=TRUE)

    test_that("Jacobian specified but sensitivity not specified.", {
        expect_true(jac$calcJac);
        expect_false(jac$calcSens);
    })

    sens <- RxODE("
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
", calcSens=TRUE)

    test_that("Jacobian and sensitivity specified.", {
        expect_true(sens$calcJac);
        expect_true(sens$calcSens);
    })

    rxDelete(norm);
    rxDelete(jac)
    rxDelete(sens);

    norm <- RxODE("
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
");
    jac <- RxODE(jac, calcJac=TRUE);

    test_that("Jacobian specified but sensitivity not specified.", {
        expect_true(jac$calcJac);
        expect_false(jac$calcSens);
    })

    sens <- RxODE(jac, calcSens=TRUE);

    test_that("Jacobian and sensitivity specified.", {
        expect_true(sens$calcJac);
        expect_true(sens$calcSens);
    })

    rxDelete(jac);
    rxDelete(sens);
    rxDelete(norm);

    sens <- RxODE("
d/dt(y)  = dy
d/dt(dy) = mu*(1-y^2)*dy - y
## Initial conditions
y(0) = 2
dy(0) = 0
## mu
mu = 1 ## nonstiff; 10 moderately stiff; 1000 stiff
", calcSens=TRUE)

    norm <- RxODE(sens, calcSens=FALSE)

    jac <- RxODE(sens, calcJac=TRUE)

    test_that("Jacobian and sensitivity specified.", {
        expect_false(norm$calcJac);
        expect_false(norm$calcSens);
        expect_true(sens$calcJac);
        expect_true(jac$calcJac);
        expect_false(jac$calcSens);
        expect_true(sens$calcJac);
        expect_true(sens$calcSens);
    })

    test_that("Unsupported derivatives will thow an error.", {
        expect_error(rxFromSymPy("E1(4)"))
    })
})

