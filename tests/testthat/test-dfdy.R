require(RxODE);
context("Test Jacobian parsing")

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
})
