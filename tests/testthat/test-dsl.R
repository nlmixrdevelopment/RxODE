context("Test DSL rxToSymPy")

test_that("d/dt(x) parsing", {
    expect_equal(rxToSymPy(d/dt(matt)), "rx__d_dt_matt__")
    expect_equal(rxToSymPy(d / dt( matt )), "rx__d_dt_matt__")
})

test_that("df(x)/dy(x) parsing", {
    expect_equal(rxToSymPy(df(matt)/dy(ruth)), "rx__df_matt_dy_ruth__")
    expect_equal(rxToSymPy(df(matt) / dy( ruth )), "rx__df_matt_dy_ruth__")
})

test_that("bessel functions", {
    expect_equal(rxToSymPy(bessel_i(x, nu, 1)), "besseli(nu,x)")
    expect_equal(rxToSymPy(bessel_i(x, nu, 2)), "exp(-x) * besseli(nu,x)")
    expect_error(rxToSymPy(bessel_i(x, nu, 1.5)))
    expect_error(rxToSymPy(bessel_i(x, nu)))
    expect_equal(rxToSymPy(bessel_j(x, nu)), "besselj(nu,x)")
    expect_error(rxToSymPy(bessel_j(x, nu, 1)))
    expect_equal(rxToSymPy(bessel_k(x, nu, 1)), "besselk(nu,x)")
    expect_equal(rxToSymPy(bessel_k(x, nu, 2)), "exp(x) * besselk(nu,x)")
    expect_error(rxToSymPy(bessel_k(x, nu, 1.5)))
    expect_error(rxToSymPy(bessel_k(x, nu)))
    expect_equal(rxToSymPy(bessel_y(x, nu)), "bessely(nu,x)")
    expect_error(rxToSymPy(bessel_y(x, nu, 1)))
})

test_that("logspace add/subtract", {
    expect_equal(rxToSymPy(logspace_add(a, 1 + 3)), "(a) + (1 + 3)")
    expect_equal(rxToSymPy(logspace_sub(a, 1 + 3)), "(a) - (1 + 3)")
})

test_that("R math functions translations", {
    expect_equal(rxToSymPy(gammafn(a)), "gamma(a)")
    expect_equal(rxToSymPy(lgammafn(a)), "log(gamma(a))")
    expect_equal(rxToSymPy(lgammafn(a)), "log(gamma(a))")
    expect_equal(rxToSymPy(tetragamma(a)), "psigamma(a, 2)")
    expect_equal(rxToSymPy(pentagamma(a)), "psigamma(a, 3)")
    expect_equal(rxToSymPy(lbeta(a)), "log(beta(a))")
    expect_equal(rxToSymPy(lgamma1p(a)), "log(gamma((a)+1))")
    expect_equal(rxToSymPy(cospi(a)), "cos(pi * (a))")
    expect_equal(rxToSymPy(sinpi(a)), "sin(pi * (a))")
    expect_equal(rxToSymPy(tanpi(a)), "tan(pi * (a))")
    expect_equal(rxToSymPy(pow(a, b)), "(a)**(b)")
    expect_equal(rxToSymPy(R_pow(a, b)), "(a)**(b)")
    expect_equal(rxToSymPy(R_pow_di(a, b)), "(a)**(b)")
    expect_equal(rxToSymPy(log1p(a)), "log(1 + (a))")
    expect_equal(rxToSymPy(log1pmx(a)), "(log(1 + (a))-(a))")
    expect_equal(rxToSymPy(expm1(a)), "(exp(a)-1)")
    expect_equal(rxToSymPy(psigamma(z, n)), "polygamma(n, z)")
    expect_equal(rxToSymPy(choose(n, k)), "(factorial(n)/(factorial(k)*factorial((n)-(k))))")
    expect_equal(rxToSymPy(lchoose(n, k)), "(log(gamma((n)+1))-log(gamma((k)+1))-log(gamma((n)-(k)+1)))")
})

for (fn in RxODE:::sympy.equiv.f){
    test_that(sprintf("Equiv syntax '%s'", fn), {
        eval(parse(text=sprintf("expect_equal(rxToSymPy(%s(a)),\"%s(a)\")", fn, fn)))
    })
}

test_that("transit compartment translation.",{
    test_that(rxToSymPy(transit(n, mtt, bio)), "exp(log((bio) * (podo)) + (log((n) + 1) - log(mtt)) + (n) * ((log((n) + 1) - log(mtt)) + log(t)) - ((n + 1) / (mtt)) * (t) - lgamma(1 + (n)))")
    test_that(rxToSymPy(transit(n, mtt)), "exp(log((1) * (podo)) + (log((n) + 1) - log(mtt)) + (n) * ((log((n) + 1) - log(mtt)) + log(t)) - ((n + 1) / (mtt)) * (t) - lgamma(1 + (n)))")
})


tmp <- function(df = "matt"){rxToSymPy(sprintf("d/dt(%s)", df))}

tmp2 <- function(df = "matt"){df2 <- sprintf("d/dt(%s)", df);rxToSymPy(df2)}

x <- "matt";

test_that("sprintf non-standard evaluation correclty parses.", {
    expect_equal(tmp(), rxToSymPy(d / dt(matt)))
    expect_equal(tmp2(), rxToSymPy(d / dt(matt)))
    expect_equal(rxToSymPy(d/dt(ab)+sprintf("d/dt(%s)",x)),
                 rxToSymPy(d / dt(ab) + d / dt(matt)))
})

test_that("unknown functions throw errors.", {
    expect_error(rxToSymPy(matt(3)));
})

test_that("time prefers t notation", {
    expect_equal(rxToSymPy(time), "t");
})


context("Test DSL rxFromSymPy")

for (fn in RxODE:::sympy.equiv.f){
    test_that(sprintf("Equiv syntax '%s'", fn), {
        eval(parse(text=sprintf("expect_equal(rxFromSymPy(%s(a)),\"%s(a)\")", fn, fn)))
    })
}

test_that("d/dt(x) parsing", {
    expect_equal("d/dt(matt)", rxFromSymPy("rx__d_dt_matt__"))
    expect_equal("d/dt(matt)", rxFromSymPy(rx__d_dt_matt__))
})

test_that("df(x)/dy(x) parsing", {
    expect_equal("df(matt)/dy(ruth)", rxFromSymPy("rx__df_matt_dy_ruth__"))
    expect_equal("df(matt)/dy(ruth)", rxFromSymPy(rx__df_matt_dy_ruth__))
})

test_that("bessel functions & polygamma", {
    expect_equal(rxFromSymPy(besseli(x, nu)), "bessel_i(nu, x, 1.0)")
    expect_equal(rxFromSymPy(besselj(x, nu)), "bessel_j(nu, x)")
    expect_equal(rxFromSymPy(besselk(x, nu)), "bessel_k(nu, x, 1.0)")
    expect_equal(rxFromSymPy(bessely(x, nu)), "bessel_y(nu, x)")
    expect_equal(rxFromSymPy(polygamma(n, z)), "psigamma(z, n)")
})

test_that("unknown function raises an error.", {
    expect_error(rxFromSymPy(matt(3)))
})




test_that("Named component parsing", {
    tmp <- c("name"="d/dt(matt)");
    expect_equal(rxToSymPy(tmp), "rx__d_dt_matt__");
    tmp <- c("name"="rx__d_dt_matt__");
    expect_equal(rxFromSymPy(tmp), "d/dt(matt)");
})
