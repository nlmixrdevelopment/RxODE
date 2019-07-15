rxPermissive({
    context("Test symengine<->RxODE dsl")

    test_that("d/dt(x) parsing", {
        expect_equal(rxToSE(d/dt(matt)), "rx__d_dt_matt__")
        expect_equal(rxFromSE(rx__d_dt_matt__), "d/dt(matt)")
        expect_equal(rxToSE(d / dt( matt )), "rx__d_dt_matt__")
        expect_equal(rxFromSE("rx__d_dt_matt__"), "d/dt(matt)")
        tmp <- symengine::Symbol("rx__d_dt_matt__")
        expect_equal(rxFromSE(tmp), "d/dt(matt)")
    })

    test_that("df(x)/dy(x) parsing", {
        expect_equal(rxToSE(df(matt)/dy(ruth)), "rx__df_matt_dy_ruth__")
        expect_equal(rxFromSE(rx__df_matt_dy_ruth__), "df(matt)/dy(ruth)")
        expect_equal(rxToSE(df(matt) / dy( ruth )), "rx__df_matt_dy_ruth__")
        expect_equal(rxFromSE("rx__df_matt_dy_ruth__"), "df(matt)/dy(ruth)")
        tmp <- symengine::Symbol("rx__df_matt_dy_ruth__")
        expect_equal(rxFromSE(tmp), "df(matt)/dy(ruth)")
    })

    test_that("function and constant translation", {

        expect_equal(rxToSE(gammafn(a)), "gamma(a)")
        expect_error(rxToSE(gammafn(a, b)))
        expect_error(rxFromSE(gamma(a, b)))

        expect_equal(rxToSE(lgammafn(a)), "loggamma(a)")
        expect_equal(rxToSE(lgammafn(a)), "loggamma(a)")
        expect_equal(rxFromSE(loggamma(a)), "lgamma(a)")

        expect_equal(rxToSE(digamma(a)), "polygamma(0,a)")
        expect_equal(rxFromSE(polygamma(0, a)), "digamma(a)")

        expect_equal(rxToSE(trigamma(a)), "polygamma(1,a)")
        expect_equal(rxFromSE(polygamma(1, a)), "trigamma(a)")

        expect_equal(rxToSE(tetragamma(a)), "polygamma(2,a)")
        expect_equal(rxFromSE(polygamma(2,a)), "tetragamma(a)")

        expect_equal(rxToSE(pentagamma(a)), "polygamma(3,a)")
        expect_equal(rxFromSE(polygamma(3, a)), "pentagamma(a)")

        expect_equal(rxToSE(lbeta(a, b)), "log(beta(a,b))")
        expect_equal(rxFromSE(log(beta(a, b))), "lbeta(a,b)")

        expect_equal(rxToSE(lgamma1p(a)), "loggamma(a+1)")
        expect_equal(rxFromSE(loggamma(a + 1)), "lgamma1p(a)")
        expect_equal(rxFromSE(loggamma(1 + a)), "lgamma1p(a)")
        expect_equal(rxFromSE(loggamma(a + 1 + b)), "lgamma1p(a+b)")

        expect_equal(rxToSE(log1p(a)), "log(1+a)")
        expect_equal(rxFromSE(log(a + 1)), "log1p(a)")

        expect_equal(rxToSE(cospi(a)), "cos(pi*(a))")
        expect_equal(rxFromSE(cos(pi * a)), "cospi(a)")
        expect_equal(rxFromSE(cos(pi * (a))), "cospi(a)")

        expect_equal(rxFromSE(cos(b * pi * a)), "cospi(b*a)")

        expect_equal(rxToSE(sinpi(a)), "sin(pi*(a))")
        expect_equal(rxFromSE(sin(pi * a)), "sinpi(a)")

        expect_equal(rxToSE(tanpi(a)), "tan(pi*(a))")
        expect_equal(rxFromSE(tan(pi * a)), "tanpi(a)")
        expect_equal(rxFromSE("tan(pi/2)"), "tanpi(1/2)")

        expect_equal(rxToSE(log1pmx(a)), "(log(1+a)-(a))")
        expect_equal(rxToSE(expm1(a)), "(exp(a)-1)")
        expect_equal(rxToSE(pow(a, b)), "(a)^(b)")
        expect_error(rxToSE(pow(a, b,c)))
        expect_error(rxToSE(pow(a)))
        expect_equal(rxToSE(R_pow(a, b)), "(a)^(b)")
        expect_equal(rxToSE(Rx_pow_di(a, b)), "(a)^(b)")
        expect_equal(rxToSE(Rx_pow(a, b)), "(a)^(b)")
        expect_equal(rxToSE(R_pow_di(a, b)), "(a)^(b)")
        expect_equal(rxToSE(factorial(n)),"gamma(n+1)")

        expect_equal(rxToSE(beta(a, b)), "beta(a,b)")
        expect_error(rxToSE(beta(a)))

        expect_equal(rxToSE(choose(n, k)),"gamma(n+1)/(gamma(k+1)*gamma(n-(k)+1))")
        expect_equal(rxToSE(lchoose(n, k)), "(loggamma(n+1)-loggamma(k+1)-loggamma(n-(k)+1))")

        expect_equal(rxFromSE(log(-x + 1)), "log1p(-x)");
        expect_equal(rxFromSE(log(1 - x)), "log1p(-x)");

        expect_equal(rxToSE(log1pexp(x)), "log(1+exp(x))")
        expect_equal(rxFromSE(log(1 + exp(x))), "log1pexp(x)")

        ## expect_equal(rxFromSymPy(log((1 + x)^2)), "log1p(Rx_pow_di(x, 2)+2 * x)")
        ## expect_equal(rxFromSymPy(log((0.75 + x)^2)), "log1p(Rx_pow_di(x, 2)+1.5 * x-0.4375)")

        expect_equal(rxFromSE(E), "M_E")
        expect_equal(rxToSE(M_E), "E")

        expect_equal(rxToSE(exp(1)), "E")
        expect_equal(rxFromSE(exp(1)), "M_E")

        expect_equal(rxFromSE(log(2)), "M_LN2")
        expect_equal(rxFromSE(log(10)), "M_LN10")

        expect_equal(rxToSE("M_LN2"), "log(2)")
        expect_equal(rxFromSE("log(2)"), "M_LN2")

        expect_equal(rxToSE("M_LN10"), "log(10)")
        expect_equal(rxFromSE("log(10)"), "M_LN10")

        expect_equal(rxToSE("M_LN_SQRT_PI"), "log(sqrt(pi))")
        expect_equal(rxFromSE("log(sqrt(pi))"), "M_LN_SQRT_PI")
        expect_equal(rxFromSE("log(pi**0.5)"), "M_LN_SQRT_PI")
        expect_equal(rxFromSE("log(pi^0.5)"), "M_LN_SQRT_PI")
        expect_equal(rxFromSE("log(pi^(1/2))"), "M_LN_SQRT_PI")
        expect_equal(rxFromSE("log(pi**(1/2))"), "M_LN_SQRT_PI")

        expect_equal(rxToSE("M_LN_SQRT_PId2"), "log(sqrt(pi/2))")
        expect_equal(rxFromSE(log(sqrt(pi/2))), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log(sqrt(pi * 0.5))), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log(sqrt(0.5 * pi))), "M_LN_SQRT_PId2")

        expect_equal(rxFromSE(log((pi/2) ^ 0.5)), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log((pi * 0.5) ^ 0.5)), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log((0.5 * pi) ^ 0.5)), "M_LN_SQRT_PId2")

        expect_equal(rxFromSE(log((pi/2) ** 0.5)), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log((pi * 0.5) ** 0.5)), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log((0.5 * pi) ** 0.5)), "M_LN_SQRT_PId2")

        expect_equal(rxFromSE(log((pi/2) ** (1 / 2))), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log((pi * 0.5) ** (1 / 2))), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log((0.5 * pi) ** (1 / 2))), "M_LN_SQRT_PId2")

        expect_equal(rxFromSE(log((pi/2) ^ (1 / 2))), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log((pi * 0.5) ^ (1 / 2))), "M_LN_SQRT_PId2")
        expect_equal(rxFromSE(log((0.5 * pi) ^ (1 / 2))), "M_LN_SQRT_PId2")

        expect_equal(rxToSE("M_LN_SQRT_2PI"), "log(sqrt(2*pi))")
        expect_equal(rxFromSE(log((pi * 2) ^ 0.5)), "M_LN_SQRT_2PI")
        expect_equal(rxFromSE(log((2 * pi) ^ 0.5)), "M_LN_SQRT_2PI")

        expect_equal(rxFromSE(log((pi * 2) ** 0.5)), "M_LN_SQRT_2PI")
        expect_equal(rxFromSE(log((2 * pi) ** 0.5)), "M_LN_SQRT_2PI")

        expect_equal(rxFromSE(log((pi * 2) ** (1 / 2))), "M_LN_SQRT_2PI")
        expect_equal(rxFromSE(log((2 * pi) ** (1 / 2))), "M_LN_SQRT_2PI")

        expect_equal(rxFromSE(log((pi * 2) ^ (1 / 2))), "M_LN_SQRT_2PI")
        expect_equal(rxFromSE(log((2 * pi) ^ (1 / 2))), "M_LN_SQRT_2PI")


        expect_equal(rxToSE("M_SQRT_3"), "sqrt(3)")
        expect_equal(rxFromSE("sqrt(3)"), "M_SQRT_3")

        expect_equal(rxToSE("M_SQRT_2"), "sqrt(2)")
        expect_equal(rxFromSE("sqrt(2)"), "M_SQRT_2")

        expect_equal(rxToSE("M_SQRT_32"), "sqrt(32)")
        expect_equal(rxFromSE("sqrt(32)"), "M_SQRT_32")

        expect_equal(rxToSE("M_SQRT_PI"), "sqrt(pi)")
        expect_equal(rxFromSE("sqrt(pi)"), "M_SQRT_PI")

        expect_equal(rxToSE("M_SQRT_2dPI"), "sqrt(2/pi)")
        expect_equal(rxFromSE("sqrt(2/pi)"), "M_SQRT_2dPI")
        expect_equal(rxFromSE("(2/pi)^0.5"), "M_SQRT_2dPI")
        expect_equal(rxFromSE("(2/pi)^(1/2)"), "M_SQRT_2dPI")

        expect_equal(rxToSE("M_PI_2"), "pi/2")
        expect_equal(rxFromSE("pi*0.5"), "M_PI_2")
        expect_equal(rxFromSE("0.5*pi"), "M_PI_2")
        expect_equal(rxFromSE("pi/2"), "M_PI_2")

        expect_equal(rxToSE("M_PI_4"), "pi/4")
        expect_equal(rxFromSE("pi*0.25"), "M_PI_4")
        expect_equal(rxFromSE("0.25*pi"), "M_PI_4")
        expect_equal(rxFromSE("pi/4"), "M_PI_4")

        expect_equal(rxFromSE("1/pi"), "M_1_PI")
        expect_equal(rxToSE("M_1_PI"), "1/pi")

        expect_equal(rxFromSE("2/pi"), "M_2_PI")
        expect_equal(rxToSE("M_2_PI"), "2/pi")

        expect_equal(rxToSE("M_2_SQRTPI"), "2/sqrt(pi)")
        expect_equal(rxFromSE("2/sqrt(pi)"), "M_2_SQRTPI")
        expect_equal(rxFromSE("2/(pi^0.5)"), "M_2_SQRTPI")

        expect_equal(rxToSE(M_1_SQRT_2PI), "1/sqrt(2*pi)")
        expect_equal(rxFromSE("1/sqrt(2*pi)"), "M_1_SQRT_2PI")
        expect_equal(rxFromSE("1/sqrt(pi*2)"), "M_1_SQRT_2PI")
        expect_equal(rxFromSE("1/((pi*2)^(1/2))"), "M_1_SQRT_2PI")
        expect_equal(rxFromSE("1/((pi*2)^0.5)"), "M_1_SQRT_2PI")

        expect_equal(rxToSE(log10(a)), "log(a)/log(10)")
        expect_equal(rxToSE(log2(a)), "log(a)/log(2)")
        ## FIXME log10 log2? fromSE?


        ## expect_equal(rxFromSymPy("3 + 4*3+2+2*matt*pi"), "3 + 4 * 3 + 2 + matt * M_2PI")
        ## expect_equal(rxFromSymPy("3 + 4*3+2+pi*matt*2"), "3 + 4 * 3 + 2 + matt * M_2PI")
    })

    test_that("transit compartment translation.",{
        expect_equal(rxToSE(transit(n, mtt, bio)),
                  "exp(log((bio)*(podo))+log(n + 1)-log(mtt)+(n)*((log(n+1)-log(mtt))+log(t))-((n+1)/(mtt))*(t)-loggamma(1+n))")
        expect_equal(rxToSE(transit(n, mtt)),
                  "exp(log(podo)+(log(n+1)-log(mtt))+(n)*((log(n+1)-log(mtt))+ log(t))-((n + 1)/(mtt))*(t)-loggamma(1+n))")
    })

    test_that("unknown functions throw errors. rxToSE", {
        expect_error(rxToSE(matt(3)));
    })

    test_that("Theta/eta conversion rxToSE", {
        expect_equal(rxToSE(THETA[1]), "THETA_1_")
        expect_equal(rxToSE(ETA[1]), "ETA_1_")
        expect_equal(rxToSE(df(matt) / dy(THETA[1])), "rx__df_matt_dy_THETA_1___")
        expect_equal(rxToSE(df(matt) / dy(ETA[1])), "rx__df_matt_dy_ETA_1___")
        expect_error(rxToSE(THETA[0]))
        expect_error(rxToSE(ETA[0]))
        expect_error(rxToSE(THETA[0.5]))
        expect_error(rxToSE(ETA[0.5]))
        expect_error(rxToSE(THETA[a]))
        expect_error(rxToSE(ETA[a]))
        expect_error(rxToSE(THETA["b"]))
        expect_error(rxToSE(ETA["b"]))
    })

    test_that("Extra Derivative Table Test", {
        expect_equal(rxFromSE("Derivative(rxTBS(a, b, c), a)"), "rxTBSd(a,b,c)")
        expect_equal(rxFromSE("Derivative(rxTBSd(a, b, c), a)"), "rxTBSd2(a,b,c)")
        expect_error(rxFromSE("Derivative(rxTBSd2(a, b, c), a)", unknownDerivatives="error"))
        expect_error(rxFromSE("Derivative(rxTBS(a, b, c), d)", unknownDerivatives="error") )
        expect_equal(rxFromSE("(2*a + b)*Subs(Derivative(rxTBS(_xi_1, b, c), _xi_1), (_xi_1), (a*b + a^2))"),
                     "(2*a+b)*rxTBSd(a*b+a^2,b,c)")
    })

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
        expect_equal(rxToSymPy(Rx_pow_di(a, b)), "(a)**(b)")
        expect_equal(rxToSymPy(Rx_pow(a, b)), "(a)**(b)")
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
    context("Test DSL THETA/ETA")
    test_that("Theta/eta conversion", {
        expect_equal(rxToSymPy(THETA[1]), "THETA_1_")
        expect_equal(rxToSymPy(ETA[1]), "ETA_1_")
        expect_equal(rxToSymPy(df(matt) / dy(THETA[1])), "rx__df_matt_dy_THETA_1___")
        expect_equal(rxToSymPy(df(matt) / dy(ETA[1])), "rx__df_matt_dy_ETA_1___")
        expect_error(rxToSymPy(THETA[0]))
        expect_error(rxToSymPy(ETA[0]))
        expect_error(rxToSymPy(THETA[0.5]))
        expect_error(rxToSymPy(ETA[0.5]))
        expect_error(rxToSymPy(THETA[a]))
        expect_error(rxToSymPy(ETA[a]))
        expect_error(rxToSymPy(THETA["b"]))
        expect_error(rxToSymPy(ETA["b"]))
        expect_equal(rxFromSymPy(rx__df_matt_dy_THETA_1___), "df(matt)/dy(THETA[1])")
        expect_equal(rxFromSymPy(rx__df_matt_dy_ETA_1___), "df(matt)/dy(ETA[1])")
        expect_equal(rxFromSymPy(ETA_1_), "ETA[1]")
        expect_equal(rxFromSymPy(THETA_1_), "THETA[1]")
        expect_equal(rxFromSymPy(BY_THETA_1_), "BY_THETA_1_")
    })

    context("Test factor expansion by `rxSplitPlusQ'")
    test_that("rxSplitPlusQ", {
        expect_equal(rxSplitPlusQ(quote(a*exp(b+c)+d*log(e-f)-g*f)), c("a * exp(b + c)", "d * log(e - f)", "- g * f"))
        expect_equal(rxSplitPlusQ(quote(-a*exp(b+c)+d*log(e-f)-g*f)), c("-a * exp(b + c)", "d * log(e - f)", "- g * f"))
        expect_equal(rxSplitPlusQ(quote( + a*exp(b+c)+d*log(e-f)-g*f)), c("+a * exp(b + c)", "d * log(e - f)", "- g * f"))
        expect_equal(rxSplitPlusQ(quote( + a*exp(b+c))), "+a * exp(b + c)")
        expect_equal(rxSplitPlusQ(quote(center)), "center");
        expect_equal(rxSplitPlusQ(quote(0)), "0");
    })

    context("Choose numerically accurate functions")
    test_that("Promote to more numerically accurate functions", {
        expect_equal(rxFromSymPy(log(1 + x)), "log1p(x)");
        expect_equal(rxFromSymPy(log(x + 1)), "log1p(x)");
        expect_equal(rxFromSymPy(log(-x + 1)), "log1p(- x)");
        expect_equal(rxFromSymPy(log(1 - x)), "log1p(- x)");
        expect_equal(rxFromSymPy(log(1 + exp(x))), "log1pexp(x)")
        expect_equal(rxFromSymPy(log((1 + x)^2)), "log1p(Rx_pow_di(x, 2)+2 * x)")
        expect_equal(rxFromSymPy(log((0.75 + x)^2)), "log1p(Rx_pow_di(x, 2)+1.5 * x-0.4375)")
        expect_equal(rxFromSymPy(E), "M_E")
        expect_equal(rxToSymPy(M_E), "E")
        expect_equal(rxFromSymPy(log(2)), "M_LN2")
        expect_equal(rxFromSymPy(log(10)), "M_LN10")
        expect_equal(rxToSymPy("M_LN2"), "log(2)")
        expect_equal(rxToSymPy("M_LN10"), "log(10)")
        expect_equal(rxFromSymPy("log(sqrt(pi))"), "M_LN_SQRT_PI")
        expect_equal(rxFromSymPy("log(pi**0.5)"), "M_LN_SQRT_PI")
        expect_equal(rxToSymPy("M_LN_SQRT_PI"), "log(sqrt(pi))")
        expect_equal(rxFromSymPy("log(sqrt(pi/2))"), "M_LN_SQRT_PId2")
        expect_equal(rxToSymPy("M_LN_SQRT_PId2"), "log(sqrt(pi/2))")
        expect_equal(rxFromSymPy("log(sqrt(2*pi))"), "M_LN_SQRT_2PI")
        expect_equal(rxFromSymPy("log(sqrt(pi*2))"), "M_LN_SQRT_2PI")
        expect_equal(rxFromSymPy("log((pi*2)**0.5)"),"M_LN_SQRT_2PI")
        expect_equal(rxToSymPy("M_LN_SQRT_2PI"), "log(sqrt(2*pi))")
        expect_equal(rxFromSymPy("log((pi/2)^0.5)"), "M_LN_SQRT_PId2")
        expect_equal(rxFromSymPy("log((pi/2)^(0.5))"),"M_LN_SQRT_PId2")
        expect_equal(rxFromSymPy("sqrt(3)"), "M_SQRT_3")
        expect_equal(rxToSymPy("M_SQRT_3"), "sqrt(3)")
        expect_equal(rxFromSymPy("sqrt(2)"), "M_SQRT_2")
        expect_equal(rxToSymPy("M_SQRT_2"), "sqrt(2)")
        expect_equal(rxFromSymPy("sqrt(32)"), "M_SQRT_32")
        expect_equal(rxToSymPy("M_SQRT_32"), "sqrt(32)")
        expect_equal(rxFromSymPy("sqrt(pi)"), "M_SQRT_PI")
        expect_equal(rxToSymPy("M_SQRT_PI"), "sqrt(pi)")
        expect_equal(rxFromSymPy("sqrt(2/pi)"), "M_SQRT_2dPI")
        expect_equal(rxFromSymPy("2*pi*matt"), "M_2PI * matt")
        expect_equal(rxFromSymPy("pi*2*matt"), "M_2PI * matt")
        expect_equal(rxFromSymPy("3 + 4*3+2+2*matt*pi"), "3 + 4 * 3 + 2 + matt * M_2PI")
        expect_equal(rxFromSymPy("3 + 4*3+2+pi*matt*2"), "3 + 4 * 3 + 2 + matt * M_2PI")
        expect_equal(rxFromSymPy("pi/2"), "M_PI_2")
        expect_equal(rxFromSymPy("pi/4"), "M_PI_4")
        expect_equal(rxFromSymPy("1/pi"), "M_1_PI")
        expect_equal(rxFromSymPy("2/pi"), "M_2_PI")
        expect_equal(rxFromSymPy("2/sqrt(pi)"), "M_2_SQRTPI")
        expect_equal(rxFromSymPy("1/sqrt(2*pi)"), "M_1_SQRT_2PI")
        expect_equal(rxFromSymPy("1/sqrt(pi*2)"), "M_1_SQRT_2PI")
        ## Test conversions
        expect_equal(rxFromSymPy("log(matt)/log(10)"), "log10(matt)")
        expect_equal(rxFromSymPy("log(matt)/log(2)"), "log2(matt)")
        expect_equal(rxToSymPy("matt*M_LOG10E"), "matt * log10(E)")
        expect_equal(rxToSymPy("matt*M_LOG2E"), "matt/log(2)")
        expect_equal(rxFromSymPy("matt*log10(E)"), "matt * M_LOG10E")
        expect_equal(rxFromSymPy("matt/log(2)"), "matt * M_LOG2E")
        expect_equal(rxFromSymPy("matt/log(10)"), "matt * M_LOG10E")
        expect_equal(rxFromSymPy("matt/log(2)"), "matt * M_LOG2E")
        expect_equal(rxFromSymPy("1/log(2)"), "M_LOG2E")
        expect_equal(rxFromSymPy("1/log(10)"), "M_LOG10E")
        expect_equal(rxToSymPy("log2(exp(1))"), "1/log(2)")
        expect_equal(rxToSymPy("log2(M_E)"), "1/log(2)")
        expect_equal(rxToSymPy("log10(exp(1))"), "1/log(10)")
        expect_equal(rxToSymPy("log10(M_E)"), "1/log(10)")
        expect_equal(rxFromSymPy("sin(pi)"), "sinpi(1)")
        expect_equal(rxFromSymPy("sin(2*pi)"), "sinpi(2)")
        expect_equal(rxFromSymPy("cos(2*pi)"), "cospi(2)")
        expect_equal(rxFromSymPy("tan(2*pi)"), "tanpi(2)")
        expect_equal(rxFromSymPy("tan(pi/2)"), "tanpi(1/2)")
        expect_equal(rxFromSymPy("tan((pi/2)^2+pi)"), "tanpi(M_PI_2 + 1)")
        expect_equal(rxFromSymPy("tan((pi/2)^2+pi+1)"), "tan(Rx_pow_di(M_PI_2, 2) + M_PI + 1)")
    })

    context("Test Error DSLs")

    pk <- function(){
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V = exp(THETA[3] + ETA[2])
        return(1)
    }

    pred <- function(){
        if (cmt == 1){
            return(centr);
        } else if (cmt == 2){
            return(depot);
        } else {
            return(perip);
        }
    }


    pred.for <- function(){
        for (i in 1:10){
            pred = 1;
        }
    }

    err <- function(){
        add(0.2) + prop(0.3)
    }

    test_that("PK/Pred/Error function parsing.", {
        expect_error(rxParsePk(pk));
        expect_error(rxParsePred(pred));
        expect_error(rxParsePred(pred.for))
        expect_error(rxParsePk(err));
        expect_error(rxParsePred(err));
    })


    pk2 <- function(){
        KA = exp(THETA[1])
        CL = exp(THETA[2] + ETA[1])
        V = exp(THETA[3] + ETA[2])
    }

}, on.validate=TRUE)

