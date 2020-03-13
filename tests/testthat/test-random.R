rxPermissive({

    context("normal random variables")

    test_that("rnorm", {

        set.seed(1024)

        rx <- RxODE({
            x1 <- rnorm()
            x2 <- rxnorm(a)
            x3 <- rnorm(b, c)
            d/dt(x0) = 0
        })

        ev <- et(1, id=1:30000)

        f <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=2)

        expect_equal(mean(f$x1), 0, tol=1e-2)
        expect_equal(sd(f$x1), 1, tol=1e-2)

        expect_equal(mean(f$x2), 3, tol=1e-2)
        expect_equal(sd(f$x1), 1, tol=1e-2)


        expect_equal(mean(f$x3), 5, tol=1e-2)
        expect_equal(sd(f$x3), 2, tol=1e-2)

        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        expect_equal(mean(f2$x1), 0, tol=1e-2)
        expect_equal(sd(f2$x1), 1, tol=1e-2)

        expect_equal(mean(f2$x2), 3, tol=1e-2)
        expect_equal(sd(f2$x1), 1, tol=1e-2)

        expect_equal(mean(f2$x3), 5, tol=1e-2)
        expect_equal(sd(f2$x3), 2, tol=1e-2)

        expect_error(RxODE({
            x4 <- rnorm(a, b, c, d)
        }))

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))


        x <- rxnorm(n=1e5)
        expect_equal(mean(x), 0, tol=0.01)

        expect_equal(sd(x), 1, tol=0.01)


    })

    context("rnormV")
    test_that("rnormV", {

        set.seed(1024)

        rx <- RxODE({
            x1 <- rnormV()
            x2 <- rxnormV(a)
            x3 <- rnormV(b, c)
            d/dt(x0) = 0
        })

        expect_error(RxODE({
            x4 <- rnormV(a, b, c, d)
        }))

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))


        x <- rxnorm(n=1e4)

        expect_equal(mean(x), 0, tol=0.1)



    })

    context("binomial tests")
    test_that("rbinom", {

        rx <- RxODE({
            x1 <- rbinom(4, 0.5)
            x2 <- rxbinom(10, 0.75)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        expect_equal(max(f$x1), 4)
        expect_equal(min(f$x1), 0)
        expect_true(all(round(f$x1) == f$x1))
        expect_equal(mean(f$x1), 4 * 0.5, tol=1e-2)
        expect_equal(sd(f$x1), sqrt(4 * 0.5 * 0.5), tol=1e-2)

        expect_equal(max(f$x2), 10)
        expect_true(min(f$x2) > 0)
        expect_true(all(round(f$x2) == f$x2))

        expect_equal(mean(f$x2), 10 * 0.75, tol=1e-2)
        expect_equal(sd(f$x2), sqrt(10 * 0.75 * 0.25), tol=1e-2)

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x1 <- rbinom()
        }))

        expect_error(RxODE({
            x1 <- rbinom(a)
        }))

        expect_error(RxODE({
            x1 <- rbinom(a, b, c)
        }))

    })

    context("Cauchy random numbers")

    test_that("rcauchy", {

        set.seed(1024)

        rx <- RxODE({
            x1 <- rcauchy()
            x2 <- rxcauchy(a)
            x3 <- rcauchy(b, c)
            d/dt(x0) = 0
        })

        ev <- et(1, id=1:100)

        f <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=2)


        ## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x4 <- rcauchy(a, b, c, d)
        }))

    })

    context("rchisq tests")
    test_that("rchisq", {


        rx <- RxODE({
            x1 <- rchisq(15)
            x2 <- rxchisq(20)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        expect_equal(mean(f$x1), 15, tol=0.1)
        expect_equal(sd(f$x1), sqrt(2 * 15), tol=0.1)

        expect_equal(mean(f$x2), 20, tol=0.1)
        expect_equal(sd(f$x2), sqrt(2 * 20), tol=0.1)

        ## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))


        expect_error(RxODE({
            x1 <- rchisq()
        }))

        expect_error(RxODE({
            x1 <- rchisq(a, b)
        }))
    })

    context("rexp tests")

    test_that("rexp tests", {

        rx <- RxODE({
            x1 <- rexp(0.5)
            x2 <- rxexp()
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        expect_equal(mean(f$x1), 2, tol=0.1)
        expect_equal(sd(f$x1), sqrt(1 / (0.5 * 0.5)), tol=0.1)

        expect_equal(mean(f$x2), 1, tol=0.1)
        expect_equal(sd(f$x2), 1, tol=0.1)

        ## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x1 <- rexp(a, b)
        }))

    })

    context("rf tests")

    test_that("rf tests", {

        rx <- RxODE({
            x1 <- rf(10, 20)
            x2 <- rxf(30, 40)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        sf <- function(d1, d2){
            sqrt((2 * d2 ^ 2 * (d1 + d2 - 2)) / (d1 * (d2 - 2) ^ 2 * (d2 - 4)))
        }

        mf <- function(d2){
            return(d2 / (d2 - 2))
        }

        expect_equal(mean(f$x1), mf(20), tol=0.01)
        expect_equal(sd(f$x1), sf(10, 20), tol=0.01)

        expect_equal(mean(f$x2), mf(40), tol=0.01)
        expect_equal(sd(f$x2), sf(30, 40), tol=0.01)

        ## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x1 <- rf(a, b, c)
        }))

        expect_error(RxODE({
            x1 <- rf(a)
        }))

        expect_error(RxODE({
            x1 <- rf()
        }))

    })

    context("rgamma tests")

    test_that("rgamma tests", {

        rx <- RxODE({
            x1 <- rgamma(9, 0.5)
            x2 <- rxgamma(7.5)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        sgamma <- function(k, theta=1){
            sqrt(k / (theta ^ 2))
        }

        expect_equal(sd(f$x1), sgamma(9, 0.5), tol=0.01)

        expect_equal(sd(f$x2), sgamma(7.5), tol=0.01)

        ## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x1 <- rgamma(a, b, c)
        }))

        expect_error(RxODE({
            x1 <- rgamma()
        }))

    })

    context("rbeta tests")

    test_that("rbeta tests", {

        rx <- RxODE({
            x1 <- rbeta(2, 5)
            x2 <- rxbeta(2, 2)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)


        mbeta <- function(a, b){
            return(a / (a + b))
        }
        sbeta <- function(a, b){
            sqrt(a * b / ((a + b) ^ 2 * (a + b + 1)))
        }

        expect_equal(mean(f$x1), mbeta(2, 5), tol=0.01)
        expect_equal(sd(f$x1), sbeta(2, 5), tol=0.01)

        expect_equal(mean(f$x2), mbeta(2, 2), tol=0.01)
        expect_equal(sd(f$x2), sbeta(2, 2), tol=0.01)

        ## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x1 <- rbeta(a, b, c)
        }))

        expect_error(RxODE({
            x1 <- rbeta(a)
        }))

        expect_error(RxODE({
            x1 <- rbeta()
        }))

    })

    context("rgeom tests")

    test_that("rgeom tests", {

        rx <- RxODE({
            x1 <- rgeom(0.5)
            x2 <- rxgeom(0.1)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        expect_equal(median(f$x1), -ceiling(1/log2(1-0.5)))
        expect_equal(median(f$x2), -ceiling(1/log2(1-0.1)))

        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x1 <- rgeom()
        }))

        expect_error(RxODE({
            x1 <- rgeom(a, b)
        }))

    })

    context("rpois tests")

    test_that("rpois", {

        rx <- RxODE({
            x1 <- rpois(1)
            x2 <- rxpois(2)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        expect_equal(mean(f$x1), 1, tol=0.01)
        expect_equal(sd(f$x1), 1, tol=0.01)

        expect_equal(mean(f$x2), 2, tol=0.01)
        expect_equal(sd(f$x2), sqrt(2), tol=0.01)
## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x1 <- rpois()
        }))

        expect_error(RxODE({
            x1 <- rxpois(a, b)
        }))

    })

    context("rt tests")

    test_that("rt", {

        rx <- RxODE({
            x1 <- rt(15)
            x2 <- rxt(20)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        expect_equal(mean(f$x1), 0, tol=0.1)
        expect_equal(sd(f$x1), sqrt(15 / (15 - 2)), tol=0.1)

        expect_equal(mean(f$x2), 0, tol=0.1)
        expect_equal(sd(f$x2), sqrt(20 / (20 - 2)), tol=0.1)

        ## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))


        expect_error(RxODE({
            x1 <- rt()
        }))

        expect_error(RxODE({
            x1 <- rt(a, b)
        }))

    })

    context("runif tests");

    test_that("runif", {

        set.seed(1024)

        rx <- RxODE({
            x1 <- runif()
            x2 <- rxunif(a)
            x3 <- runif(b, c)
            d/dt(x0) = 0
        })

        ev <- et(1, id=1:30000)

        f <- rxSolve(rx, ev, c(a=0.5, b=0.25, c=0.75), cores=2)

        expect_equal(mean(f$x1), 0.5, tol=1e-2)
        expect_equal(sd(f$x1), sqrt(1 / 12), tol=1e-2)

        expect_equal(mean(f$x2), 0.5 * (0.5 + 1), tol=1e-2)
        expect_equal(sd(f$x2), sqrt((1 - 0.5) ^ 2 / 12), tol=1e-2)

        expect_equal(mean(f$x3), 0.5 * (0.25 + 0.75), tol=1e-2)
        expect_equal(sd(f$x3), sqrt((0.75 - 0.25) ^ 2 / 12), tol=1e-2)

        f2 <- rxSolve(rx, ev, c(a=0.5, b=0.25, c=0.75), cores=1)

        expect_equal(mean(f2$x1), 0.5, tol=1e-2)
        expect_equal(sd(f2$x1), sqrt(1 / 12), tol=1e-2)

        expect_equal(mean(f2$x2), 0.5 * (0.5 + 1), tol=1e-2)
        expect_equal(sd(f2$x2), sqrt((1 - 0.5) ^ 2 / 12), tol=1e-2)

        expect_equal(mean(f2$x3), 0.5 * (0.25 + 0.75), tol=1e-2)
        expect_equal(sd(f2$x3), sqrt((0.75 - 0.25) ^ 2 / 12), tol=1e-2)

        expect_error(RxODE({
            x4 <- runif(a, b, c, d)
        }))

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, c(a=3, b=5, c=2), cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

    })


    context("rweibull tests")

    test_that("rweibull tests", {

        rx <- RxODE({
            x1 <- rweibull(9, 0.5)
            x2 <- rxweibull(7.5)
        })

        ev <- et(1, id=1:30000)

        set.seed(1024)
        f <- rxSolve(rx, ev, cores=2)

        mweibull <- function(shape, scale=1){
            lambda <- scale;  k <- shape
            lambda * gamma(1 + 1 / k)
        }

        sweibull <- function(shape, scale=1){
            lambda <- scale;  k <- shape
            sqrt(lambda ^ 2 * (gamma(1 + 2 / k)
            -(gamma(1 + 1 / k)) ^ 2))
        }

        expect_equal(mean(f$x1), mweibull(9, 0.5), tol=0.01)
        expect_equal(sd(f$x1), sweibull(9, 0.5), tol=0.01)

        expect_equal(mean(f$x2), mweibull(7.5), tol=0.01)
        expect_equal(sd(f$x2), sweibull(7.5), tol=0.01)

        ## Seed tests

        ## Make sure seeds are reproducible
        ev <- et(1, id=1:10)

        set.seed(1)
        f <- rxSolve(rx, ev, cores=1)

        set.seed(1)
        f2 <- rxSolve(rx, ev, cores=1)
        expect_equal(as.data.frame(f), as.data.frame(f2))

        ## Make sure different seed value gives different result
        set.seed(2)
        f2 <- rxSolve(rx, ev, cores=1)

        expect_false(isTRUE(all.equal(as.data.frame(f), as.data.frame(f2))))

        expect_error(RxODE({
            x1 <- rweibull(a, b, c)
        }))

        expect_error(RxODE({
            x1 <- rweibull()
        }))

    })

}, test="norm")
