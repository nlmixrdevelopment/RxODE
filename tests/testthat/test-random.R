rxPermissive({

    context("normal random variables")
    test_that("rnorm", {

        set.seed(1024)

        rx <- RxODE({
            x1 <- rnorm()
            x2 <- rnorm(a)
            x3 <- rnorm(b, c)
            d/dt(x0) = 0
        })

        ev <- et(1, id=1:20000)

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


    })

})
