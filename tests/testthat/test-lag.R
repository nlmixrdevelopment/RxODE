rxPermissive({

    et <- et(1:10)
    et$b <- 1:10

    context("lag() tests")
    test_that("lag()", {

        expect_error(RxODE({
            a = lag()
        }))

        expect_error(RxODE({
            a = lag(b, c)
        }))

        m1 <- RxODE({
            c = b + 2
            a = lag(b, 3)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(NA, NA, NA, 1, 2, 3, 4, 5, 6, 7))

        m1 <- RxODE({
            a = lag(b, 3)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(NA, NA, NA, 1, 2, 3, 4, 5, 6, 7))

        m1 <- RxODE({
            a = lag(b, -1)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(2, 3, 4, 5, 6, 7, 8, 9, 10, NA))

        m1 <- RxODE({
            a = lag(b, 0)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, 1:10)

        m1 <- RxODE({
            a = lag(b)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))

        m1 <- RxODE({
            c = b + 3
            a = lag(b)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))

    })

    context("lead()")

    test_that("lead()", {

        expect_error(RxODE({
            a = lead()
        }))

        expect_error(RxODE({
            a = lead(b, c)
        }))

        m1 <- RxODE({
            c = b + 2
            a = lead(b, 3)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(4, 5, 6, 7, 8, 9, 10, NA, NA, NA))

        m1 <- RxODE({
            a = lead(b, 3)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(4, 5, 6, 7, 8, 9, 10, NA, NA, NA))

        m1 <- RxODE({
            a = lead(b, -1)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(NA, 1, 2, 3, 4, 5, 6, 7, 8, 9))

        m1 <- RxODE({
            a = lead(b, 0)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, 1:10)


        m1 <- RxODE({
            a = lead(b)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(2:10, NA))

        m1 <- RxODE({
            c = b + 3
            a = lead(b)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_equal(x1$a, c(2:10, NA))

    })


    context("first()")

    test_that("first()", {

        expect_error(RxODE({
            a = first()
        }))

        expect_error(RxODE({
            a = first(b, 1)
        }))

        expect_error(RxODE({
            a = first(b, 1, 2)
        }))

        m1 <- RxODE({
            a = first(b)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_true(all(x1$a == 1))

        m1 <- RxODE({
            c = b + 3
            a = first(b)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_true(all(x1$a == 1))

    })

    context("last()")

    test_that("last()", {

        expect_error(RxODE({
            a = last()
        }))

        expect_error(RxODE({
            a = last(b, 1)
        }))

        expect_error(RxODE({
            a = last(b, 1, 2)
        }))

        m1 <- RxODE({
            a = last(b)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_true(all(x1$a == 10))

        m1 <- RxODE({
            c = b + 3
            a = last(b)
        })

        expect_true(inherits(m1, "RxODE"))

        x1 <- m1 %>% rxSolve(et)

        expect_true(all(x1$a == 10))

    })

})
