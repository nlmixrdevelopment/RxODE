context("Conditional statements")
rxPermissive({

    ## else if is actually already supported...
    test_that("else if", {
        m <- RxODE({
            if (cnd <= 1){
                a = 1.0
            } else if (cnd <= 2){
                a = 2.0
            } else if (cnd <= 3){
                a = 3.
            } else {
                a = 100;
            }
            tmp = cnd
        })

        ## The prefered syntax is only if / else but it still works...
        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)

    })

    test_that("ifelse", {
        m <- RxODE({
            a = ifelse(cnd <= 1, 1.0, ifelse(cnd <= 2, 2, ifelse(cnd <= 3, 3, 100)))
            tmp = cnd
        })
        ## The prefered syntax is only if / else but it still works...
        tmp <- rxSolve(m, c(cnd=1), et(0.1))
        expect_equal(tmp$tmp, 1)
        expect_equal(tmp$a, 1)

        tmp <- rxSolve(m, c(cnd=2), et(0.1))
        expect_equal(tmp$tmp, 2)
        expect_equal(tmp$a, 2)

        tmp <- rxSolve(m, c(cnd=3), et(0.1))
        expect_equal(tmp$tmp, 3)
        expect_equal(tmp$a, 3)

        tmp <- rxSolve(m, c(cnd=4), et(0.1))
        expect_equal(tmp$tmp, 4)
        expect_equal(tmp$a, 100)
    })

    

    

    ##
})
