rxPermissive({
    context("Testing solving with ID(s) in the dataset")

    test_that("simple solving with ID(s) in the dataset", {

        source("theoSd.R")
        d <- theoSd[,names(theoSd) != "EVID"];
        d <- d[d$ID != 10, ];

        mod  <- RxODE({
            tka = 1
            tcl <- 2
            tv <- 3
            ka <- exp(tka)
            cl <- exp(tcl)
            v <- exp(tv)
            cp <- linCmt()
        })

        tmp <- rxSolve(mod, d)
        expect_true(inherits(tmp$id, "factor"))
        ## Notice 10 is missing.
        expect_equal(levels(tmp$id), c("1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "12"))

        tmp <- rxSolve(mod, d, idFactor=FALSE)
        expect_false(inherits(tmp$id, "factor"))

        d2 <- d
        d2$ID <- factor(d2$ID, c(1:10, 12), letters[c(1:10, 12)])

        tmp <- rxSolve(mod, d2)

        tmp <- rxSolve(mod, d, idFactor=FALSE)
        expect_true(inherits(tmp$id, "factor"))
        expect_equal(levels(tmp$id), letters[c(1:10, 12)])

    })


})
