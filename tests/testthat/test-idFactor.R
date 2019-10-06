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
        expect_false(inherits(tmp$id, "factor"))

    })

    test_that("Test giving IDs to data-frames", {

        source("theoSd.R")
        d <- theoSd;

        mod  <- RxODE({
            ka <- exp(tka)
            cl <- exp(tcl)
            v <- exp(tv)
            cp <- linCmt()
        })

        set.seed(100)
        parData <- data.frame(id=1:12, tka=1 + rnorm(12, sd=0.01),
                              tcl=2 + rnorm(12, sd=0.01),
                              tv=3 + rnorm(12, sd=0.01))

        parData2 <- parData[order(-parData$id), ]

        tmp1 <- rxSolve(mod, d, parData)

        tmp2 <- rxSolve(mod, d, parData2)

        expect_equal(data.frame(tmp1),data.frame(tmp2))
        expect_equal(data.frame(tmp1$params),data.frame(tmp2$params))

        expect_equal(levels(tmp1$id), levels(tmp1$params$id))
        expect_equal(levels(tmp2$id), levels(tmp2$params$id))

        ## Now drop an id from the dataset
        d <- d[d$ID != 10, ];

        tmp1 <- expect_warning(rxSolve(mod, d, parData))

        tmp2 <- expect_warning(rxSolve(mod, d, parData2))

        ## Now add an id to the dataset

        source("theoSd.R")
        d <- theoSd;

        parData <- data.frame(id=1:13, tka=1 + rnorm(13, sd=0.01),
                              tcl=2 + rnorm(13, sd=0.01),
                              tv=3 + rnorm(13, sd=0.01))

        parData2 <- parData[order(-parData$id), ]

        tmp1 <- expect_warning(rxSolve(mod, d, parData))

        tmp2 <- expect_warning(rxSolve(mod, d, parData2))

        expect_equal(data.frame(tmp1), data.frame(tmp2))
        expect_equal(data.frame(tmp1$params),data.frame(tmp2$params))

        expect_equal(levels(tmp1$id), levels(tmp1$params$id))
        expect_equal(levels(tmp2$id), levels(tmp2$params$id))

        ## Now try letters
        source("theoSd.R")
        d <- theoSd;
        d$ID <- letters[d$ID]

        parData <- data.frame(id=1:13, tka=1 + rnorm(13, sd=0.01),
                              tcl=2 + rnorm(13, sd=0.01),
                              tv=3 + rnorm(13, sd=0.01))

        parData2 <- parData[order(-parData$id), ]

        expect_error(rxSolve(mod, d, parData))

        expect_error(rxSolve(mod, d, parData2))

        parData$id <- letters[parData$id]
        parData2$id <- letters[parData2$id]

        tmp1 <- expect_warning(rxSolve(mod, d, parData))

        tmp2 <- expect_warning(rxSolve(mod, d, parData2))

        expect_equal(data.frame(tmp1), data.frame(tmp2))
        expect_equal(data.frame(tmp1$params),data.frame(tmp2$params))

        expect_equal(levels(tmp1$id), levels(tmp1$params$id))
        expect_equal(levels(tmp2$id), levels(tmp2$params$id))

        ## Last case, though rare

        parData <- data.frame(id=1:12, tka=1 + rnorm(12, sd=0.01),
                              tcl=2 + rnorm(12, sd=0.01),
                              tv=3 + rnorm(12, sd=0.01))

        parData$ID <- letters[parData$id]

        tmp1 <- expect_warning(rxSolve(mod, d, parData))

        parData2 <- parData[order(-parData$id), ];

        tmp2 <- expect_warning(rxSolve(mod, d, parData2))

        expect_false(all(tmp1$params$tka == tmp2$params$tka))

    })



})
