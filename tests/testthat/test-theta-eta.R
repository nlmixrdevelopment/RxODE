context("Test RxODE THETA/ETA support")
library(digest)
rxPermissive({

    rigid <- RxODE({
        y1(0)    = 1
        y2(0)    = 0
        y3(0)    = 0.9
        a1       = theta[1]
        a2       = theta[2]
        a3       = theta[3] + eta[1]
        d/dt(y1) = a1*y2*y3
        d/dt(y2) = a2*y1*y3
        d/dt(y3) = a3*y1*y2
    })

    et <- eventTable();
    et$add.sampling(seq(0,20,by=0.01))

    out <- solve(rigid,et, theta=c(-2, 1.25, -0.5), eta=c(0))

    test_that("Test rigid body example",{
        expect_equal(digest(round(as.data.frame(out),3)),
                     "14ac439006ceb32455ba00c7817f9163")
    })

    out <- solve(rigid,et, theta=c(-2, 1.25, -0.5), eta=c(1))

    test_that("Test rigid body example",{
        expect_equal(digest(round(as.data.frame(out),3)),
                     "c7cffaa650a47e2b28e4cba99c603dde")
    })

    theta <- RxODE({
        a = cos(theta)
    })

    eta <- RxODE({
        a = cos(eta)
    })

    test_that("theta/eta only parsing works", {
        expect_equal(class(theta), "RxODE")
        expect_equal(class(eta), "RxODE")
    })


}, silent=TRUE)
