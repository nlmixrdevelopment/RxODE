library("RxODE")


context("Test addl->ss solving")

rxPermissive({

    rx <- RxODE({
        ## RxODE-style differential equations are supported
        KA <- 2
        Cl <- 3.5
        Vc <- 80
        d/dt(depot) ~ -KA*depot;
        d/dt(centr) ~ KA*depot-(Cl/Vc)*centr;
        ## Concentration is calculated
        cp = centr / Vc
    })

    et1 <- et(amt=50, until=800, ii=4) %>% et(0, 900, length.out=100)

    et1 <- etTrans(et1, rx)

    s1 <- rxSolve(rx, et1, addlSS=TRUE)

    s2 <- rxSolve(rx, et1, addlSS=FALSE)

    expect_equal(as.data.frame(s1), as.data.frame(s2))

    ## plot(microbenchmark::microbenchmark(rxSolve(rx, et1, addlSS=TRUE), rxSolve(rx, et1, addlSS=FALSE)), log="y")

    ## plot(s1)

})
