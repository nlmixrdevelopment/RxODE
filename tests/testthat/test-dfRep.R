rxPermissive({
    context("tests the internal df repetition routines")
    expect_equal(RxODE:::.vecDf(c(a=1,b=1,c=3),3),
                 structure(list(a = c(1, 1, 1),
                                b = c(1, 1, 1),
                                c = c(3, 3, 3)),
                           row.names = c(NA, -3L),
                           class = "data.frame"))
    expect_error(RxODE:::.vecDf(c(a=1,b=1,c=3),0))
}, test="cran")
